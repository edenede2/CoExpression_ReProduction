library(WGCNA)
options(stringsAsFactors = FALSE)
suppressMessages(WGCNA::allowWGCNAThreads())  # multi-thread WGCNA when available
.extract_donor <- function(x) sub("^([^-]+-[^-]+).*", "\\1", x)

# Average duplicate samples per donor (rows) within a tissue
.aggregate_by_donor <- function(mat) {
    d <- .extract_donor(rownames(mat))
    split_idx <- split(seq_len(nrow(mat)), d)
    res <- do.call(rbind, lapply(split_idx, function(ix) {
        colMeans(mat[ix, , drop = FALSE], na.rm = TRUE)
    }))
    rownames(res) <- names(split_idx)
    storage.mode(res) <- "double"
    res
}

.make_CT_map <- function(tissues, beta) {
  if (length(tissues) < 2) return(setNames(numeric(0), character(0)))
  pairs <- t(combn(tissues, 2))
  keys  <- paste(pairs[,1], pairs[,2], sep = "||")
  setNames(rep.int(beta, length(keys)), keys)
}

checkScaleFree_logbin <- function(k, nBreaks = 10, removeFirst = FALSE) {
  kpos <- k[k > 0]                      
  br <- unique(10^seq(log10(min(kpos)), log10(max(kpos)), length.out = nBreaks + 1))
  h  <- hist(kpos, breaks = br, plot = FALSE, right = TRUE)

  dk   <- h$mids                        
  p.dk <- h$counts / sum(h$counts)

  log.dk  <- log10(dk)
  if (removeFirst) { log.dk <- log.dk[-1]; p.dk <- p.dk[-1] }
  log.p   <- log10(p.dk + 1e-9)

  lm1 <- lm(log.p ~ log.dk)
  lm2 <- lm(log.p ~ log.dk + I(10^log.dk)) 

  data.frame(
    Rsquared.SFT = summary(lm1)$r.squared,
    slope.SFT    = coef(lm1)[2],
    truncatedExponentialAdjRsquared = summary(lm2)$adj.r.squared
  )
}

checkScaleFree <- function (k, nBreaks = 10, removeFirst = FALSE) 
{
  discretized.k = cut(k, nBreaks)
  dk = tapply(k, discretized.k, mean)
  p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
  breaks1 = seq(from = min(k), to = max(k), length = nBreaks + 1)
  hist1 = hist(k, breaks = breaks1, plot = FALSE, right = TRUE)
  dk2 = hist1$mids
  dk = ifelse(is.na(dk), dk2, dk)
  dk = ifelse(dk == 0, dk2, dk)
  p.dk = ifelse(is.na(p.dk), 0, p.dk)
  log.dk = as.vector(log10(dk))
  if (removeFirst) {
    p.dk = p.dk[-1]
    log.dk = log.dk[-1]
  }
  log.p.dk = as.numeric(log10(p.dk + 1e-09))
  lm1 = lm(log.p.dk ~ log.dk)
  lm2 = lm(log.p.dk ~ log.dk + I(10^log.dk))
  datout = data.frame(Rsquared.SFT = summary(lm1)$r.squared, slope.SFT = summary(lm1)$coefficients[2, 1], truncatedExponentialAdjRsquared = summary(lm2)$adj.r.squared)
  datout
}
# ======================= helpers =======================
.get_col <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(nm)
  NA_character_
}

.first_crossing_beta <- function(df, thr = 0.80, require_nonpos_slope = TRUE) {
  if (!nrow(df)) return(NA_integer_)
  r2_col    <- .get_col(df, c("SFT.R.sq", "Rsquared.SFT"))
  slope_col <- .get_col(df, c("slope", "slope.SFT"))
  if (is.na(r2_col)) return(NA_integer_)
  df <- df[order(df$Power), , drop = FALSE]
  ok <- is.finite(df[[r2_col]]) & df[[r2_col]] >= thr
  if (require_nonpos_slope && !is.na(slope_col)) {
    ok <- ok & is.finite(df[[slope_col]]) & (df[[slope_col]] <= 0)
  }
  idx <- which(ok)
  if (length(idx)) as.integer(df$Power[idx[1]]) else NA_integer_
}

# return named vectors/lists of crossings
.compute_crossings_from_fit_curves <- function(beta_info, tissues, thr = 0.80, require_nonpos_slope = TRUE) {
  ts_df <- beta_info$TS_fit_curves
  ct_df <- beta_info$CT_fit_curves

  # TS crossings per tissue
  TS_cross <- setNames(rep(NA_integer_, length(tissues)), tissues)
  if (!is.null(ts_df) && nrow(ts_df)) {
    for (t in tissues) {
      cols <- c("Power",
          .get_col(ts_df, c("SFT.R.sq","Rsquared.SFT")),
                .get_col(ts_df, c("slope","slope.SFT")))
      cols <- cols[!is.na(cols) & cols %in% names(ts_df)]
      sub  <- ts_df[ts_df$tissue == t, cols, drop = FALSE]
      TS_cross[t] <- .first_crossing_beta(sub, thr = thr, require_nonpos_slope = require_nonpos_slope)
    }
  }

  # CT crossings per pair ("A||B" names)
  CT_cross <- numeric(0)
  if (!is.null(ct_df) && nrow(ct_df)) {
    pairs <- unique(ct_df$pair)
    for (pr in pairs) {
      cols <- c("Power",
          .get_col(ct_df, c("Rsquared.SFT","SFT.R.sq")),
          .get_col(ct_df, c("slope.SFT","slope")))
      cols <- cols[!is.na(cols) & cols %in% names(ct_df)]
      sub  <- ct_df[ct_df$pair == pr, cols, drop = FALSE]
      CT_cross[pr] <- .first_crossing_beta(sub, thr = thr, require_nonpos_slope = require_nonpos_slope)
    }
  }
  list(TS = TS_cross, CT = CT_cross)
}

# Build beta matrix strictly from *crossings* (diag = TS, off-diag = CT)
build_beta_matrix_from_crossings <- function(beta_info, tissues, thr = 0.80, require_nonpos_slope = TRUE) {
  crosses <- .compute_crossings_from_fit_curves(beta_info, tissues, thr, require_nonpos_slope)
  B <- matrix(NA_real_, nrow = length(tissues), ncol = length(tissues),
              dimnames = list(tissues, tissues))

  # diag from TS crossings
  for (t in tissues) B[t, t] <- crosses$TS[[t]]

  # off-diagonals from CT crossings
  if (!is.null(crosses$CT) && length(crosses$CT)) {
    for (nm in names(crosses$CT)) {
      ab <- strsplit(nm, "\\|\\|")[[1]]
      if (length(ab) != 2) next
      a <- ab[1]; b <- ab[2]
      if (a %in% tissues && b %in% tissues) {
        B[a, b] <- crosses$CT[[nm]]
        B[b, a] <- crosses$CT[[nm]]
      }
    }
  }
  B
}


.require_or_stop <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = TRUE, quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse=", "),
                         "\nInstall: install.packages(c(", 
                         paste(sprintf('"%s"', miss), collapse=", "), "))")
}

.tissue_palette <- function(tissues) {
  .require_or_stop(c("RColorBrewer"))
  pal <- RColorBrewer::brewer.pal(9, "Set1")[c(1:5,7:9)]
  if (length(tissues) <= length(pal)) {
    setNames(pal[seq_along(tissues)], tissues)
  } else {
    setNames(rep(pal, length.out = length(tissues)), tissues)
  }
}

.split_pair <- function(x) strsplit(x, "\\|\\|", fixed = FALSE)

.ct_partner <- function(pair, t) {
  ab <- .split_pair(pair)[[1]]
  if (length(ab) != 2) return(NA_character_)
  if (identical(ab[1], t)) ab[2] else if (identical(ab[2], t)) ab[1] else NA_character_
}

.signed_R2_TS <- function(df) {
  slope_col <- if ("slope" %in% names(df)) "slope" else if ("slope.SFT" %in% names(df)) "slope.SFT" else NA
  if (is.na(slope_col)) return(df$SFT.R.sq)
  -sign(df[[slope_col]]) * df$SFT.R.sq
}

.signed_R2_CT <- function(df) {
  r2_col <- if ("Rsquared.SFT" %in% names(df)) "Rsquared.SFT" else if ("SFT.R.sq" %in% names(df)) "SFT.R.sq" else NA
  slope_col <- if ("slope.SFT" %in% names(df)) "slope.SFT" else if ("slope" %in% names(df)) "slope" else NA
  if (is.na(r2_col)) return(rep(NA_real_, nrow(df)))
  if (is.na(slope_col)) return(df[[r2_col]])
  -sign(df[[slope_col]]) * df[[r2_col]]
}

# ======================= 1) Beta series per tissue =======================
plot_beta_series_pdf <- function(beta_info, tissues, out_file = "plots/beta_series.pdf",
                                 vline_TS = NULL, vline_CT = NULL,
                                 layout_rows = 4, layout_cols = 2,
                                 thr = 0.80, require_nonpos_slope = TRUE) {
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  cols <- .tissue_palette(tissues)

  ts_df <- beta_info$TS_fit_curves
  ct_df <- beta_info$CT_fit_curves

  # find first β crossing per TS and per CT pair (R²≥thr & slope≤0)
  crosses <- .compute_crossings_from_fit_curves(beta_info, tissues, thr, require_nonpos_slope)

  # helpers to extract correct columns
  .r2_col_ts <- .get_col(ts_df, c("SFT.R.sq","Rsquared.SFT"))
  .slope_col_ts <- .get_col(ts_df, c("slope","slope.SFT"))
  .r2_col_ct <- .get_col(ct_df, c("Rsquared.SFT","SFT.R.sq"))

  pdf(out_file, height = 9, width = 5)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(layout_rows, layout_cols))

  for (t in tissues) {
    ts_t <- if (!is.null(ts_df) && nrow(ts_df)) ts_df[ts_df$tissue == t, , drop = FALSE] else ts_df[0,]
    if (!nrow(ts_t)) { plot.new(); title(t); next }

    # TS line (black)
    plot(ts_t$Power, ts_t[[.r2_col_ts]], type = "l", lwd = 2, col = "black",
         ylim = c(0, 1), main = t,
         xlab = expression(beta), ylab = expression(R^2))
    points(ts_t$Power, ts_t[[.r2_col_ts]], pch = 16, cex = 0.8, col = "black")

    # add CT partner lines (colored)
    if (!is.null(ct_df) && nrow(ct_df)) {
      keep <- grepl(paste0("^", t, "\\|\\|"), ct_df$pair) | grepl(paste0("\\|\\|", t, "$"), ct_df$pair)
      ct_sub <- ct_df[keep, , drop = FALSE]
      if (nrow(ct_sub)) {
        partners <- vapply(ct_sub$pair, function(p) {
          ab <- strsplit(p, "\\|\\|")[[1]]; if (ab[1] == t) ab[2] else ab[1]
        }, character(1))
        split_list <- split(ct_sub, partners)
        for (p in names(split_list)) {
          df <- split_list[[p]]
          lines(df$Power, df[[.r2_col_ct]], col = cols[[p]], lwd = 1.5)
          points(df$Power, df[[.r2_col_ct]], col = cols[[p]], pch = 16, cex = 0.6)
        }
      }
    }

    # horizontal line at threshold for visual aid
    abline(h = thr, col = "grey70", lty = 3)

    # dashed verticals at crossing β*, colored like the curve they belong to
    b_ts <- crosses$TS[[t]]
    if (is.finite(b_ts)) abline(v = b_ts, col = "black", lty = 2, lwd = 1.5)

    if (!is.null(ct_df) && nrow(ct_df)) {
      # for every partner, find the pair β* and color the dashed line with that partner's color
      partners_all <- unique(unlist(strsplit(names(crosses$CT), "\\|\\|")))
      partners_all <- setdiff(partners_all, NA_character_)
      for (p in setdiff(partners_all, t)) {
        nm1 <- paste0(t, "||", p); nm2 <- paste0(p, "||", t)
        b_ct <- if (!is.na(crosses$CT[nm1])) crosses$CT[nm1] else crosses$CT[nm2]
        if (is.finite(b_ct)) abline(v = b_ct, col = cols[[p]], lty = 2, lwd = 1.2)
      }
    }
  }

  # color legend (partners)
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center", legend = tissues, col = cols[tissues], pch = 15, cex = 0.95, ncol = 2, bty = "n")
}

# ======================= 2) Slope & connectivity summaries =======================
plot_tissue_summaries_pdf <- function(beta_info, tissues,
                                      out_file = "plots/beta_series_tissue_summaries.pdf") {
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  cols <- .tissue_palette(tissues)
  ts_df <- beta_info$TS_fit_curves
  if (is.null(ts_df) || !nrow(ts_df)) {
    warning("TS_fit_curves is empty; skipping summaries plot.")
    return(invisible(NULL))
  }

  ts_df$signed_R2 <- .signed_R2_TS(ts_df)
  slope_col <- if ("slope" %in% names(ts_df)) "slope" else if ("slope.SFT" %in% names(ts_df)) "slope.SFT" else NA
  meank_col <- if ("mean.k" %in% names(ts_df)) "mean.k" else if ("mean.k." %in% names(ts_df)) "mean.k." else NA

  pdf(out_file, height = 16, width = 5)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(4,1))

  # 1) signed R²
  yr <- range(ts_df$signed_R2, na.rm = TRUE)
  i <- 1
  ts_t <- ts_df[ts_df$tissue == tissues[i], ]
  plot(ts_t$Power, ts_t$signed_R2, type = "b", pch = 16, col = cols[[tissues[i]]],
       ylim = yr, xlab = expression(beta), ylab = expression(R^2))
  if (length(tissues) > 1) {
    for (i in 2:length(tissues)) {
      ts_t <- ts_df[ts_df$tissue == tissues[i], ]
      points(ts_t$Power, ts_t$signed_R2, type = "b", pch = 16, col = cols[[tissues[i]]])
    }
  }

  # 2) slope
  if (!is.na(slope_col)) {
    yr <- range(ts_df[[slope_col]], -1, na.rm = TRUE)
    i <- 1
    ts_t <- ts_df[ts_df$tissue == tissues[i], ]
    plot(ts_t$Power, ts_t[[slope_col]], type = "b", pch = 16, col = cols[[tissues[i]]],
         ylim = yr, xlab = expression(beta), ylab = "Slope")
    if (length(tissues) > 1) {
      for (i in 2:length(tissues)) {
        ts_t <- ts_df[ts_df$tissue == tissues[i], ]
        points(ts_t$Power, ts_t[[slope_col]], type = "b", pch = 16, col = cols[[tissues[i]]])
      }
    }
  } else {
    plot.new(); title("Slope (missing in TS_fit_curves)")
  }

  # 3) mean connectivity
  if (!is.na(meank_col)) {
    yr <- range(ts_df[[meank_col]], na.rm = TRUE)
    i <- 1
    ts_t <- ts_df[ts_df$tissue == tissues[i], ]
    plot(ts_t$Power, ts_t[[meank_col]], type = "b", pch = 16, col = cols[[tissues[i]]],
         ylim = yr, xlab = expression(beta), ylab = "Mean connectivity")
    if (length(tissues) > 1) {
      for (i in 2:length(tissues)) {
        ts_t <- ts_df[ts_df$tissue == tissues[i], ]
        points(ts_t$Power, ts_t[[meank_col]], type = "b", pch = 16, col = cols[[tissues[i]]])
      }
    }
    # 4) log10(mean connectivity)
    yr <- range(log10(ts_df[[meank_col]]), na.rm = TRUE)
    i <- 1
    ts_t <- ts_df[ts_df$tissue == tissues[i], ]
    plot(ts_t$Power, log10(ts_t[[meank_col]]), type = "b", pch = 16, col = cols[[tissues[i]]],
         ylim = yr, xlab = expression(beta), ylab = "Mean connectivity (log10)")
    if (length(tissues) > 1) {
      for (i in 2:length(tissues)) {
        ts_t <- ts_df[ts_df$tissue == tissues[i], ]
        points(ts_t$Power, log10(ts_t[[meank_col]]), type = "b", pch = 16, col = cols[[tissues[i]]])
      }
    }
  } else {
    plot.new(); title("Mean connectivity (missing)"); plot.new(); title("Mean connectivity (log10) (missing)")
  }

  legend("topright", legend = tissues, col = cols[tissues], pch = 16, cex = 1.0, bty = "n")
}

# ======================= 3) Beta matrix heatmap + CSV =======================
build_beta_matrix <- function(beta_info, tissues) {
  B <- matrix(NA_real_, nrow = length(tissues), ncol = length(tissues),
              dimnames = list(tissues, tissues))
  ts_map <- beta_info$TS_power_map
  if (is.null(ts_map) || all(is.na(ts_map))) {
    diag(B) <- beta_info$TS_power
  } else {
    for (t in tissues) B[t, t] <- as.numeric(ts_map[[t]])
  }
  ct_map <- beta_info$CT_power_map
  if (!is.null(ct_map) && length(ct_map)) {
    for (nm in names(ct_map)) {
      ab <- .split_pair(nm)[[1]]
      if (length(ab) != 2) next
      a <- ab[1]; b <- ab[2]
      if (a %in% tissues && b %in% tissues) {
        B[a, b] <- as.numeric(ct_map[[nm]])
        B[b, a] <- as.numeric(ct_map[[nm]])
      }
    }
  } else {
    for (i in seq_along(tissues)) for (j in seq_along(tissues)) if (i != j) B[i,j] <- beta_info$CT_power
  }
  B
}

plot_beta_matrix_heatmap <- function(beta_info, tissues,
                                     out_pdf = "plots/betas_heatmap.pdf",
                                     out_csv = "output/beta_matrix.csv",
                                     thr = 0.8, require_nonpos_slope = TRUE,
                                     palette_fun = function() grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Spectral")))(21)) {
  .require_or_stop(c("gplots", "RColorBrewer"))
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)

  B <- build_beta_matrix_from_crossings(beta_info, tissues, thr = thr, require_nonpos_slope = require_nonpos_slope)
  # --- cluster on an imputed copy, plot the original ---
  if (!any(is.finite(B))) stop("All entries in beta matrix are NA.")
  B_imp <- B
  med <- stats::median(B[is.finite(B)], na.rm = TRUE)
  B_imp[!is.finite(B_imp)] <- med  # simple, safe imputation for clustering

  # precompute dendrograms so heatmap.2 doesn't try to cluster NA's
  rowv <- as.dendrogram(hclust(dist(B_imp), method = "average"))
  colv <- as.dendrogram(hclust(dist(t(B_imp)), method = "average"))

  gplots::heatmap.2(
    B,
    Rowv = rowv, Colv = colv,      # <- use precomputed trees
    trace = "none",
    col = palette_fun(),
    cellnote = ifelse(is.finite(B), round(B, 2), NA),
    notecol = "black",
    margins = c(12, 14),
    key = TRUE,
    density.info = "none",
    main = bquote("First " * beta * " crossing (R"^2 * " ≥ " * .(thr) * ", slope ≤ 0)"),
    cexRow = 1.0, cexCol = 1.0,
    srtCol = 45,
    na.color = "white"
  )
  dev.flush()

  # Generate pdf file for the heatmap
  ggsave(out_pdf, width = 10, height = 8)

  write.csv(B, file = out_csv, row.names = TRUE)
  invisible(list(mat = B, pdf = out_pdf, csv = out_csv))

}

make_all_beta_plots <- function(beta_info, tissues,
                                out_prefix = "xwgcna",
                                series_pdf = file.path("plots", paste0(out_prefix, "_beta_series.pdf")),
                                summaries_pdf = file.path("plots", paste0(out_prefix, "_beta_series_tissue_summaries.pdf")),
                                heatmap_pdf = file.path("plots", paste0(out_prefix, "_betas_heatmap.pdf")),
                                heatmap_csv = file.path("output", paste0(out_prefix, "_beta_matrix.csv")),
                                thr = 0.80) {
  # Per-tissue series PDF (now with dashed β* at crossings, colored)
  plot_beta_series_pdf(beta_info, tissues, out_file = series_pdf, thr = thr, require_nonpos_slope = TRUE)

  # Summaries (unchanged)
  plot_tissue_summaries_pdf(beta_info, tissues, out_file = summaries_pdf)

  # Heatmap strictly from crossings, bigger page/margins
  plot_beta_matrix_heatmap(beta_info, tissues,
                           out_pdf = heatmap_pdf, out_csv = heatmap_csv,
                           thr = thr, require_nonpos_slope = TRUE)

  message("✓ Plots written:\n  - ", series_pdf,
          "\n  - ", summaries_pdf,
          "\n  - ", heatmap_pdf,
          "\n✓ Beta matrix CSV: ", heatmap_csv)
}

R2_connectivity <- function(adj_mat) {
  #--------------------------For scale Free-----------------------------
  k <- rowSums(adj_mat)
  r2 <- checkScaleFree_logbin(k)$Rsquared.SFT
  r2_conn <- list(r2=r2, mean_conn=mean(k), median_conn=median(k), max_conn=max(k), min_conn=min(k))
  
  #------------------For connectons-----------------------
  tissue_indexed <- c(0, as.numeric(table(unlist(lapply(strsplit(colnames(adj_mat), split = '_'), function(x) {return(x[1])})))))
  for(i in 2:length(tissue_indexed)) tissue_indexed[i] <- tissue_indexed[i] <- sum(tissue_indexed[i:(i-1)])
  
  TS_conn_mat <- matrix(NA, nrow=(length(tissue_indexed)-1), ncol=4)
  colnames(TS_conn_mat) <- c('mean', 'median', 'max', 'min')
  
  CT_conn_mat <- matrix(NA, nrow=((((length(tissue_indexed)-1)*(length(tissue_indexed)-1))-(length(tissue_indexed)-1))/2), ncol=4)
  colnames(CT_conn_mat) <- c('mean', 'median', 'max', 'min')
  CT_counter <- 1
  
  for(i in 1:(length(tissue_indexed)-1)) {
    temp_adj_mat <- adj_mat[(tissue_indexed[i]+1):tissue_indexed[i+1], (tissue_indexed[i]+1):tissue_indexed[i+1]]
    k <- rowSums(temp_adj_mat)
    TS_conn_mat[i, ] <- c(mean(k), median(k), max(k), min(k))
    
    if(i < length(tissue_indexed)-1) {
      for(j in 1:(length(tissue_indexed)-1-i)) {
        temp_adj_mat <- adj_mat[(tissue_indexed[i]+1):tissue_indexed[i+1], (tissue_indexed[i+j]+1):tissue_indexed[i+j+1]]
        k <- rowSums(temp_adj_mat)
        CT_conn_mat[CT_counter, ] <- c(mean(k), median(k), max(k), min(k))
        CT_counter <- CT_counter + 1
      }
    }
    
  }
  
  
  scale_free_conn_list <- list(r2_conn=r2_conn, TS_conn_mat=TS_conn_mat, CT_conn_mat=CT_conn_mat)
  # save(scale_free_conn_list, file='Scale_free_conn_list.RData')
  
  
  R2_Conn <- list(r2=r2_conn$r2, mean_conn=r2_conn$mean_conn, TS_mean_conn=paste(TS_conn_mat [ ,'mean'], collapse = ', '), mean_TS_mean_conn=mean(TS_conn_mat [ ,'mean']),
                  max_TS_mean_conn=max(TS_conn_mat [ ,'mean']), min_TS_mean_conn=min(TS_conn_mat [ ,'mean']), CT_mean_conn=paste(CT_conn_mat[ ,'mean'], collapse = ', '),
                  mean_CT_mean_conn=mean(CT_conn_mat[ ,'mean']), max_CT_mean_conn=max(CT_conn_mat[ ,'mean']), min_CT_mean_conn=min(CT_conn_mat[ ,'mean']))
  
  return(R2_Conn)
}

LoadExprData<-function(tissue_name, tissue_file_name, 
                       sd_quantile = 0.00,
                       max_genes_per_tissue = 5000
                       ) {
  
    if (grepl("\\.csv$", tissue_file_name)) {
        datExpr <- read.csv(tissue_file_name, check.names = FALSE)
        rownames(datExpr) <- datExpr[,1]
        datExpr <- as.matrix(datExpr[,-1,drop=FALSE])
    } else stop("Unsupported input file !!!")
  

    sds <- apply(datExpr, 2, sd, na.rm = TRUE)
    thr <- stats::quantile(sds, probs = sd_quantile, na.rm = TRUE)
    keep_sd <- sds >= thr
    datExpr <- datExpr[, keep_sd, drop = FALSE]
  
  
    if (ncol(datExpr) > max_genes_per_tissue) {
        o <- order(sds[colnames(datExpr)], decreasing = TRUE)
        datExpr <- datExpr[, o[seq_len(max_genes_per_tissue)], drop = FALSE]
    }
  
    colnames(datExpr) <- paste0(tissue_name, "_", colnames(datExpr))
    storage.mode(datExpr) <- "double"
    datExpr
}
  

AdjacencyFromExpr <- function(
    tissue_names = NULL,
    tissue_expr_file_names = NULL,
    sd_quantile = 0.00,
    max_genes_per_tissue = 5000,
    cor_method = "pearson",
    TS_power_map = NULL,
    CT_power_map = NULL,
    default_TS = 6L,
    default_CT = 3L
){
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)
  
  expr_list   <- vector("list", T)
  donor_mat <- vector("list", T)
  idx <- integer(T+1); rc_names <- character(0)

  
  # Load + filter each tissue; also precompute donor-aggregated matrices
  for (i in seq_len(T)) {
    X <- LoadExprData(
      tissue_names[i], tissue_expr_file_names[i],
      sd_quantile = sd_quantile,
      max_genes_per_tissue = max_genes_per_tissue
    )
    expr_list[[i]] <- X
    donor_mat[[i]] <- .aggregate_by_donor(X)
    idx[i + 1] <- idx[i] + ncol(X)
    rc_names <- c(rc_names, colnames(X))
  }
  
  # Allocate adjacency
  A <- matrix(0, nrow = idx[T + 1],
                    ncol = idx[T + 1],
                    dimnames = list(rc_names, rc_names))
  
  # Within-tissue blocks (TS)
  for (i in 1:T) {
    rows <- (idx[i] + 1):idx[i + 1]
    pow_i <- if(!is.null(TS_power_map) && !is.na(TS_power_map[tissue_names[i]])) {
      as.integer(TS_power_map[tissue_names[i]])
    } else {
      default_TS
    }
    Sii <- abs(cor(
      expr_list[[i]], use = "pairwise.complete.obs", method = cor_method
    ))
    A[rows, rows] <- Sii^pow_i
  }
  
  # Cross-tissue blocks (CT) via donor intersection
  if (T >= 2) {
    for (i in 1:(T - 1)) {
      rows <- (idx[i] + 1):idx[i + 1]
      di <- rownames(donor_mat[[i]])
      for (j in (i + 1):T) {
        cols <- (idx[j] + 1):idx[j + 1]
        pair_ij <- paste(tissue_names[i], tissue_names[j], sep = "||")
        pair_ji <- paste(tissue_names[j], tissue_names[i], sep = "||")
        pow_ij <- if(!is.null(CT_power_map)) {
          if(!is.na(CT_power_map[pair_ij])){
            as.integer(CT_power_map[pair_ij])
          } else if(!is.na(CT_power_map[pair_ji])){
            as.integer(CT_power_map[pair_ji])
          } else {
            default_CT
          }
        } else {
          default_CT
        }

        dj <- rownames(donor_mat[[j]])
        common <- intersect(di, dj)

        if (length(common) < 3) {
          msg <- sprintf("No/too-few common donors between %s and %s (|common|=%d).",
                         tissue_names[i], tissue_names[j], length(common))
          stop(msg, " Consider lowering filters or using consensus WGCNA.")
          
        }

        Mi <- donor_mat[[i]][common, , drop = FALSE]
        Mj <- donor_mat[[j]][common, , drop = FALSE]
        Sij <- abs(cor(Mi, Mj, use = "pairwise.complete.obs", method = cor_method))
        C <- Sij^pow_ij
        A[rows, cols] <- C
        A[cols, rows] <- t(C)
      }
    }
  }

  saveRDS(A, "adj_mat.rds")
  A
}

network_heatmap_gplots <- function(
  TOM_mat, dynamicColors, restGenes = NULL,
  out_file = "TOM_heatmap_gplots.png", tom_power = 4,
  palette_fun = function() {
    gplots::colorpanel(250, 'red','orange','lemonchiffon')
  }
) {
  if (!requireNamespace("gplots", quietly = TRUE))
    stop("Package 'gplots' is required but not installed.")

  if (is.null(restGenes)) {
    keep_idx <- seq_len(ncol(TOM_mat))
  } else if (is.logical(restGenes)){
    keep_idx <- which(restGenes)
  } else {
    keep_idx <- restGenes
  }
  if (length(keep_idx) < 4L) {
    warning("Too few genes to plot (<4).")
    return(invisible(NULL))
  }

  subTOM <- TOM_mat[keep_idx, keep_idx, drop = FALSE]
  colors_sub <- dynamicColors[keep_idx]

  diss <- 1 - subTOM
  hc <- hclust(as.dist(diss), method = "average")
  rowv <- as.dendrogram(hc)
  colv <- as.dendrogram(hc)

  plotM <- (1 - subTOM)^ tom_power
  diag(plotM) <- NA

  png(out_file, width = 1200, height = 1200, res = 140)
  gplots::heatmap.2(plotM, Rowv = rowv, Colv = colv, dendrogram = "both",
                    trace = "none", density.info = "none", scale = "none",
                    col = palette_fun(), key = TRUE, key.title = "Scale", 
                    key.xlab = "Dissimilarity^power", symm = TRUE, symkey = FALSE,
                    margins = c(6,6), labRow = FALSE, labCol = FALSE,
                    RowSideColors = colors_sub, ColSideColors = colors_sub,
                    lhei = c(0.15, 0.85),
                    lwid = c(0.15, 0.85)
  )
  legend("topright", legend = unique(colors_sub), fill = unique(colors_sub), border = NA,
         bty = "n", cex = 0.8, title = "Modules")
  dev.off()

  invisible(out_file)

}

network_heatmap_WGCNA <- function(
  TOM_mat, dynamicColors, restGenes = NULL,
  out_file = "TOM_heatmap_WGCNA.png", tom_power = 2,
  palette_fun = function() {
    gplots::colorpanel(250, 'red', 'orange', 'lemonchiffon')
  }
){
  if (!requireNamespace("gplots", quietly = TRUE))
    stop("Package 'gplots' is required but not installed.")
  
  if (is.null(restGenes)) {
    keep_idx <- seq_len(ncol(TOM_mat))
  } else if (is.logical(restGenes)){
    keep_idx <- which(restGenes)
  } else {
    keep_idx <- restGenes
  }

  if (length(keep_idx) < 4L) {
    warning("Too few genes to plot (<4).")
    return(invisible(NULL))
  }

  subTOM <- TOM_mat[keep_idx, keep_idx, drop = FALSE]
  colors_sub <- dynamicColors[keep_idx]

  diss <- 1 - subTOM
  hc <- hclust(as.dist(diss), method = "average")

  plotTOM <- diss^ tom_power
  diag(plotTOM) <- NA

  png(out_file, width = 1200, height = 1200, res = 140)
  WGCNA::TOMplot(plotTOM, dendro = hc, moduleColors = colors_sub, main = "Network heatmap (WGCNA::TOMplot)",
  col = palette_fun())
  dev.off()

  invisible(out_file)
}

network_heatmap <- function(restGenes, dynamicColors, tissue_names, tissue_expr_file_names) { 

  library(gplots)
  myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
  adj_mat_filter <- AdjacencyFromExpr(tissue_names, tissue_expr_file_names, 0.0)

  TOM_mat_filter <- Cross_Tissue_TOM(adj_mat_filter[restGenes, restGenes])
  rm(adj_mat_filter)
  gc()
  
  dist_TOM_mat_filter <- 1-TOM_mat_filter
  rm(TOM_mat_filter)
  gc()
  
  h_TOM_filter <- hclust(as.dist(dist_TOM_mat_filter), method="average")
  plotTOM = dist_TOM_mat_filter^4; 
  # Set diagonal to NA for a nicer plot 
  diag(plotTOM) = NA; 
  
  # Call the plot function 
  sizeGrWindow(9,9)
  png(filename="TOMplot.png")
  TOMplot(plotTOM, h_TOM_filter, as.character(dynamicColors[restGenes]),main = "Network heatmap plot, module genes", col=myheatcol)
  
  #plot()
  TOMplot
  dev.off()
  
}



Cross_Tissue_TOM <- function(adj_mat, block_size = 2000L) {
    message("Computing TOM from adjacency matrix...")
    n <- nrow(adj_mat)
    rn <- rownames(adj_mat)
    cn <- colnames(adj_mat)
    k <- rowSums(adj_mat)
    message("Number of genes: ", n, "; Number of samples: ", ncol(adj_mat))
    TOM_mat <- matrix(0, n, n)
    for (i in seq(1, n, by = block_size)){
        message("Processing block ", i, " to ", min(i + block_size - 1, n))
        i2 <- min(i + block_size - 1, n)
        Ai <- adj_mat[i:i2, , drop = FALSE]
        ki <- k[i:i2]

        for (j in seq(i, n, by = block_size)){
            j2 <- min(j + block_size - 1, n)
            Aj <- adj_mat[j:j2, , drop = FALSE]
            kj <- k[j:j2]

            num <- Ai %*% t(Aj)
            Aij <- adj_mat[i:i2, j:j2, drop = FALSE]
            denom <- outer(ki, kj, pmin) + 1 - Aij

            block <- (Aij + num) / denom
            TOM_mat[i:i2, j:j2] <- block
            if (j > i) TOM_mat[j:j2, i:i2] <- t(block)
            rm(Aj, num, Aij, denom, block)
            gc(FALSE)
        }
        message("Finished block ", i, " to ", i2)
        rm(Ai)
        gc(FALSE)
    }
    dimnames(TOM_mat) <- list(rn, cn)
    message("TOM calculation done.")
    TOM_mat
}

filter_plot_inputs_topVar <- function(
  TOM_mat, dynamicColors,
  tissue_names, tissue_expr_file_names,
  n_top = 5000,
  sd_quantile = 0.00,
  max_genes_per_tissue = 5000,
  aggregate_by_donor = FALSE
){
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  if (length(dynamicColors) != ncol(TOM_mat))
    stop("dynamicColors length must equal ncol(TOM_mat).")

  all_vars <- numeric(0)
  for (i in seq_along(tissue_names)) {
    X <- LoadExprData(
      tissue_name = tissue_names[i],
      tissue_file_name = tissue_expr_file_names[i],
      sd_quantile = sd_quantile,
      max_genes_per_tissue = max_genes_per_tissue
    )
    if (aggregate_by_donor) {
      X <- .aggregate_by_donor(X)  
    }
    v <- apply(X, 2, stats::var, na.rm = TRUE)
    names(v) <- colnames(X)       
    all_vars <- c(all_vars, v)
  }
  all_vars <- all_vars[is.finite(all_vars)]
  if (!length(all_vars)) stop("No genes found for variance ranking.")

  ord <- order(all_vars, decreasing = TRUE)
  top_names <- names(all_vars)[ord[seq_len(min(n_top, length(ord)))]]
  keep_names <- intersect(colnames(TOM_mat), top_names)
  keep_idx   <- match(keep_names, colnames(TOM_mat))

  if (length(keep_idx) < 4L)
    warning("Fewer than 4 genes kept — plots may not render well.")

  TOM_sub <- TOM_mat[keep_idx, keep_idx, drop = FALSE]
  colors_sub <- dynamicColors[keep_idx]

  invisible(list(
    TOM_mat      = TOM_sub,
    dynamicColors= colors_sub,
    keep_idx     = keep_idx,
    keep_names   = keep_names
  ))
}

Clusters_Table <- function(TOM_mat, minClusterSize = 30, plot_heatmap = FALSE, tissue_names, tissue_expr_file_names, group) {
    library(dynamicTreeCut)

    gene_names <- colnames(TOM_mat)
    if (is.null(gene_names)) gene_names <- rownames(TOM_mat)

    dist_TOM_mat <- 1 - TOM_mat

    h_TOM <- hclust(as.dist(dist_TOM_mat), method = "average")

    sizeGrWindow(12, 9)
    plot(h_TOM, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", 
    labels = FALSE, hang = 0.04)

    dynamicMods <- cutreeDynamic(h_TOM, method = "tree", deepSplit = TRUE, 
                                minClusterSize = minClusterSize)
    dynamicColors <- labels2colors(dynamicMods)
    
    restGenes <- (dynamicColors != "grey")
    message("Plot heatmap", sum(restGenes))


    network_heatmap_WGCNA(
      TOM_mat = TOM_mat,
      dynamicColors = dynamicColors,
      restGenes = restGenes,
      out_file = paste0("TOM_heatmap_WGCNA_", group, ".png"),
      tom_power = 2
    )

    module_labels <- sort(setdiff(unique(as.integer(dynamicMods)), 0L))
    if (length(module_labels) == 0L) {
        warning("No non-grey modules detected.")
        empty <- matrix(nrow = 0, ncol = 3)
        colnames(empty) <- c("Cluster ID", "Tissue", "Gene Symbol")
        return(empty)
    }

    split_pos <- regexpr("_", gene_names, fixed = TRUE)
    tissue_vec <- ifelse(split_pos > 0, 
                        substr(gene_names, 1, split_pos - 1), 
                        "Unknown")
    gene_vec <- ifelse(split_pos > 0, 
                       substr(gene_names, split_pos + 1, nchar(gene_names)), 
                       gene_names)

    df_list <- vector("list", length(module_labels))
    for (k in seq_along(module_labels)) {
        lab <- module_labels[k]
        idx <- which(dynamicMods == lab)
        if (length(idx) == 0L) next

        df_list[[k]] <- data.frame(
            "Cluster ID" = rep.int(k, length(idx)),
            "Tissue" = tissue_vec[idx],
            "Gene Symbol" = gene_vec[idx],
            stringAsFactors = FALSE,
            check.names = FALSE
        )
    }

    cluster_table <- do.call(rbind, Filter(Negate(is.null), df_list))
    cluster_table <- as.matrix(cluster_table)
    return(cluster_table)
}

Clusters_Details <- function(clusters_table, cluster_type_thr = 0.95) {
  if (is.null(clusters_table) || nrow(clusters_table) == 0) {
    warning("clusters_table is empty; returning empty details.")
    # Build an empty details table with expected columns
    return(matrix(nrow = 0, ncol = 0))
  }
  
  # Normalize column names: turn dots into spaces
  cn <- colnames(clusters_table)
  colnames(clusters_table) <- sub("\\.", " ", cn)
  
  required <- c("Cluster ID", "Tissue", "Gene Symbol")
  missing  <- setdiff(required, colnames(clusters_table))
  if (length(missing) > 0) {
    stop("clusters_table missing columns: ", paste(missing, collapse = ", "))
  }
  
  tissue_names   <- names(table(clusters_table[, "Tissue"]))
  total_clusters <- max(as.numeric(clusters_table[, "Cluster ID"]))
  
  clusters_table_details <- matrix(0, nrow = total_clusters, ncol = (5 + length(tissue_names)))
  colnames(clusters_table_details) <- c('Cluster ID', 'Cluster Size', 'Cluster Type',
                                        'Cluster Tissues', tissue_names, 'Dominant Tissue')
  
  for (i in 1:total_clusters) {
    temp_cluster <- clusters_table[as.numeric(clusters_table[, "Cluster ID"]) == i, , drop = FALSE]
    
    clusters_table_details[i, "Cluster ID"]   <- i
    clusters_table_details[i, "Cluster Size"] <- nrow(temp_cluster)
    
    if (nrow(temp_cluster) == 0) {
      clusters_table_details[i, "Cluster Type"] <- 'NA'
      next
    }
    
    prop_tbl <- round(table(temp_cluster[, "Tissue"]) / nrow(temp_cluster), 2)
    clusters_table_details[i, "Cluster Type"] <- ifelse(max(prop_tbl) >= cluster_type_thr, 'TS', 'CT')
    
    clusters_table_details[i, "Cluster Tissues"] <- paste(names(table(temp_cluster[, "Tissue"])), collapse = ',')
    clusters_table_details[i, names(table(temp_cluster[, "Tissue"]))] <- table(temp_cluster[, "Tissue"])
    clusters_table_details[i, "Dominant Tissue"] <- names(which.max(table(temp_cluster[, "Tissue"])))
  }
  
  clusters_table_details
}



# ====================== STARNET-style soft-threshold search (TS & CT) ======================

.cor_mat_TS <- function(X, cor_method) {
  cor(X, use = "pairwise.complete.obs", method = cor_method)
}

.cor_mat_CT <- function(Mi, Mj, cor_method) {
  cor(Mi, Mj, use = "pairwise.complete.obs", method = cor_method)
}

pickSoftThreshold_withinTissue <- function(
  expr,                        
  powerVector = c(1:10, seq(12, 20, 2)),
  cor_method  = "pearson",
  unsigned    = TRUE,        
  nBreaks     = 10,
  removeFirst = FALSE,
  use_signed_R2 = FALSE       
){
  S <- .cor_mat_TS(expr, cor_method)
  if (unsigned) S <- abs(S)

  fit_df <- data.frame(
    Power    = powerVector,
    SFT.R.sq = NA_real_,
    slope    = NA_real_,
    mean.k   = NA_real_,
    median.k = NA_real_,
    max.k    = NA_real_,
    stringsAsFactors = FALSE
  )

  for (ix in seq_along(powerVector)) {
    b <- powerVector[ix]
    B <- S ^ b
    diag(B) <- 0                    
    k <- rowSums(B)                

    cf <- checkScaleFree_logbin(k, nBreaks = nBreaks, removeFirst = removeFirst)
    r2    <- cf$Rsquared.SFT[1]
    slope <- cf$slope.SFT[1]

    fit_df$slope[ix]    <- slope
    fit_df$SFT.R.sq[ix] <- if (use_signed_R2) sign(slope) * r2 else r2
    fit_df$mean.k[ix]   <- mean(k)
    fit_df$median.k[ix] <- stats::median(k)
    fit_df$max.k[ix]    <- max(k)
  }
  list(fitIndices = fit_df)
}
.first_bin_mass <- function(k, nBreaks) {
  br <- seq(min(k, na.rm=TRUE), max(k, na.rm=TRUE), length.out = nBreaks + 1L)
  if (!is.finite(br[1]) || !is.finite(br[2])) return(0)
  mean(k >= br[1] & k <= br[2], na.rm = TRUE)
}
pickSoftThreshold_crossTissue <- function(
  Mi, Mj,                       
  powerVector = c(1:10, seq(12, 20, 2)),
  cor_method  = "pearson",
  unsigned    = TRUE,
  nBreaks     = 10,
  removeFirst = TRUE,
  use_signed_R2 = FALSE
){
  stopifnot(nrow(Mi) == nrow(Mj))
  S <- .cor_mat_CT(Mi, Mj, cor_method)
  if (unsigned) S <- abs(S)

  fit_df <- data.frame(
    Power    = powerVector,
    SFT.R.sq = NA_real_,
    slope    = NA_real_,
    mean.k   = NA_real_,
    median.k = NA_real_,
    max.k    = NA_real_,
    stringsAsFactors = FALSE
  )

  for (ix in seq_along(powerVector)) {
    b <- powerVector[ix]
    B <- S ^ b
    B[!is.finite(B)] <- 0
    # k <- c(rowSums(B), colSums(B))
    k <- c(rowSums(B) / ncol(B), colSums(B) / nrow(B))  # normalize by number of genes in other tissue
    k <- k * min(nrow(B), ncol(B))  # rescale to original scale
    k[!is.finite(k)] <- 0

    # rmFirst <- removeFirst || (.first_bin_mass(k, nBreaks) > 0.4)
    prop_zero <- mean(k < 1e-6, na.rm=TRUE)
    message(sprintf("[β=%a] mean(k)=%.3f, median(k)=%.3f, max(k)=%.3f, zeros=%.1f%%\n",
              b, mean(k), median(k), max(k), 100*prop_zero))
    cf <- checkScaleFree_logbin(k, nBreaks = nBreaks, removeFirst = removeFirst)
    r2    <- cf$Rsquared.SFT[1]
    slope <- cf$slope.SFT[1]

    fit_df$slope[ix]    <- slope
    fit_df$SFT.R.sq[ix] <- if (use_signed_R2) sign(slope) * r2 else r2
    fit_df$mean.k[ix]   <- mean(k)
    fit_df$median.k[ix] <- stats::median(k)
    fit_df$max.k[ix]    <- max(k)
    prop_zero <- mean(k < 1e-6, na.rm=TRUE)
    message(sprintf("[β=%a] mean(k)=%.3f, median(k)=%.3f, max(k)=%.3f, zeros=%.1f%%\n",
              b, mean(k), median(k), max(k), 100*prop_zero))
  }
  list(fitIndices = fit_df)
}


# =============================== auto_pick_powers ===============================
.choose_power_from_pickSoft <- function(sft, targetR2 = 0.80) {
  df <- sft$fitIndices
  ok <- which(df$SFT.R.sq >= targetR2 & df$slope < 0)
  if (length(ok) > 0) {
    return(min(df$Power[ok]))
  } else {
    return(df$Power[which.max(df$SFT.R.sq)])
  }
}

auto_pick_powers <- function(
  tissue_names,
  tissue_expr_file_names,
  sd_quantile = 0.90,
  max_genes_per_tissue = 4000,
  cor_method = "pearson",                
  beta_grid = c(1:10, seq(12, 20, 2)),
  targetR2 = 0.80,
  unsigned        = TRUE,
  nBreaks         = 10,
  removeFirst     = TRUE,
  use_signed_R2_TS = FALSE,
  use_signed_R2_CT = FALSE
){
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)

  expr_list   <- vector("list", T)
  donors_list <- vector("list", T)
  for (i in seq_len(T)) {
    X <- LoadExprData(
      tissue_name      = tissue_names[i],
      tissue_file_name = tissue_expr_file_names[i],
      sd_quantile      = sd_quantile,
      max_genes_per_tissue = max_genes_per_tissue
    )
    expr_list[[i]]   <- X
    donors_list[[i]] <- .aggregate_by_donor(X)
  }

  TS_per_tissue <- integer(T)
  TS_fit_curves <- vector("list", T)

  TS_power_map <- setNames(integer(T), tissue_names)
  for (i in seq_len(T)) {
    sft_ts <- tryCatch(
      pickSoftThreshold_withinTissue(
        expr            = expr_list[[i]],
        powerVector     = beta_grid,
        cor_method      = cor_method,
        unsigned        = unsigned,
        nBreaks         = nBreaks,
        removeFirst     = removeFirst,
        use_signed_R2   = use_signed_R2_TS
      ),
      error = function(e) { message("[TS beta pick] ", tissue_names[i], " error: ", e$message); NULL }
    )
    if (!is.null(sft_ts)) {
      fi <- sft_ts$fitIndices
      fi$tissue <- tissue_names[i]
      TS_fit_curves[[i]] <- fi
      TS_per_tissue[i] <- .choose_power_from_pickSoft(sft_ts, targetR2 = targetR2)
      TS_power_map[tissue_names[i]] <- TS_per_tissue[i]
    } else {
      TS_fit_curves[[i]] <- data.frame(Power = beta_grid,
                                       SFT.R.sq = NA_real_,
                                       slope = NA_real_,
                                       mean.k = NA_real_,
                                       median.k = NA_real_,
                                       max.k = NA_real_,
                                       tissue = tissue_names[i])
      TS_per_tissue[i] <- 6L
    }
  }
  TS_fit_curves_df <- do.call(rbind, TS_fit_curves)
  TS_power <- as.integer(stats::median(TS_per_tissue, na.rm = TRUE))

  CT_betas    <- integer(0)
  CT_fit_list <- list()

  CT_power_map <- numeric(0)
  if (T >= 2) {
    for (i in 1:(T - 1)) {
      for (j in (i + 1):T) {
        di <- rownames(donors_list[[i]]); dj <- rownames(donors_list[[j]])
        common <- intersect(di, dj)
        pair_id <- paste(tissue_names[i], tissue_names[j], sep = "||")

        if (length(common) < 3) {
          message(sprintf("[CT beta pick] %s: |common donors|=%s → skip", pair_id, length(common)))
          next
        }
        Mi <- donors_list[[i]][common, , drop = FALSE]
        Mj <- donors_list[[j]][common, , drop = FALSE]

        sft_ct <- tryCatch(
          pickSoftThreshold_crossTissue(
            Mi, Mj,
            powerVector   = beta_grid,
            cor_method    = cor_method,
            unsigned      = unsigned,
            nBreaks       = nBreaks,
            removeFirst   = removeFirst,
            use_signed_R2 = use_signed_R2_CT
          ),
          error = function(e) { message("[CT beta pick] ", pair_id, " error: ", e$message); NULL }
        )

        if (!is.null(sft_ct)) {
          best_beta <- .choose_power_from_pickSoft(sft_ct, targetR2 = targetR2)
          CT_betas <- c(CT_betas, best_beta)
          CT_power_map[pair_id] <- best_beta
          message(sprintf("[CT beta pick] %s: beta=%s", pair_id, best_beta))

          fi <- sft_ct$fitIndices
          CT_fit_list[[length(CT_fit_list) + 1]] <- data.frame(
            pair = pair_id,
            Power = fi$Power,
            Rsquared.SFT = fi$SFT.R.sq,
            stringsAsFactors = FALSE
          )
        } else {
          CT_fit_list[[length(CT_fit_list) + 1]] <-
            data.frame(pair = pair_id, Power = beta_grid, Rsquared.SFT = NA_real_)
        }
      }
    }
  }

  CT_power <- if (length(CT_betas) == 0) {
    message("[CT beta pick] no valid CT pairs → CT_power=3 (fallback)")
    3L
  } else {
    as.integer(stats::median(CT_betas, na.rm = TRUE))
  }

  CT_fit_curves_df <- if (length(CT_fit_list)) do.call(rbind, CT_fit_list) else
    data.frame(pair = character(), Power = integer(), Rsquared.SFT = numeric())

  list(
    TS_power       = TS_power,
    CT_power       = CT_power,
    TS_per_tissue  = TS_per_tissue,
    CT_per_pair    = CT_betas,
    TS_fit_curves  = TS_fit_curves_df,
    CT_fit_curves  = CT_fit_curves_df,
    TS_power_map   = TS_power_map,
    CT_power_map   = CT_power_map
  )
}


plot_beta_curves_per_tissue <- function(
  beta_info,
  out_prefix = "xwgcna",
  targetR2   = NULL,   
  save_png   = TRUE,
  width = 8, height = 5, dpi = 150
){
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; skipping plots.")
    return(invisible(NULL))
  }
  library(ggplot2)

  ts_df <- beta_info$TS_fit_curves   
  ct_df <- beta_info$CT_fit_curves  

  if (is.null(ts_df) || !nrow(ts_df)) {
    warning("Empty TS_fit_curves; nothing to plot.")
    return(invisible(NULL))
  }

  tissues <- unique(ts_df$tissue)
  plots <- vector("list", length(tissues))
  names(plots) <- tissues

  .partner_of <- function(pair, t){
    a_b <- strsplit(pair, "\\|\\|", fixed = FALSE)
    sapply(a_b, function(x){
      if (length(x) != 2) return(NA_character_)
      if (identical(x[1], t)) x[2] else if (identical(x[2], t)) x[1] else NA_character_
    })
  }

  for (t in tissues) {
    ts_t <- ts_df[ts_df$tissue == t, c("Power","SFT.R.sq")]
    ct_sub <- NULL
    if (!is.null(ct_df) && nrow(ct_df)) {
      keep <- grepl(paste0("^", t, "\\|\\|"), ct_df$pair) | grepl(paste0("\\|\\|", t, "$"), ct_df$pair)
      ct_sub <- ct_df[keep, , drop = FALSE]
      if (nrow(ct_sub)) {
        ct_sub$partner <- .partner_of(ct_sub$pair, t)
      }
    }

    p <- ggplot() +
      geom_line(data = ts_t, aes(x = Power, y = SFT.R.sq), color = "black", linewidth = 1.1) +
      geom_point(data = ts_t, aes(x = Power, y = SFT.R.sq), color = "black", size = 1.6)

    if (!is.null(ct_sub) && nrow(ct_sub)) {
      p <- p +
        geom_line(data = ct_sub, aes(x = Power, y = Rsquared.SFT, color = partner), linewidth = 0.9, alpha = 0.9) +
        geom_point(data = ct_sub, aes(x = Power, y = Rsquared.SFT, color = partner), size = 1.2, alpha = 0.9)
    } else {
      p <- p + annotate("text",
                        x = min(ts_t$Power, na.rm = TRUE),
                        y = min(ts_t$SFT.R.sq, na.rm = TRUE),
                        label = "No CT pairs for this tissue",
                        hjust = 0, vjust = 0, size = 3.2)
    }

    if (!is.null(targetR2)) {
      p <- p + geom_hline(yintercept = targetR2, linetype = "dashed")
    }

    p <- p +
      labs(
        title = paste0("Scale-free fit (R²) vs β — ", t),
        x = "Soft-threshold (β)",
        y = "Scale-free fit R²",
        color = "CT partner"
      ) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.key.size = grid::unit(0.6, "lines"))

    if (save_png) {
      fn <- paste0(out_prefix, "_beta_curves_", gsub("[^A-Za-z0-9]+", "_", t), ".png")
      ggplot2::ggsave(fn, p, width = width, height = height, dpi = dpi)
    }
    plots[[t]] <- p
  }

  return(plots)
}

# =========================
# WGCNA-based β (soft-threshold) pickers for TS & CT
# Drop-in additions that mirror the STARNET script logic while
# fitting into your existing XWGCNA_* pipeline.
#
# Key points:
# - Uses WGCNA::pickSoftThreshold() **on expression matrices** (not adjacency),
#   exactly like the STARNET snippet.
# - Chooses the **lowest β** whose scale-free fit R^2 crosses a target (default 0.80);
#   otherwise takes the β with maximal R^2.
# - Returns the same structures your current `auto_pick_powers()` returns,
#   so you can swap methods with a single flag.
# - Respects your existing helpers: LoadExprData(), .aggregate_by_donor(), .extract_donor().
# - Preserves the original STARNET parameter style: power vector seq(0.5, 20, length.out = 20),
#   verbose = 5, networkType derived from TOMType.
#
# Suggested wiring:
# - Add parameter `beta_method = c("custom", "wgcna")` to XWGCNA_Clusters_autoBeta().
# - When auto_beta && beta_method == "wgcna", call `wgcna_auto_pick_powers()` below.
# - Everything else in your pipeline (AdjacencyFromExpr, TOM, clustering, plots)
#   can remain unchanged.
# =========================

suppressMessages(WGCNA::allowWGCNAThreads())

.networkType_from_TOMType <- function(TOMType) {
  tt <- tolower(TOMType %||% "unsigned")
  if (startsWith(tt, "signed")) {
    if (grepl("hybrid", tt)) return("signed hybrid")
    return("signed")
  }
  "unsigned"
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# --------------------------- chooser ---------------------------
.wgcna_choose_beta_min_above <- function(sft, targetR2 = 0.80, require_neg_slope = TRUE) {
  stopifnot(!is.null(sft), !is.null(sft$fitIndices))
  df <- sft$fitIndices
  if (!all(c("Power","SFT.R.sq") %in% names(df))) return(NA_integer_)
  ok <- !is.na(df$SFT.R.sq) & (df$SFT.R.sq >= targetR2)
  if (require_neg_slope && "slope" %in% names(df)) ok <- ok & (df$slope < 0)
  if (any(ok)) {
    return(as.integer(min(df$Power[ok])))
  }
  # fallback: argmax R^2
  i <- which.max(df$SFT.R.sq)
  as.integer(df$Power[i])
}

# --------------------------- TS picker ---------------------------
wgcna_pick_TS <- function(expr_mat,
                          powerVector = seq(0.5, 20, length.out = 20),
                          TOMType = "unsigned",
                          cor_method = "pearson",
                          verbose = 5,
                          targetR2 = 0.80,
                          require_neg_slope = TRUE,
                          ...) {
  stopifnot(is.matrix(expr_mat) || is.data.frame(expr_mat))
  # WGCNA expects samples x genes
  # LoadExprData() already returns samples as rows, genes as columns -> OK
  corFnc <- if (tolower(cor_method) == "pearson") "cor" else "bicor"
  corOptions <- if (identical(corFnc, "cor"))
    list(use = "pairwise.complete.obs")
  else
    list(use = "pairwise.complete.obs", maxPOutliers = 1, robustY = FALSE)

  sft <- WGCNA::pickSoftThreshold(
    data = expr_mat,
    powerVector = powerVector,
    networkType = .networkType_from_TOMType(TOMType),
    corFnc = corFnc,
    corOptions = corOptions,
    verbose = verbose,
    ...
  )

  chosen <- .wgcna_choose_beta_min_above(sft, targetR2 = targetR2,
                                         require_neg_slope = require_neg_slope)
  # Ensure a tidy fitIndices like your plotter expects
  fi <- sft$fitIndices
  if (!"SFT.R.sq" %in% names(fi) && "SFT.R.sq" %in% names(sft)) fi$SFT.R.sq <- sft$SFT.R.sq
  list(beta = chosen, sft = sft, fit_df = fi)
}

# --------------------------- CT picker ---------------------------

wgcna_pick_CT_new <- function(
  expr_A, expr_B,
  aggregate_by_donor = FALSE,
  powerVector = c(1:10, seq(12, 20, 2)),
  TOMType     = "unsigned",        # "unsigned" | "signed" | "signed hybrid"
  cor_method  = "pearson",         # "pearson" | "bicor"
  bicor_maxPOutliers = 1,
  bicor_robustY = FALSE,
  targetR2   = 0.80,
  require_neg_slope = TRUE,
  nBreaks    = 50,
  removeFirst = TRUE,
  min_common = 3L,
  verbose    = 1
){
  stopifnot(is.matrix(expr_A) || is.data.frame(expr_A))
  stopifnot(is.matrix(expr_B) || is.data.frame(expr_B))
  X_A <- expr_A; X_B <- expr_B

  if (aggregate_by_donor) {
    X_A <- .aggregate_by_donor(X_A)
    X_B <- .aggregate_by_donor(X_B)
  }

  common <- intersect(rownames(X_A), rownames(X_B))
  if (length(common) < min_common) {
    if (verbose) message(sprintf(
      "[wgcna_pick_CT] Too few common donors: %d < %d", length(common), min_common))
    return(list(
      beta = NA_integer_,
      fitIndices = data.frame(Power = powerVector, SFT.R.sq = NA_real_,
                              slope.SFT = NA_real_, mean.k = NA_real_,
                              median.k = NA_real_, max.k = NA_real_)
    ))
  }
  Mi <- X_A[common, , drop = FALSE]
  Mj <- X_B[common, , drop = FALSE]

  # קורלציה בין-רקמתית (genes_A x genes_B)
  if (tolower(cor_method) == "bicor") {
    S <- WGCNA::bicor(Mi, Mj, use = "pairwise.complete.obs",
                      maxPOutliers = bicor_maxPOutliers, robustY = bicor_robustY)
  } else {
    S <- stats::cor(Mi, Mj, use = "pairwise.complete.obs", method = "pearson")
  }
  S[!is.finite(S)] <- 0

  # המרת קורלציה ל-adjacency לפי רשת חתומה/לא חתומה
  .adj_from_cor <- function(S, beta, TOMType) {
    tt <- tolower(TOMType)
    if (startsWith(tt, "signed")) {
      A <- ((S + 1)/2) ^ beta         # שומר סימן (בטווח [0,1])
    } else {
      A <- abs(S) ^ beta               # unsigned
    }
    A[!is.finite(A)] <- 0
    A
  }

  fit_df <- data.frame(
    Power    = powerVector,
    SFT.R.sq = NA_real_,
    slope    = NA_real_,
    mean.k   = NA_real_,
    median.k = NA_real_,
    max.k    = NA_real_,
    stringsAsFactors = FALSE
  )

  for (ix in seq_along(powerVector)) {
    b <- powerVector[ix]
    A <- .adj_from_cor(S, b, TOMType)
    # דרגות משני הצדדים (דו-חלקי): שורות = גנים ברקמה A, עמודות = גנים ברקמה B
    k <- c(rowSums(A, na.rm = TRUE), colSums(A, na.rm = TRUE))
    k[!is.finite(k)] <- 0

    sf <- WGCNA::scaleFreeFitIndex(k, nBreaks = nBreaks, removeFirst = removeFirst)
    fit_df$SFT.R.sq[ix] <- sf$Rsquared.SFT
    fit_df$slope[ix]    <- sf$slope.SFT
    fit_df$mean.k[ix]   <- mean(k)
    fit_df$median.k[ix] <- stats::median(k)
    fit_df$max.k[ix]    <- max(k)
    if (verbose > 1) {
      message(sprintf("[β=%s] R2=%.3f, slope=%.3f, mean(k)=%.2f", b, sf$Rsquared.SFT, sf$slope.SFT, mean(k)))
    }
  }

  # בחירת β: הנמוך החוצה סף R^2 ובעל שיפוע שלילי (אם מבוקש); אחרת argmax R^2
  ok <- which(!is.na(fit_df$SFT.R.sq) & fit_df$SFT.R.sq >= targetR2)
  if (require_neg_slope) ok <- ok[fit_df$slope[ok] < 0]
  beta_chosen <- if (length(ok)) min(fit_df$Power[ok]) else fit_df$Power[which.max(fit_df$SFT.R.sq)]

  list(beta = as.integer(beta_chosen), fitIndices = fit_df)
}
# Builds a paired expression matrix (samples-in-common x [genes_A | genes_B])
# and runs WGCNA::pickSoftThreshold() on it, exactly like the STARNET snippet.
# wgcna_pick_CT <- function(expr_A, expr_B,
#                           aggregate_by_donor = FALSE,
#                           powerVector = seq(0.5, 20, length.out = 20),
#                           TOMType = "unsigned",
#                           cor_method = "pearson",
#                           verbose = 5,
#                           targetR2 = 0.80,
#                           require_neg_slope = TRUE,
#                           min_common = 3L,
#                           bicor_maxPOutliers = bicor_maxPOutliers,
#                           bicor_robustY = bicor_robustY,
#                           nBreaks = 12,
#                           removeFirst = TRUE,
#                           min_common = 3L,
#                           verbose = 5,
#                           ...) {
#   stopifnot(is.matrix(expr_A) || is.data.frame(expr_A))
#   stopifnot(is.matrix(expr_B) || is.data.frame(expr_B))

#   X_A <- expr_A; X_B <- expr_B
#   if (aggregate_by_donor) {
#     X_A <- .aggregate_by_donor(X_A)
#     X_B <- .aggregate_by_donor(X_B)
#   }

#   common <- intersect(rownames(X_A), rownames(X_B))
#   if (length(common) < min_common)
#     return(list(beta = NA_integer_, sft = NULL,
#                 fit_df = data.frame(Power = powerVector, SFT.R.sq = NA_real_)))

#   paired <- cbind(X_A[common, , drop = FALSE], X_B[common, , drop = FALSE])

#   corFnc <- if (tolower(cor_method) == "pearson") "cor" else "bicor"
#   corOptions <- if (identical(corFnc, "cor"))
#     list(use = "pairwise.complete.obs")
#   else
#     list(use = "pairwise.complete.obs", maxPOutliers = 1, robustY = FALSE)

#   sft <- WGCNA::pickSoftThreshold(
#     data = paired,
#     powerVector = powerVector,
#     networkType = .networkType_from_TOMType(TOMType),
#     corFnc = corFnc,
#     corOptions = corOptions,
#     verbose = verbose,
#     ...
#   )
#   chosen <- .wgcna_choose_beta_min_above(sft, targetR2 = targetR2,
#                                          require_neg_slope = require_neg_slope)
#   fi <- sft$fitIndices
#   if (!"SFT.R.sq" %in% names(fi) && "SFT.R.sq" %in% names(sft)) fi$SFT.R.sq <- sft$SFT.R.sq
#   list(beta = chosen, sft = sft, fit_df = fi)
# }

wgcna_pick_CT <- function(
  expr_A, expr_B,
  aggregate_by_donor = FALSE,
  powerVector = seq(0.5, 20, length.out = 20),
  TOMType = "unsigned",
  cor_method = "pearson",
  verbose = 5,
  targetR2 = 0.80,
  require_neg_slope = TRUE,
  min_common = 3L,
  bicor_maxPOutliers = 1,
  bicor_robustY = FALSE,
  ...
){
  stopifnot(is.matrix(expr_A) || is.data.frame(expr_A))
  stopifnot(is.matrix(expr_B) || is.data.frame(expr_B))

  XA <- expr_A; XB <- expr_B
  if (aggregate_by_donor) { XA <- .aggregate_by_donor(XA); XB <- .aggregate_by_donor(XB) }

  common <- intersect(rownames(XA), rownames(XB))
  if (length(common) < min_common)
    return(list(beta = NA_integer_, sft = NULL,
                fit_df = data.frame(Power = powerVector, SFT.R.sq = NA_real_)))

  paired <- cbind(XA[common, , drop = FALSE], XB[common, , drop = FALSE])

  corFnc <- if (tolower(cor_method) == "pearson") "cor" else "bicor"
  corOptions <- if (identical(corFnc, "cor")) {
    list(use = "pairwise.complete.obs")
  } else {
    list(use = "pairwise.complete.obs",
         maxPOutliers = bicor_maxPOutliers, robustY = bicor_robustY)
  }

  sft <- WGCNA::pickSoftThreshold(
    data = paired,
    powerVector = powerVector,
    networkType = .networkType_from_TOMType(TOMType),
    corFnc = corFnc,
    corOptions = corOptions,
    verbose = verbose
  )

  chosen <- .wgcna_choose_beta_min_above(sft, targetR2 = targetR2,
                                         require_neg_slope = require_neg_slope)
  fi <- sft$fitIndices
  if (!"SFT.R.sq" %in% names(fi) && "SFT.R.sq" %in% names(sft)) fi$SFT.R.sq <- sft$SFT.R.sq
  list(beta = as.integer(chosen), sft = sft, fitIndices = fi, fit_df = fi)
}



# --------------------------- Auto picker (TS & CT) ---------------------------
# Mirrors your current auto_pick_powers() return shape so you can plug it in.
wgcna_auto_pick_powers_new <- function(
 tissue_names,
  tissue_expr_file_names,
  sd_quantile = 0.50,
  max_genes_per_tissue = 5000,
  TOMType = "unsigned",
  cor_method = "pearson",                # "pearson" | "bicor"
  powerVector = seq(0.5, 20, length.out = 20),
  targetR2 = 0.80,
  require_neg_slope = TRUE,
  verbose = 5,
  # ---- CT-specific knobs ----
  aggregate_by_donor_CT = FALSE,
  min_common_CT = 3L,
  ct_nBreaks = 50,
  ct_removeFirst = TRUE,
  bicor_maxPOutliers = 1,
  bicor_robustY = FALSE
) {
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)

  # Load & prefilter expression per tissue (samples x genes)
  expr_list <- vector("list", T)
  names(expr_list) <- tissue_names
  message("Loading expression data for ", T, " tissues…")
  for (i in seq_len(T)) {
    X <- LoadExprData(
      tissue_name = tissue_names[i],
      tissue_file_name = tissue_expr_file_names[i],
      sd_quantile = sd_quantile,
      max_genes_per_tissue = max_genes_per_tissue
    )
    expr_list[[i]] <- X
  }

  # ---- TS per tissue ----
  TS_per_tissue <- integer(T)
  TS_power_map  <- setNames(integer(T), tissue_names)
  TS_fit_curves <- vector("list", T)

  for (i in seq_len(T)) {
    ts_fit <- tryCatch(
      wgcna_pick_TS(
        expr_mat          = expr_list[[i]],
        powerVector       = powerVector,
        TOMType           = TOMType,
        cor_method        = cor_method,
        verbose           = verbose,
        targetR2          = targetR2,
        require_neg_slope = require_neg_slope
      ),
      error = function(e) { message("[TS] ", tissue_names[i], ": ", e$message); NULL }
    )
    if (!is.null(ts_fit)) {
      TS_per_tissue[i] <- ts_fit$beta
      fi <- ts_fit$fit_df
      fi$tissue <- tissue_names[i]
      TS_fit_curves[[i]] <- fi
    } else {
      TS_per_tissue[i] <- NA_integer_
      TS_fit_curves[[i]] <- data.frame(Power = powerVector, SFT.R.sq = NA_real_,
                                       tissue = tissue_names[i])
    }
    TS_power_map[tissue_names[i]] <- TS_per_tissue[i]
  }
  TS_power <- as.integer(stats::median(TS_per_tissue, na.rm = TRUE))
  TS_fit_curves_df <- do.call(rbind, TS_fit_curves)

  # ---- CT per pair (using the NEW wgcna_pick_CT) ----
  CT_betas <- integer(0)
  CT_power_map <- numeric(0)
  CT_fit_list <- list()

  if (T >= 2) {
    for (i in 1:(T - 1)) {
      for (j in (i + 1):T) {
        pair_id <- paste(tissue_names[i], tissue_names[j], sep = "||")
        ct_fit <- tryCatch(
          wgcna_pick_CT_new(
            expr_A = expr_list[[i]], expr_B = expr_list[[j]],
            aggregate_by_donor = aggregate_by_donor_CT,
            powerVector = powerVector,
            TOMType     = TOMType,
            cor_method  = cor_method,
            bicor_maxPOutliers = bicor_maxPOutliers,
            bicor_robustY      = bicor_robustY,
            targetR2    = targetR2,
            require_neg_slope = require_neg_slope,
            nBreaks     = ct_nBreaks,
            removeFirst = ct_removeFirst,
            min_common  = min_common_CT,
            verbose     = verbose
          ),
          error = function(e) { message("[CT] ", pair_id, ": ", e$message); NULL }
        )

        if (!is.null(ct_fit) && !is.na(ct_fit$beta)) {
          CT_betas <- c(CT_betas, ct_fit$beta)
          CT_power_map[pair_id] <- ct_fit$beta

          fi <- ct_fit$fitIndices
          # שימור שם העמודה כדי להתאים לפונקציות שרטוט קיימות
          CT_fit_list[[length(CT_fit_list) + 1]] <- data.frame(
            pair = pair_id,
            Power = fi$Power,
            Rsquared.SFT = fi$SFT.R.sq,
            slope.SFT    = if ("slope" %in% names(fi)) fi$slope else fi$slope.SFT,
            mean.k       = fi$mean.k,
            median.k     = fi$median.k,
            max.k        = fi$max.k,
            stringsAsFactors = FALSE
          )
        } else {
          CT_fit_list[[length(CT_fit_list) + 1]] <-
            data.frame(pair = pair_id, Power = powerVector,
                       Rsquared.SFT = NA_real_, slope.SFT = NA_real_,
                       mean.k = NA_real_, median.k = NA_real_, max.k = NA_real_)
        }
      }
    }
  }

  CT_power <- if (length(CT_betas)) as.integer(stats::median(CT_betas, na.rm = TRUE)) else 3L
  CT_fit_curves_df <- if (length(CT_fit_list)) do.call(rbind, CT_fit_list) else
    data.frame(pair = character(), Power = numeric(), Rsquared.SFT = numeric())

  list(
    method         = "wgcna+scaleFreeFitIndex(CT)",
    TS_power       = TS_power,
    CT_power       = CT_power,
    TS_per_tissue  = TS_per_tissue,
    CT_per_pair    = as.integer(CT_betas),
    TS_fit_curves  = TS_fit_curves_df,
    CT_fit_curves  = CT_fit_curves_df,
    TS_power_map   = TS_power_map,
    CT_power_map   = CT_power_map
  )
}


wgcna_auto_pick_powers <- function(
  tissue_names,
  tissue_expr_file_names,
  sd_quantile = 0.50,
  max_genes_per_tissue = 5000,
  TOMType = "unsigned",
  cor_method = "pearson",                # "pearson" | "bicor"
  powerVector = seq(0.5, 20, length.out = 20),
  targetR2 = 0.80,
  require_neg_slope = TRUE,
  verbose = 5,
  # ---- CT-specific knobs ----
  aggregate_by_donor_CT = FALSE,
  min_common_CT = 3L,
  ct_nBreaks = 50,
  ct_removeFirst = TRUE,
  bicor_maxPOutliers = 1,
  bicor_robustY = FALSE
) {
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)

  # Load & prefilter expression per tissue (samples x genes)
  expr_list <- vector("list", T)
  names(expr_list) <- tissue_names
  message("Loading expression data for ", T, " tissues…")
  for (i in seq_len(T)) {
    X <- LoadExprData(
      tissue_name = tissue_names[i],
      tissue_file_name = tissue_expr_file_names[i],
      sd_quantile = sd_quantile,
      max_genes_per_tissue = max_genes_per_tissue
    )
    expr_list[[i]] <- X
  }

  # ---- TS per tissue ----
  TS_per_tissue <- integer(T)
  TS_power_map  <- setNames(integer(T), tissue_names)
  TS_fit_curves <- vector("list", T)

  for (i in seq_len(T)) {
    ts_fit <- tryCatch(
      wgcna_pick_TS(
        expr_mat          = expr_list[[i]],
        powerVector       = powerVector,
        TOMType           = TOMType,
        cor_method        = cor_method,
        verbose           = verbose,
        targetR2          = targetR2,
        require_neg_slope = require_neg_slope
      ),
      error = function(e) { message("[TS] ", tissue_names[i], ": ", e$message); NULL }
    )
    if (!is.null(ts_fit)) {
      TS_per_tissue[i] <- ts_fit$beta
      fi <- ts_fit$fit_df
      fi$tissue <- tissue_names[i]
      TS_fit_curves[[i]] <- fi
    } else {
      TS_per_tissue[i] <- NA_integer_
      TS_fit_curves[[i]] <- data.frame(Power = powerVector, SFT.R.sq = NA_real_,
                                       tissue = tissue_names[i])
    }
    TS_power_map[tissue_names[i]] <- TS_per_tissue[i]
  }
  TS_power <- as.integer(stats::median(TS_per_tissue, na.rm = TRUE))
  TS_fit_curves_df <- do.call(rbind, TS_fit_curves)

  # ---- CT per pair (using the NEW wgcna_pick_CT) ----
  CT_betas <- integer(0)
  CT_power_map <- numeric(0)
  CT_fit_list <- list()

  if (T >= 2) {
    for (i in 1:(T - 1)) {
      for (j in (i + 1):T) {
        pair_id <- paste(tissue_names[i], tissue_names[j], sep = "||")
        ct_fit <- tryCatch(
          wgcna_pick_CT(
            expr_A = expr_list[[i]], expr_B = expr_list[[j]],
            aggregate_by_donor = aggregate_by_donor_CT,
            powerVector = powerVector,
            TOMType     = TOMType,
            cor_method  = cor_method,
            bicor_maxPOutliers = bicor_maxPOutliers,
            bicor_robustY      = bicor_robustY,
            targetR2    = targetR2,
            require_neg_slope = require_neg_slope,
            nBreaks     = ct_nBreaks,
            removeFirst = ct_removeFirst,
            min_common  = min_common_CT,
            verbose     = verbose
          ),
          error = function(e) { message("[CT] ", pair_id, ": ", e$message); NULL }
        )

        if (!is.null(ct_fit) && !is.na(ct_fit$beta)) {
          CT_betas <- c(CT_betas, ct_fit$beta)
          CT_power_map[pair_id] <- ct_fit$beta

          fi <- ct_fit$fitIndices %||% ct_fit$fit_df
          if (!is.null(fi)) {
            if (!"Rsquared.SFT" %in% names(fi) && "SFT.R.sq" %in% names(fi)) fi$Rsquared.SFT <- fi$SFT.R.sq
            if (!"slope.SFT"    %in% names(fi) && "slope"    %in% names(fi)) fi$slope.SFT    <- fi$slope
          }
          # שימור שם העמודה כדי להתאים לפונקציות שרטוט קיימות
          CT_fit_list[[length(CT_fit_list) + 1]] <- data.frame(
            pair = pair_id,
            Power = fi$Power,
            Rsquared.SFT = fi$SFT.R.sq,
            slope.SFT    = if ("slope" %in% names(fi)) fi$slope else fi$slope.SFT,
            mean.k       = fi$mean.k,
            median.k     = fi$median.k,
            max.k        = fi$max.k,
            stringsAsFactors = FALSE
          )
        } else {
          CT_fit_list[[length(CT_fit_list) + 1]] <-
            data.frame(pair = pair_id, Power = powerVector,
                       Rsquared.SFT = NA_real_, slope.SFT = NA_real_,
                       mean.k = NA_real_, median.k = NA_real_, max.k = NA_real_)
        }
      }
    }
  }

  CT_power <- if (length(CT_betas)) as.integer(stats::median(CT_betas, na.rm = TRUE)) else 3L
  CT_fit_curves_df <- if (length(CT_fit_list)) do.call(rbind, CT_fit_list) else
    data.frame(pair = character(), Power = numeric(), Rsquared.SFT = numeric())

  list(
    method         = "wgcna+scaleFreeFitIndex(CT)",
    TS_power       = TS_power,
    CT_power       = CT_power,
    TS_per_tissue  = TS_per_tissue,
    CT_per_pair    = as.integer(CT_betas),
    TS_fit_curves  = TS_fit_curves_df,
    CT_fit_curves  = CT_fit_curves_df,
    TS_power_map   = TS_power_map,
    CT_power_map   = CT_power_map
  )
}

# ---- TS: degrees from within-tissue adjacency ----
.compute_TS_degrees <- function(expr_mat, beta, cor_method = "pearson", unsigned = TRUE) {
  S <- cor(expr_mat, use = "pairwise.complete.obs", method = cor_method)
  if (unsigned) S <- abs(S)
  A <- S^beta
  diag(A) <- 0
  k <- rowSums(A, na.rm = TRUE)
  k[is.finite(k) & k > 0]
}

# ---- CT: degrees from cross-tissue bipartite adjacency ----
.compute_CT_degrees <- function(expr_A, expr_B, beta,
                                cor_method = "pearson",
                                unsigned = TRUE,
                                aggregate_by_donor = FALSE,
                                min_common = 3L) {
  XA <- expr_A; XB <- expr_B
  if (aggregate_by_donor) {
    XA <- .aggregate_by_donor(XA)
    XB <- .aggregate_by_donor(XB)
  }
  common <- intersect(rownames(XA), rownames(XB))
  if (length(common) < min_common) return(numeric(0))
  Mi <- XA[common, , drop = FALSE]
  Mj <- XB[common, , drop = FALSE]

  S <- cor(Mi, Mj, use = "pairwise.complete.obs", method = cor_method)
  if (unsigned) S <- abs(S)
  S[!is.finite(S)] <- 0
  A <- S^beta

  # כמו ב-pickSoftThreshold_crossTissue שלך (נרמול ואז rescale)
  k <- c(rowSums(A) / ncol(A), colSums(A) / nrow(A))
  k <- k * min(nrow(A), ncol(A))
  k <- k[is.finite(k) & k > 0]
  k
}

plot_scaleFree_TS_CT <- function(
  beta_info,
  tissue_names, tissue_expr_file_names,
  mode = c("chosen", "all"),                 
  cor_method = "pearson",
  unsigned = TRUE,
  nBreaks = 10, removeFirst = TRUE,
  aggregate_by_donor_CT = FALSE,
  min_common_CT = 3L,
  out_dir = "scaleFreePlots"
){
  mode <- match.arg(mode)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)
  expr_list <- setNames(vector("list", T), tissue_names)
  for (i in seq_len(T)) {
    expr_list[[i]] <- LoadExprData(
      tissue_name = tissue_names[i],
      tissue_file_name = tissue_expr_file_names[i],
      sd_quantile = 0.00, max_genes_per_tissue = 1e9
    )
  }

  # ================= TS =================
  if (!is.null(beta_info$TS_fit_curves) && nrow(beta_info$TS_fit_curves)) {
    ts_map <- beta_info$TS_power_map
    for (t in tissue_names) {
      if (mode == "chosen") {
        betas <- as.numeric(ts_map[[t]])
        betas <- betas[is.finite(betas)]
      } else {
        betas <- unique(beta_info$TS_fit_curves$Power[beta_info$TS_fit_curves$tissue == t])
      }
      if (!length(betas)) next

      pdf(file.path(out_dir, sprintf("scaleFree_TS_%s.pdf", gsub("[^A-Za-z0-9]+","_", t))),
          width = 5, height = 5 * length(betas))
      par(mfrow = c(length(betas), 1), mar = c(5,5,3,2))
      for (b in betas) {
        k <- .compute_TS_degrees(expr_list[[t]], b, cor_method, unsigned)
        if (!length(k)) { plot.new(); title(sprintf("%s | TS | beta=%s [no k]", t, b)); next }
        WGCNA::scaleFreePlot(k, main = sprintf("%s | TS | \u03B2 = %s", t, b),
                             nBreaks = nBreaks, removeFirst = removeFirst)
        fit <- WGCNA::scaleFreeFitIndex(k, nBreaks = nBreaks, removeFirst = removeFirst)
        mtext(sprintf("R^2 = %.3f | slope = %.3f", fit$Rsquared.SFT, fit$slope.SFT),
              side = 3, adj = 1, cex = 0.8)
      }
      dev.off()
    }
  }

  # ================= CT =================
  if (!is.null(beta_info$CT_fit_curves) && nrow(beta_info$CT_fit_curves)) {
    ct_map <- beta_info$CT_power_map
    pairs <- combn(tissue_names, 2, simplify = FALSE)
    for (pr in pairs) {
      a <- pr[1]; b <- pr[2]
      pair_id_ab <- paste(a, b, sep = "||")
      pair_id_ba <- paste(b, a, sep = "||")

      if (mode == "chosen" && length(ct_map)) {
        chosen <- if (!is.na(ct_map[pair_id_ab])) ct_map[pair_id_ab] else ct_map[pair_id_ba]
        betas <- as.numeric(chosen)
        betas <- betas[is.finite(betas)]
        if (!length(betas)) next
      } else if (mode == "chosen") {
        betas <- as.numeric(beta_info$CT_power)
      } else {
        keep <- beta_info$CT_fit_curves$pair %in% c(pair_id_ab, pair_id_ba)
        betas <- unique(beta_info$CT_fit_curves$Power[keep])
      }
      if (!length(betas)) next

      pdf(file.path(out_dir, sprintf("scaleFree_CT_%s__%s.pdf",
                                     gsub("[^A-Za-z0-9]+","_", a),
                                     gsub("[^A-Za-z0-9]+","_", b))),
          width = 5, height = 5 * length(betas))
      par(mfrow = c(length(betas), 1), mar = c(5,5,3,2))
      for (beta in betas) {
        k <- .compute_CT_degrees(expr_list[[a]], expr_list[[b]], beta,
                                 cor_method, unsigned,
                                 aggregate_by_donor = aggregate_by_donor_CT,
                                 min_common = min_common_CT)
        if (!length(k)) { plot.new(); title(sprintf("%s↔%s | CT | beta=%s [no k]", a, b, beta)); next }
        WGCNA::scaleFreePlot(k, main = sprintf("%s \u2194 %s | CT | \u03B2 = %s", a, b, beta),
                             nBreaks = nBreaks, removeFirst = removeFirst)
        fit <- WGCNA::scaleFreeFitIndex(k, nBreaks = nBreaks, removeFirst = removeFirst)
        mtext(sprintf("R^2 = %.3f | slope = %.3f", fit$Rsquared.SFT, fit$slope.SFT),
              side = 3, adj = 1, cex = 0.8)
      }
      dev.off()
    }
  }

  message("✓ scaleFreePlot PDFs ", normalizePath(out_dir))
}




# =========================
# Integration hook into your main function
# =========================
# In XWGCNA_Clusters_autoBeta(), add parameters:
#   beta_method = c("custom", "wgcna"),
#   wgcna_powerVector = seq(0.5, 20, length.out = 20),
#   aggregate_by_donor_CT = FALSE,
# and replace the auto-beta block with:
#
#   if (auto_beta) {
#     if (match.arg(beta_method) == "wgcna") {
#       message("Auto-picking TS/CT betas (WGCNA)…")
#       beta_info <- wgcna_auto_pick_powers(
#         tissue_names, tissue_expr_file_names,
#         sd_quantile = sd_quantile,
#         max_genes_per_tissue = max_genes_per_tissue,
#         TOMType = TOMType,
#         cor_method = cor_method,
#         powerVector = wgcna_powerVector,
#         targetR2 = targetR2,
#         aggregate_by_donor_CT = aggregate_by_donor_CT
#       )
#     } else {
#       message("Auto-picking TS/CT betas (custom)…")
#       beta_info <- auto_pick_powers(
#         tissue_names, tissue_expr_file_names,
#         sd_quantile = sd_quantile,
#         max_genes_per_tissue = max_genes_per_tissue,
#         cor_method = cor_method,
#         beta_grid = beta_grid,
#         targetR2 = targetR2
#       )
#     }
#     TS_power <- beta_info$TS_power
#     TS_map   <- beta_info$TS_power_map
#     CT_power <- beta_info$CT_power
#     CT_map   <- beta_info$CT_power_map
#   } else {
#     TS_map <- setNames(rep.int(TS_power, length(tissue_names)), tissue_names)
#     CT_map <- .make_CT_map(tissue_names, CT_power)
#   }
#
# Your downstream plotting function `plot_beta_curves_per_tissue()` will
# work as-is on `beta_info` from either method.




# =============================== kME analysis ===============================
# ========= Helpers =========
.safe_cor <- function(x, y, method = "pearson", min_common = 3L) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_common) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], use = "pairwise.complete.obs", method = method))
}

.pc1 <- function(X) {
  # X: donors x genes (center+scale → PC1 scores, length = donors)
  X <- as.matrix(X)
  if (ncol(X) == 1L) {
    sx <- scale(X[,1], center = TRUE, scale = TRUE)
    return(as.numeric(sx))
  }
  pr <- suppressWarnings(prcomp(X, center = TRUE, scale. = TRUE))
  as.numeric(pr$x[,1])
}

.make_gene_id <- function(tissue, gene_symbol) paste0(tissue, "_", gene_symbol)

load_donor_mats <- function(tissue_names, tissue_expr_file_names,
                            sd_quantile = 0.00, max_genes_per_tissue = 5000) {
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)
  donor_mats <- vector("list", T); names(donor_mats) <- tissue_names
  for (i in seq_len(T)) {
    X <- LoadExprData(tissue_names[i], tissue_expr_file_names[i],
                      sd_quantile = sd_quantile,
                      max_genes_per_tissue = max_genes_per_tissue)
    donor_mats[[i]] <- .aggregate_by_donor(X)  
  }
  donor_mats
}

build_gene_metadata <- function(clusters_table) {
  ct <- as.data.frame(clusters_table, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(ct) <- sub("\\.", " ", colnames(ct))
  req <- c("Cluster ID", "Tissue", "Gene Symbol")
  stopifnot(all(req %in% colnames(ct)))
  ct$`Cluster ID` <- as.integer(ct$`Cluster ID`)
  ct$module_id    <- paste0("M", ct$`Cluster ID`)
  ct$gene_id      <- .make_gene_id(ct$Tissue, ct$`Gene Symbol`)
  ct[, c("gene_id", "Tissue", "Gene Symbol", "module_id", "Cluster ID")]
}

compute_module_eigengenes_per_tissue <- function(gene_meta, donor_mats,
                                                 min_genes_for_ME = 1L) {
  tissues <- names(donor_mats)
  modules <- sort(unique(gene_meta$module_id))
  MEs_by_tissue <- lapply(tissues, function(tt) {
    donors_tt <- rownames(donor_mats[[tt]])
    out <- matrix(NA_real_, nrow = length(donors_tt), ncol = 0,
                  dimnames = list(donors_tt, NULL))
    out
  })
  names(MEs_by_tissue) <- tissues

  for (m in modules) {
    gm <- gene_meta[gene_meta$module_id == m, , drop = FALSE]
    for (tt in tissues) {
      ids_tt <- gm$gene_id[gm$Tissue == tt]
      ids_tt <- intersect(ids_tt, colnames(donor_mats[[tt]]))
      if (length(ids_tt) >= min_genes_for_ME) {
        X <- donor_mats[[tt]][, ids_tt, drop = FALSE]
        me <- .pc1(X)
        cur <- MEs_by_tissue[[tt]]
        MEs_by_tissue[[tt]] <- cbind(cur, setNames(data.frame(me), m))
        colnames(MEs_by_tissue[[tt]])[ncol(MEs_by_tissue[[tt]])] <- m
      }
    }
  }
  MEs_by_tissue <- lapply(MEs_by_tissue, function(M) {
    if (is.null(colnames(M)) || ncol(M) == 0L) {
      data.frame(row.names = rownames(M))
    } else {
      as.data.frame(M, check.names = FALSE)
    }
  })
  MEs_by_tissue
}

compute_kME_tables <- function(gene_meta, donor_mats, MEs_by_tissue,
                               cor_method = "pearson",
                               unsigned_kME = TRUE,
                               min_common = 3L) {
  tissues <- names(donor_mats)
  modules <- sort(unique(gene_meta$module_id))
  genes   <- gene_meta$gene_id
  kME_max  <- matrix(NA_real_, nrow = length(genes), ncol = length(modules),
                     dimnames = list(genes, modules))
  kME_self <- matrix(NA_real_, nrow = length(genes), ncol = length(modules),
                     dimnames = list(genes, modules))
  best_ME_tissue <- matrix(NA_character_, nrow = length(genes), ncol = length(modules),
                           dimnames = list(genes, modules))

  gene2tissue <- setNames(gene_meta$Tissue, gene_meta$gene_id)

  for (g in genes) {
    tg <- gene2tissue[[g]]
    vg <- donor_mats[[tg]][, g, drop = TRUE]
    donors_g <- rownames(donor_mats[[tg]])

    for (m in modules) {
      if (!is.null(MEs_by_tissue[[tg]]) && (m %in% colnames(MEs_by_tissue[[tg]]))) {
        me_tg   <- MEs_by_tissue[[tg]][, m, drop = TRUE]
        donors_c <- intersect(donors_g, rownames(MEs_by_tissue[[tg]]))
        k_self <- .safe_cor(vg[donors_c], me_tg[donors_c],
                            method = cor_method, min_common = min_common)
        if (unsigned_kME) k_self <- abs(k_self)
        kME_self[g, m] <- k_self
      }

      best_val <- NA_real_; best_tt <- NA_character_
      for (tt in tissues) {
        if (is.null(MEs_by_tissue[[tt]]) || !(m %in% colnames(MEs_by_tissue[[tt]]))) next
        me_tt   <- MEs_by_tissue[[tt]][, m, drop = TRUE]
        donors_tt <- rownames(MEs_by_tissue[[tt]])
        donors_c  <- intersect(donors_g, donors_tt)
        r <- .safe_cor(vg[donors_c], me_tt[donors_c],
                       method = cor_method, min_common = min_common)
        if (is.na(r)) next
        v <- if (unsigned_kME) abs(r) else r
        if (is.na(best_val) || v > best_val) {
          best_val <- v; best_tt <- tt
        }
      }
      kME_max[g, m] <- best_val
      best_ME_tissue[g, m] <- best_tt
    }
  }

  list(kME_max = kME_max,
       kME_self = kME_self,
       best_ME_tissue = best_ME_tissue)
}

classify_core_bridge <- function(gene_meta, kME_max, best_ME_tissue,
                                 alpha = 2, tau = 0.5,
                                 S_core = 0.70, H_core = 0.40, 
                                 b0_bridge = 0.35                
){
  stopifnot(identical(rownames(kME_max), gene_meta$gene_id))
  modules <- colnames(kME_max)
  mod_self <- setNames(gene_meta$module_id, gene_meta$gene_id)
  tissue_self <- setNames(gene_meta$Tissue, gene_meta$gene_id)

  genes <- rownames(kME_max)
  k_self   <- mapply(function(g, m) kME_max[g, m], g = genes, m = mod_self[genes])
  k_other_max <- apply(kME_max, 1, function(v) {
    mself <- mod_self[names(v)[1]] 
    v[names(which.max(v*0 + 1))] 
  })
  k_other_max <- vapply(genes, function(g) {
    mself <- mod_self[[g]]
    v <- kME_max[g, setdiff(modules, mself), drop = TRUE]
    if (!length(v)) return(NA_real_)
    max(v, na.rm = TRUE)
  }, numeric(1))

  delta <- k_self - k_other_max
  W <- kME_max^alpha
  W <- sweep(W, 1, rowSums(W, na.rm = TRUE) + 1e-12, "/")
  S  <- apply(W, 1, max, na.rm = TRUE)
  H  <- -rowSums(W * log(pmax(W, 1e-12)), na.rm = TRUE) / log(ncol(W))
  n_above_tau <- rowSums(kME_max >= tau, na.rm = TRUE)

  second_idx <- apply(kME_max, 1, function(v){
    o <- order(v, decreasing = TRUE); if (length(o) >= 2) o[2] else NA_integer_
  })
  second_mod <- ifelse(is.na(second_idx), NA_character_, colnames(kME_max)[second_idx])
  second_val <- ifelse(is.na(second_idx), NA_real_, kME_max[cbind(seq_along(second_idx), second_idx)])
  second_tt  <- ifelse(is.na(second_idx), NA_character_,
                       best_ME_tissue[cbind(seq_along(second_idx), second_idx)])
  cross_tissue_bridge <- second_tt != tissue_self[genes]

  class <- ifelse(!is.na(k_self) &
                    (k_self >= max(S_core, 0)) &
                    (delta >= 0.30 | S >= S_core) &
                    (H <= H_core) &
                    (n_above_tau <= 2),
                  "core",
                  ifelse(!is.na(second_val) &
                           (second_val >= b0_bridge | n_above_tau >= 3),
                         "bridge", "other"))

  out <- data.frame(
    gene_id = genes,
    tissue_self = tissue_self[genes],
    module_self = mod_self[genes],
    kME_self = k_self,
    kME_second = second_val,
    module_second = second_mod,
    ME_tissue_second = second_tt,
    delta = delta,
    S = S,
    H = H,
    n_above_tau = n_above_tau,
    cross_tissue_bridge = cross_tissue_bridge,
    class = class,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  rownames(out) <- NULL
  out
}

export_kME_and_classes <- function(kME_list, classes_df, out_prefix = "xwgcna") {
  kME_max  <- kME_list$kME_max
  kME_self <- kME_list$kME_self
  best_tt  <- kME_list$best_ME_tissue

  kME_max_df <- data.frame(gene_id = rownames(kME_max), kME_max, check.names = FALSE)
  write.table(kME_max_df, file = paste0(out_prefix, "_kME_gene_by_module.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  kME_self_df <- data.frame(gene_id = rownames(kME_self), kME_self, check.names = FALSE)
  write.table(kME_self_df, file = paste0(out_prefix, "_kME_selfTissue.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  best_tt_df <- data.frame(gene_id = rownames(best_tt), best_tt, check.names = FALSE)
  write.table(best_tt_df, file = paste0(out_prefix, "_kME_best_ME_tissue.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  write.table(classes_df, file = paste0(out_prefix, "_gene_core_bridge.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(list(
    kME_gene_by_module = paste0(out_prefix, "_kME_gene_by_module.tsv"),
    kME_selfTissue     = paste0(out_prefix, "_kME_selfTissue.tsv"),
    kME_best_ME_tissue = paste0(out_prefix, "_kME_best_ME_tissue.tsv"),
    gene_core_bridge   = paste0(out_prefix, "_gene_core_bridge.tsv")
  ))
}

run_kME_and_core_bridge <- function(
  tissue_names, tissue_expr_file_names,
  clusters_table,
  sd_quantile = 0.00, max_genes_per_tissue = 5000,
  cor_method = "pearson", unsigned_kME = TRUE,
  min_common = 3L,
  alpha = 2, tau = 0.5, S_core = 0.70, H_core = 0.40, b0_bridge = 0.35,
  out_prefix = "xwgcna"
){
  donor_mats <- load_donor_mats(tissue_names, tissue_expr_file_names,
                                sd_quantile = sd_quantile,
                                max_genes_per_tissue = max_genes_per_tissue)
  gene_meta  <- build_gene_metadata(clusters_table)
  MEs_by_tissue <- compute_module_eigengenes_per_tissue(gene_meta, donor_mats)

  kME_list <- compute_kME_tables(
    gene_meta = gene_meta,
    donor_mats = donor_mats,
    MEs_by_tissue = MEs_by_tissue,
    cor_method = cor_method,
    unsigned_kME = unsigned_kME,
    min_common = min_common
  )

  classes_df <- classify_core_bridge(
    gene_meta = gene_meta,
    kME_max = kME_list$kME_max,
    best_ME_tissue = kME_list$best_ME_tissue,
    alpha = alpha, tau = tau, S_core = S_core, H_core = H_core, b0_bridge = b0_bridge
  )

  export_kME_and_classes(kME_list, classes_df, out_prefix = out_prefix)
  list(kME = kME_list, classes = classes_df)
}


XWGCNA_Clusters_autoBeta <- function(
    tissue_names = NULL,
    tissue_expr_file_names = NULL,
    sd_quantile = 0.00,
    max_genes_per_tissue = 5000,
    TS_power = 6,
    CT_power = 3,
    cor_method = "pearson",
    TOMType = "unsigned",
    minClusterSize = 30,
    cluster_type_thr = 0.95,
    out_prefix = "xwgcna",
    save_intermediates = TRUE,
    plot_heatmap = TRUE,
    auto_beta = TRUE,
    beta_method = c("custom", "wgcna", "wgcna_new"),
    targetR2 = 0.80,
    beta_grid = c(1:10, seq(12, 20, 2)),
    wgcna_powerVector = seq(0.5, 20, length.out = 20),
    aggregate_by_donor_CT = FALSE,
    plot_beta_curves = TRUE,
    blockwise_TOM = TRUE,
    TOM_block_size = 2000L,
    group = "young",
    scaleFree_plots = "all",
    scaleFree_nBreaks = 12,
    scaleFree_removeFirst = TRUE
){
    stopifnot(length(tissue_names) == length(tissue_expr_file_names))

    TS_map <- NULL
    CT_map <- NULL
    beta_info <- NULL
    if (auto_beta) {
      if (beta_method == "wgcna_new"){
        message("Auto-picking TS/CT betas (WGCNA new)…")
        beta_info <- wgcna_auto_pick_powers_new(
          tissue_names, tissue_expr_file_names,
          sd_quantile = sd_quantile,
          max_genes_per_tissue = max_genes_per_tissue,
          TOMType = TOMType,
          cor_method = cor_method,
          powerVector = wgcna_powerVector,
          targetR2 = targetR2,
          require_neg_slope = TRUE,
          verbose = 5,
          aggregate_by_donor_CT = aggregate_by_donor_CT,
          min_common_CT = 3L,
          ct_nBreaks = 12
        )
      } else if (beta_method == "wgcna") {
        message("Auto-picking TS/CT betas (WGCNA)…")
        beta_info <- wgcna_auto_pick_powers(
          tissue_names, tissue_expr_file_names,
          sd_quantile = sd_quantile,
          max_genes_per_tissue = max_genes_per_tissue,
          TOMType = TOMType,
          cor_method = cor_method,
          powerVector = wgcna_powerVector,
          targetR2 = targetR2,
          aggregate_by_donor_CT = aggregate_by_donor_CT
        )
      } else {
        message("Auto-picking TS/CT betas ...")
        beta_info <- auto_pick_powers(
            tissue_names, tissue_expr_file_names,
            sd_quantile = sd_quantile,
            max_genes_per_tissue = max_genes_per_tissue,
            cor_method = cor_method,
            beta_grid = beta_grid,
            targetR2 = targetR2
        )
      }
        TS_power <- beta_info$TS_power
        TS_map <- beta_info$TS_power_map
        CT_power <- beta_info$CT_power
        CT_map <- beta_info$CT_power_map
        message(sprintf("Auto β selected: TS=%d (per tissue: %s); CT=%d",
                    TS_power, paste(beta_info$TS_per_tissue, collapse = ","), CT_power))
        if (plot_beta_curves) {
          plots_list <- plot_beta_curves_per_tissue(
            beta_info,
            out_prefix = paste0(out_prefix, "_per_tissue"),
            targetR2   = targetR2,
            save_png   = TRUE
            )
          tissues_vec <- tissue_names
          make_all_beta_plots(
            beta_info,
            tissues = tissues_vec,
            out_prefix = out_prefix
          )
        }
        if (scaleFree_plots != "none") {
            plot_scaleFree_TS_CT(
                beta_info = beta_info,
                tissue_names = tissue_names,
                tissue_expr_file_names = tissue_expr_file_names,
                mode = match.arg(scaleFree_plots),        
                out_dir = file.path("plots", paste0(out_prefix, "_scaleFree")),
                cor_method = cor_method,
                unsigned = (tolower(TOMType) == "unsigned"),
                nBreaks = scaleFree_nBreaks,
                removeFirst = scaleFree_removeFirst,
                aggregate_by_donor_CT = aggregate_by_donor_CT,
                min_common_CT = 3L
            )
        }
      } else {
      message("Using constant β values.")
      TS_map <- setNames(rep.int(TS_power, length(tissue_names)), tissue_names)
      CT_map <- .make_CT_map(tissue_names, CT_power)
      plot_beta_curves <- FALSE
    }
    message(sprintf(
        "Adjacency: %d tissues, up to %d genes/tissue (sd_quantile=%.2f).",
        length(tissue_names), max_genes_per_tissue, sd_quantile
    ))

    adj_mat <- AdjacencyFromExpr(
        tissue_names = tissue_names,
        tissue_expr_file_names = tissue_expr_file_names,
        sd_quantile = sd_quantile,
        max_genes_per_tissue = max_genes_per_tissue,
        cor_method = cor_method,
        TS_power_map = TS_map,
        CT_power_map = CT_map,
        default_TS = TS_power,
        default_CT = CT_power
    )
  
    if (save_intermediates) {
        saveRDS(adj_mat, file = paste0(out_prefix, "_adjacency.rds"))
    }
    gc()

    message("Computing TOM (this may take a while) ...")

    TOM_mat <- Cross_Tissue_TOM(adj_mat)

    if (save_intermediates) {
        saveRDS(TOM_mat, file = paste0(out_prefix, "_TOM.rds"))
    }
    rm(adj_mat)
    gc()

    message("Clustering on TOM...")
    clusters_table <- Clusters_Table(
        TOM_mat,
        minClusterSize = minClusterSize,
        plot_heatmap = plot_heatmap,
        tissue_names = tissue_names,
        tissue_expr_file_names = tissue_expr_file_names,
        group = group
    )

    if (!plot_heatmap) rm(TOM_mat); gc()

    clusters_details <- Clusters_Details(
        clusters_table,
        cluster_type_thr = cluster_type_thr
    )

    write.table(
        clusters_table,
        file = paste0(out_prefix, "_Cluster_table.txt"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )

    write.table(
        clusters_details,
        file = paste0(out_prefix, "_Cluster_details.txt"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )
    # res <- run_kME_and_core_bridge(
    #   tissue_names = tissue_names,
    #   tissue_expr_file_names = tissue_expr_file_names,
    #   clusters_table = clusters_table,       
    #   cor_method = "pearson", unsigned_kME = TRUE,
    #   out_prefix = "xwgcna_mycohort"
    # )
    clusters_table <- as.data.frame(clusters_table, stringsAsFactors = FALSE)
    clusters_details <- as.data.frame(clusters_details, stringsAsFactors = FALSE)

    if (auto_beta) {
        return(list(
            clusters_table   = clusters_table,
            clusters_details = clusters_details,
            TS_power         = TS_power,
            CT_power         = CT_power,
            beta_info        = beta_info
        ))
    } else {
        return(list(
            clusters_table   = clusters_table,
            clusters_details = clusters_details,
            TS_power         = TS_power,
            CT_power         = CT_power
        ))
    }
  }


