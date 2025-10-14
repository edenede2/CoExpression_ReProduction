# === Add-on: scale-free beta fitting and visualizations, per normalization method ===
# This augments your *new* compare_correlation_methods() so it will also:
#  (1) pick scale-free β (WGCNA-style) for CT blocks for each normalization method (RAW / PERM-Ref / CorShrink)
#  (2) export tidy curves (R^2 vs β) and a beta summary per pair & method
#  (3) produce overlay plots and per-method scaleFreePlot PDFs at chosen β
#  (4) improve |r| vs n plots (handles constant-n-per-pair by faceting & binning and adding null curves)
#
# Drop this file **after** your helpers are sourced. It redeclares compare_correlation_methods()
# with new arguments but keeps all your existing behavior.


# ----------------- safety: tiny utils -----------------
.require_or_stop <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = TRUE, quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse=", "),
                         "\nInstall: install.packages(",
                         paste(sprintf('"%s"', miss), collapse=", "), ")")
}

# Choose the minimal β whose R^2≥target (optionally slope≤0); else argmax R^2
.wgcna_choose_beta_min_above <- function(fit_df, targetR2 = 0.80, require_neg_slope = TRUE) {
  if (!nrow(fit_df)) return(NA_integer_)
  r2col <- if ("SFT.R.sq" %in% names(fit_df)) "SFT.R.sq" else if ("Rsquared.SFT" %in% names(fit_df)) "Rsquared.SFT" else NA
  slcol <- if ("slope" %in% names(fit_df)) "slope" else if ("slope.SFT" %in% names(fit_df)) "slope.SFT" else NA
  if (is.na(r2col)) return(NA_integer_)
  ok <- is.finite(fit_df[[r2col]]) & fit_df[[r2col]] >= targetR2
  if (require_neg_slope && !is.na(slcol)) ok <- ok & is.finite(fit_df[[slcol]]) & fit_df[[slcol]] < 0
  if (any(ok)) return(as.integer(min(fit_df$Power[ok])))
  as.integer(fit_df$Power[which.max(fit_df[[r2col]])])
}

# degrees for bipartite adjacency given a correlation matrix S and β
.ct_degrees_from_S <- function(S, beta, TOMType = "unsigned") {
  tt <- tolower(TOMType)
  if (startsWith(tt, "signed")) {
    A <- ((S + 1)/2) ^ beta
  } else {
    A <- abs(S) ^ beta
  }
  # bipartite degrees: concat row and column degrees
  k <- c(rowSums(A, na.rm = TRUE), colSums(A, na.rm = TRUE))
  k[is.finite(k) & k > 0]
}

# Fit scale-free curves for a *bipartite* correlation matrix S across powerVector
.pickSoftThreshold_CT_from_S <- function(S,
                                         powerVector = c(1:10, seq(12, 30, 2)),
                                         TOMType = "unsigned",
                                         nBreaks = 50, removeFirst = TRUE) {
  .require_or_stop(c("WGCNA"))
  out <- data.frame(Power = powerVector, SFT.R.sq = NA_real_, slope = NA_real_,
                    mean.k = NA_real_, median.k = NA_real_, max.k = NA_real_)
  for (i in seq_along(powerVector)) {
    b <- powerVector[i]
    k <- .ct_degrees_from_S(S, b, TOMType)
    if (!length(k)) next
    sf <- WGCNA::scaleFreeFitIndex(k, nBreaks = nBreaks, removeFirst = removeFirst)
    out$SFT.R.sq[i] <- sf$Rsquared.SFT
    out$slope[i]    <- sf$slope.SFT
    out$mean.k[i]   <- mean(k); out$median.k[i] <- stats::median(k); out$max.k[i] <- max(k)
  }
  out
}

# Null |r| quantile curves as a function of n (for overlays)
.r_null_curve_df <- function(n_min_max, probs = c(0.95, 0.99)) {
  if (length(n_min_max) == 1) n_min_max <- c(n_min_max, n_min_max)
  nn <- seq.int(max(4L, floor(n_min_max[1])), max(4L, ceiling(n_min_max[2])), by = 1L)
  out <- do.call(rbind, lapply(probs, function(p) {
    data.frame(n = nn,
               r = tanh(stats::qnorm((1 + p)/2) / sqrt(pmax(nn - 3, 1))),
               q = paste0(round(p * 100), "%"))
  }))
  out
}

# Improved |r| vs n plot: binned, with null quantiles, robust to constant-n-per-pair
.plot_n_vs_r_improved <- function(edges_df, out_prefix, methods = c("raw_r","perm_reff","corshrink_r"),
                                  pair_label = NULL) {
  .require_or_stop(c("ggplot2","scales"))
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  n_range <- range(edges_df$n[is.finite(edges_df$n)], na.rm = TRUE)
  null_df <- .r_null_curve_df(n_range, probs = c(0.95, 0.99))
  for (m in methods) {
    if (!m %in% names(edges_df)) next
    df <- edges_df[is.finite(edges_df[[m]]) & edges_df$n >= 3, c("n", m)]
    names(df) <- c("n","r"); if (!nrow(df)) next
    df$absr <- abs(df$r); big <- nrow(df) > 150000
    rho <- suppressWarnings(stats::cor(df$n, df$absr, method = "spearman"))
    slope <- tryCatch(as.numeric(coef(stats::lm(absr ~ n, data = df))[2]), error = function(e) NA_real_)
    p <- ggplot2::ggplot(df, ggplot2::aes(n, absr)) +
      (if (big) ggplot2::geom_bin2d(bins = 60) else ggplot2::geom_point(alpha = 0.15, size = 0.5)) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.7) +
      ggplot2::geom_line(data = null_df, ggplot2::aes(n, r, linetype = q), color = "red") +
      ggplot2::scale_linetype_manual(values = c("95%" = "dashed", "99%" = "dotdash"), name = "Null |r| quantile") +
      ggplot2::scale_x_continuous(labels = scales::comma) +
      ggplot2::labs(title = paste0(ifelse(is.null(pair_label), "", paste0("[", pair_label, "] ")), "|r| vs n — ", m),
                    x = "n (common donors)", y = "|r|",
                    caption = sprintf("Spearman=%.3f, slope=%.4f, N=%s", rho, slope, scales::comma(nrow(df)))) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(size = 8))
    fn <- file.path("plots", paste0(out_prefix, "_nvsr_", ifelse(is.null(pair_label), "all", gsub("[^A-Za-z0-9]+","_", pair_label)), "_", m, ".png"))
    ggplot2::ggsave(fn, p, width = 7.2, height = 5, dpi = 150)
  }
}

# Faceted across pairs & methods (handles constant-n per pair)
plot_all_pairs_faceted <- function(edges_df, out_prefix, methods = c("raw_r","perm_reff","corshrink_r")) {
  .require_or_stop(c("ggplot2","tidyr","dplyr","scales"))
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  keep <- c("pair","n", intersect(methods, names(edges_df)))
  df <- edges_df[, keep, drop = FALSE]
  df <- tidyr::pivot_longer(df, cols = tidyselect::any_of(intersect(methods, names(edges_df))),
                            names_to = "method", values_to = "r")
  df <- df[is.finite(df$r) & df$n >= 3, ]; if (!nrow(df)) return(invisible(NULL))
  df$absr <- abs(df$r); big <- nrow(df) > 200000
  p <- ggplot2::ggplot(df, ggplot2::aes(n, absr)) +
    (if (big) ggplot2::geom_bin2d(bins = 60) else ggplot2::geom_point(alpha = 0.10, size = 0.4)) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.5, color = "black") +
    ggplot2::facet_grid(pair ~ method) +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::labs(title = "|r| vs n across pairs and methods", x = "n (common donors)", y = "|r|") +
    ggplot2::theme_bw() + ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
  fn <- file.path("plots", paste0(out_prefix, "_faceted_pairs_methods.png"))
  hh <- 0.9 + 2.0 * max(1L, length(unique(df$pair)))
  ggplot2::ggsave(fn, p, width = 12, height = hh, dpi = 140)
}

# Δ|r| vs n curves
plot_delta_vs_n <- function(edges_df, out_prefix) {
  .require_or_stop(c("ggplot2","tidyr","dplyr","scales"))
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  if (!all(c("raw_r","corshrink_r") %in% names(edges_df))) return(invisible(NULL))
  df <- edges_df[, c("pair","n", intersect(c("raw_r","perm_reff","corshrink_r"), names(edges_df))), drop = FALSE]
  df <- df[stats::complete.cases(df[, c("n","raw_r")]), ]
  df$delta_cs_raw   <- abs(df$corshrink_r) - abs(df$raw_r)
  if ("perm_reff" %in% names(df)) df$delta_perm_raw <- abs(df$perm_reff) - abs(df$raw_r)
  dfl <- tidyr::pivot_longer(df, cols = tidyselect::any_of(c("delta_cs_raw","delta_perm_raw")),
                             names_to = "delta", values_to = "value")
  dfl <- dfl[is.finite(dfl$value) & dfl$n >= 3, ]
  if (!nrow(dfl)) return(invisible(NULL))
  p <- ggplot2::ggplot(dfl, ggplot2::aes(n, value)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_smooth(method = "loess", se = FALSE) +
    ggplot2::facet_grid(pair ~ delta, scales = "free_y") +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::labs(title = "Method-induced change vs n", x = "n (common donors)", y = "Δ|r| relative to RAW") +
    ggplot2::theme_bw()
  fn <- file.path("plots", paste0(out_prefix, "_delta_vs_n.png"))
  hh <- 0.9 + 2.0 * max(1L, length(unique(dfl$pair)))
  ggplot2::ggsave(fn, p, width = 12, height = hh, dpi = 140)
}

# Overlay R^2(β) curves per pair across methods + chosen β markers
.plot_SFT_overlay_per_pair <- function(curves_df, beta_summary_df, out_prefix) {
  .require_or_stop(c("ggplot2"))
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  pairs <- unique(curves_df$pair)
  for (pr in pairs) {
    sub <- curves_df[curves_df$pair == pr, , drop = FALSE]
    if (!nrow(sub)) next
    p <- ggplot2::ggplot(sub, ggplot2::aes(Power, Rsquared.SFT, color = method)) +
      ggplot2::geom_line(linewidth = 1) + ggplot2::geom_point(size = 1) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(title = paste0("Scale-free fit (R^2) vs β — ", pr), x = "β", y = "R^2") +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
    # annotate chosen β
    pick <- beta_summary_df[beta_summary_df$pair == pr, , drop = FALSE]
    if (nrow(pick)) {
      for (m in unique(pick$method)) {
        b <- pick$beta[pick$method == m]; if (!length(b) || !is.finite(b)) next
        p <- p + ggplot2::geom_vline(xintercept = b, linetype = "dashed")
      }
    }
    fn <- file.path("plots", paste0(out_prefix, "_SFT_overlay_", gsub("[^A-Za-z0-9]+","_", pr), ".png"))
    ggplot2::ggsave(fn, p, width = 7.2, height = 5.2, dpi = 150)
  }
}

# Per-method scaleFreePlot PDFs at *chosen* β (one PDF per pair×method)
.save_scaleFreePlot_PDFs <- function(S_list_by_pair, beta_summary_df, out_dir = "scaleFreePlots",
                                     TOMType = "unsigned", nBreaks = 50, removeFirst = TRUE) {
  .require_or_stop(c("WGCNA"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(S_list_by_pair)) {
    Slist <- S_list_by_pair[[nm]]  # list RAW / PERM / CS
    for (m in names(Slist)) {
      S <- Slist[[m]]; if (is.null(S)) next
      b <- beta_summary_df$beta[beta_summary_df$pair == nm & beta_summary_df$method == m]
      if (!length(b) || !is.finite(b)) next
      k <- .ct_degrees_from_S(S, b[1], TOMType)
      if (!length(k)) next
      fn <- file.path(out_dir, paste0("scaleFree_CT_", gsub("[^A-Za-z0-9]+","_", nm), "_", m, "_beta", b[1], ".pdf"))
      grDevices::pdf(fn, width = 5, height = 5)
      WGCNA::scaleFreePlot(k, main = paste0(nm, " | ", m, " | β=", b[1]), nBreaks = nBreaks, removeFirst = removeFirst)
      fit <- WGCNA::scaleFreeFitIndex(k, nBreaks = nBreaks, removeFirst = removeFirst)
      mtext(sprintf("R^2=%.3f | slope=%.3f", fit$Rsquared.SFT, fit$slope.SFT), side = 3, adj = 1, cex = 0.8)
      grDevices::dev.off()
    }
  }
}


# ----------------- MAIN: compare_correlation_methods() with SFT add-ons -----------------
compare_correlation_methods <- function(
  tissue_names,
  tissue_expr_file_names,
  sd_quantile = 0.00,
  max_genes_per_tissue = 2000,
  cor_method = "pearson",
  # RAW adjacency options (kept for compatibility)
  TS_beta_raw = 3L, CT_beta_raw = 2L,
  # Permutation-Z choices:
  perm_mode = c("analytic","montecarlo"),
  perm_B = 40L,
  n_ref = "max",
  # NEW: how to apply n_ref — per pair (default) or one global constant across all pairs
  n_ref_mode = c("pairwise","global"),
  n_ref_global = NULL,
  # CorShrink (ashr) options:
  ash_method = "fdr",
  # Sampling to keep RAM reasonable in summaries/plots
  edge_sample = 200000L,
  make_plots = TRUE,
  out_prefix = "corr_compare",
  report_Z = TRUE,
  # ================= NEW: scale-free beta & plots =================
  fit_scale_free = TRUE,
  sft_power_grid = c(1:10, seq(12, 30, 2)),
  sft_target_R2 = 0.80,
  sft_require_neg_slope = TRUE,
  sft_TOMType = "unsigned",
  sft_nBreaks = 50, sft_removeFirst = TRUE,
  sft_make_plots = TRUE
){
  perm_mode <- match.arg(perm_mode)
  n_ref_mode <- match.arg(n_ref_mode)
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)

  message(sprintf("[compare_correlation_methods] Start (%d tissues, perm_mode=%s, plots=%s)", T, perm_mode, make_plots))

  # ---- Load + donor-aggregate ----
  expr_list   <- vector("list", T)
  donor_list  <- vector("list", T)
  for (i in seq_len(T)) {
    message(sprintf("  - Loading tissue %d/%d: %s", i, T, tissue_names[i]))
    X <- LoadExprData(
      tissue_name = tissue_names[i],
      tissue_file_name = tissue_expr_file_names[i],
      sd_quantile = sd_quantile,
      max_genes_per_tissue = max_genes_per_tissue
    )
    expr_list[[i]]  <- X
    donor_list[[i]] <- .aggregate_by_donor(X)
    message(sprintf("    Completed loading %s (cells=%d, genes=%d)", tissue_names[i], nrow(X), ncol(X)))
  }

  # ---- Global n_ref (optional) ----
  global_n0 <- NA_integer_
  if (n_ref_mode == "global") {
    if (!is.null(n_ref_global)) {
      global_n0 <- max(as.integer(round(n_ref_global)), 4L)
    } else {
      n_sizes <- integer(0)
      if (T >= 2) {
        for (ii in 1:(T - 1)) {
          di <- rownames(donor_list[[ii]])
          for (jj in (ii + 1):T) {
            dj <- rownames(donor_list[[jj]])
            n_sizes <- c(n_sizes, length(intersect(di, dj)))
          }
        }
      }
      message(sprintf("  - Global n_ref derived from %d intersections", length(n_sizes)))
      if (!length(n_sizes)) {
        global_n0 <- 4L
      } else if (is.character(n_ref[1])) {
        global_n0 <- if (n_ref[1] == "max") max(n_sizes) else stats::median(n_sizes)
      } else {
        global_n0 <- as.integer(n_ref[1])
      }
      global_n0 <- max(as.integer(round(global_n0)), 4L)
    }
    message(sprintf("  - Using global n_ref = %d", global_n0))
  } else {
    message("  - Using pairwise n_ref mode")
  }

  all_edges <- list(); all_summ <- list()
  sft_curves <- list(); sft_beta_rows <- list()
  S_cache_for_plots <- list()  # per pair: list(RAW=, PERM=, CS=)

  if (T >= 2) {
    for (i in 1:(T - 1)) {
      for (j in (i + 1):T) {
        ti <- tissue_names[i]; tj <- tissue_names[j]
        di <- rownames(donor_list[[i]]); dj <- rownames(donor_list[[j]])
        common <- intersect(di, dj)
        pair_id <- paste(ti, tj, sep = "||")
        message(sprintf("  - Pair %s: %d common donors", pair_id, length(common)))
        if (length(common) < 3L) { message("    skip: too few common donors"); next }

        Mi <- donor_list[[i]][common, , drop = FALSE]
        Mj <- donor_list[[j]][common, , drop = FALSE]

        # ---- raw ----
        S_raw <- stats::cor(Mi, Mj, use = "pairwise.complete.obs", method = cor_method)
        S_raw[!is.finite(S_raw)] <- 0
        n_mat <- t(is.finite(as.matrix(Mi))) %*% is.finite(as.matrix(Mj))

        # ---- perm ----
        perm_Z_mat <- NULL; S_perm_eff <- NULL; n0_pair <- NA_integer_
        if (perm_mode == "analytic") {
          Z <- atanh(pmax(pmin(S_raw, 1-1e-9), -1+1e-9)) * sqrt(pmax(n_mat - 3, 0))
          if (n_ref_mode == "global" && is.finite(global_n0)) {
            n0_pair <- global_n0
          } else {
            n_vec <- as.vector(n_mat[n_mat >= 4])
            if (!length(n_vec)) n0_pair <- 4L else if (is.character(n_ref[1])) n0_pair <- if (n_ref[1]=="max") max(n_vec) else stats::median(n_vec) else n0_pair <- as.integer(n_ref[1])
            n0_pair <- max(as.integer(round(n0_pair)), 4L)
          }
          S_perm_eff <- tanh(Z / sqrt(n0_pair - 3)); dimnames(S_perm_eff) <- dimnames(S_raw)
          perm_Z_mat <- Z
        } else {
          # lightweight MC (optional path)
          B <- perm_B; pA <- ncol(Mi); pB <- ncol(Mj)
          sumM <- matrix(0, pA, pB); sumSq <- matrix(0, pA, pB)
          for (b in seq_len(B)) {
            idx <- sample(nrow(Mj)); Sb <- stats::cor(Mi, Mj[idx, , drop=FALSE], use="pairwise.complete.obs", method = cor_method)
            Sb[!is.finite(Sb)] <- 0; sumM <- sumM + Sb; sumSq <- sumSq + Sb*Sb
          }
          mu <- sumM/B; var <- pmax(sumSq/B - mu*mu, 1e-12); sd <- sqrt(var)
          Z <- (S_raw - mu)/sd; perm_Z_mat <- Z
          if (n_ref_mode == "global" && is.finite(global_n0)) n0_pair <- global_n0 else {
            n_vec <- as.vector(n_mat[n_mat >= 4]); n0_pair <- if (!length(n_vec)) 4L else if (is.character(n_ref[1])) if (n_ref[1]=="max") max(n_vec) else stats::median(n_vec) else as.integer(n_ref[1])
            n0_pair <- max(as.integer(round(n0_pair)), 4L)
          }
          S_perm_eff <- tanh(Z / sqrt(n0_pair - 3)); dimnames(S_perm_eff) <- dimnames(S_raw)
        }

        # ---- CorShrink via ashr ----
        .require_or_stop(c("ashr"))
        Z0 <- atanh(pmax(pmin(S_raw, 1-1e-9), -1+1e-9))
        se <- 1/sqrt(pmax(n_mat - 3, 1)); z_vec <- as.vector(Z0); se_vec <- as.vector(se)
        keep <- is.finite(z_vec) & is.finite(se_vec) & se_vec > 0
        fit <- ashr::ash(betahat = z_vec[keep], sebetahat = se_vec[keep], method = ash_method, mixcompdist = "normal")
        post <- rep(NA_real_, length(z_vec)); post[keep] <- ashr::get_pm(fit)
        Z_shr <- matrix(post, nrow = nrow(S_raw), ncol = ncol(S_raw)); S_cs <- tanh(Z_shr); S_cs[!is.finite(S_cs)] <- 0
        dimnames(S_cs) <- dimnames(S_raw)

        # ---- Collect sampled edges for |r|-vs-n ----
        # sample edges to keep memory small
        pA <- nrow(n_mat); pB <- ncol(n_mat); tot <- pA * pB; take <- min(edge_sample, tot)
        idx <- if (take < tot) sample.int(tot, take) else seq_len(tot)
        ia <- ((idx - 1L) %% pA) + 1L; jb <- ((idx - 1L) %/% pA) + 1L
        geneA <- rownames(n_mat)[ia]; geneB <- colnames(n_mat)[jb]
        n_ij  <- n_mat[cbind(ia, jb)]
        edf <- data.frame(tissueA = ti, tissueB = tj, geneA = geneA, geneB = geneB, n = as.integer(n_ij))
        edf$raw_r <- S_raw[cbind(ia, jb)]; edf$perm_reff <- S_perm_eff[cbind(ia, jb)]; edf$corshrink_r <- S_cs[cbind(ia, jb)]
        edf$pair <- pair_id
        all_edges[[length(all_edges) + 1L]] <- edf

        meths <- c("raw_r","perm_reff","corshrink_r")
        tmp <- lapply(meths, function(m) {
          x <- edf[[m]]; ok <- is.finite(x); x <- abs(x[ok]); n <- edf$n[ok]
          rho <- if (!length(n)) NA_real_ else suppressWarnings(stats::cor(n, x, method="spearman"))
          slope <- tryCatch(as.numeric(coef(stats::lm(x ~ n))[2]), error = function(e) NA_real_)
          data.frame(method = m, spearman = rho, lin_slope = slope)
        })
        summ <- do.call(rbind, tmp); summ$tissueA <- ti; summ$tissueB <- tj; summ$pair <- pair_id; summ$n_ref <- n0_pair
        all_summ[[length(all_summ) + 1L]] <- summ

        # ---- SFT fits per method (CT only) ----
        if (isTRUE(fit_scale_free)) {
          curves_raw  <- .pickSoftThreshold_CT_from_S(S_raw,       powerVector = sft_power_grid, TOMType = sft_TOMType,
                                                      nBreaks = sft_nBreaks, removeFirst = sft_removeFirst)
          curves_perm <- .pickSoftThreshold_CT_from_S(S_perm_eff, powerVector = sft_power_grid, TOMType = sft_TOMType,
                                                      nBreaks = sft_nBreaks, removeFirst = sft_removeFirst)
          curves_cs   <- .pickSoftThreshold_CT_from_S(S_cs,        powerVector = sft_power_grid, TOMType = sft_TOMType,
                                                      nBreaks = sft_nBreaks, removeFirst = sft_removeFirst)
          curves_raw$method <- "RAW"; curves_perm$method <- "PERM"; curves_cs$method <- "CorShrink"
          curves_raw$pair <- pair_id; curves_perm$pair <- pair_id; curves_cs$pair <- pair_id
          sft_curves[[length(sft_curves) + 1L]] <- rbind(curves_raw, curves_perm, curves_cs)

          b_raw  <- .wgcna_choose_beta_min_above(curves_raw,  targetR2 = sft_target_R2, require_neg_slope = sft_require_neg_slope)
          b_perm <- .wgcna_choose_beta_min_above(curves_perm, targetR2 = sft_target_R2, require_neg_slope = sft_require_neg_slope)
          b_cs   <- .wgcna_choose_beta_min_above(curves_cs,   targetR2 = sft_target_R2, require_neg_slope = sft_require_neg_slope)
          sft_beta_rows[[length(sft_beta_rows) + 1L]] <- data.frame(pair = pair_id, method = c("RAW","PERM","CorShrink"),
                                                                     beta = c(b_raw, b_perm, b_cs))
          # cache for PDFs
          S_cache_for_plots[[pair_id]] <- list(RAW = S_raw, PERM = S_perm_eff, CorShrink = S_cs)
        }

        # ---- Per-pair plots ----
        if (isTRUE(make_plots)) {
          .plot_n_vs_r_improved(edf, out_prefix = paste0(out_prefix, "_", ti, "__", tj),
                                methods = c("raw_r","perm_reff","corshrink_r"), pair_label = pair_id)
        }
      }
    }
  }

  edges_df <- if (length(all_edges)) do.call(rbind, all_edges) else data.frame()
  summary_df <- if (length(all_summ)) do.call(rbind, all_summ) else data.frame()

  # ---- write standard outputs ----
  dir.create("output", showWarnings = FALSE, recursive = TRUE)
  if (nrow(edges_df)) utils::write.csv(edges_df, file = file.path("output", paste0(out_prefix, "_edges_sample.csv")), row.names = FALSE)
  if (nrow(summary_df)) utils::write.csv(summary_df, file = file.path("output", paste0(out_prefix, "_n_dependence_summary.csv")), row.names = FALSE)

  # ---- SFT exports & plots ----
  sft_curves_df <- data.frame(); beta_summary_df <- data.frame()
  if (isTRUE(fit_scale_free) && length(sft_curves)) {
    sft_curves_df <- do.call(rbind, sft_curves)
    names(sft_curves_df)[names(sft_curves_df)=="SFT.R.sq"] <- "Rsquared.SFT"
    utils::write.csv(sft_curves_df, file = file.path("output", paste0(out_prefix, "_SFT_curves.csv")), row.names = FALSE)

    beta_summary_df <- do.call(rbind, sft_beta_rows)
    # also add R^2 at chosen β
    beta_summary_df$R2_at_beta <- NA_real_;
    for (ii in seq_len(nrow(beta_summary_df))) {
      pr <- beta_summary_df$pair[ii]; m  <- beta_summary_df$method[ii]; b <- beta_summary_df$beta[ii]
      if (!is.finite(b)) next
      sub <- sft_curves_df[sft_curves_df$pair == pr & sft_curves_df$method == m & sft_curves_df$Power == b, , drop = FALSE]
      if (nrow(sub)) beta_summary_df$R2_at_beta[ii] <- sub$Rsquared.SFT[1]
    }
    utils::write.csv(beta_summary_df, file = file.path("output", paste0(out_prefix, "_SFT_beta_summary.csv")), row.names = FALSE)

    if (isTRUE(sft_make_plots)) {
      .plot_SFT_overlay_per_pair(sft_curves_df, beta_summary_df, out_prefix = out_prefix)
      # Δβ vs RAW bar/dot plot
      .require_or_stop(c("ggplot2","dplyr"))
      wide <- tidyr::pivot_wider(beta_summary_df, id_cols = pair, names_from = method, values_from = beta)
      if (all(c("RAW","PERM","CorShrink") %in% names(wide))) {
        # במקום: wide <- wide %>% dplyr::mutate(...)
        wide$Delta_PERM      <- wide$PERM      - wide$RAW
        wide$Delta_CorShrink <- wide$CorShrink - wide$RAW

        p <- ggplot2::ggplot(wide, ggplot2::aes(x = pair)) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
          ggplot2::geom_point(ggplot2::aes(y = Delta_PERM), position = ggplot2::position_dodge(width = 0.4)) +
          ggplot2::geom_point(ggplot2::aes(y = Delta_CorShrink), shape = 17, position = ggplot2::position_dodge(width = 0.4)) +
          ggplot2::coord_flip() + ggplot2::labs(title = "Δβ vs RAW per pair", y = "Δβ", x = "Pair") + ggplot2::theme_bw()
        ggplot2::ggsave(file.path("plots", paste0(out_prefix, "_SFT_beta_deltas.png")), p,
                        width = 7, height = 4 + 0.3*length(unique(wide$pair)), dpi = 150)
      }
      # Per-method PDFs of scaleFreePlot at chosen β
      .save_scaleFreePlot_PDFs(S_cache_for_plots, beta_summary_df, out_dir = file.path("plots","scaleFreePlots"),
                               TOMType = sft_TOMType, nBreaks = sft_nBreaks, removeFirst = sft_removeFirst)
    }
  }

  # global multi-pair visuals
  if (isTRUE(make_plots) && nrow(edges_df)) {
    plot_all_pairs_faceted(edges_df, out_prefix)
    plot_delta_vs_n(edges_df, out_prefix)
  }

  message("[compare_correlation_methods] Done")
  list(edges = edges_df,
       summary = summary_df,
       sft_curves = sft_curves_df,
       sft_beta_summary = beta_summary_df)
}


# ================= Example call =================
# res <- compare_correlation_methods(
#   tissue_names = c("AC","MFBA9BA46","PCGBA23"),
#   tissue_expr_file_names = c(
#     "/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_AC.csv",
#     "/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_MFBA9BA46.csv",
#     "/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_PCGBA23.csv"),
#   perm_mode   = "analytic",
#   n_ref_mode  = "global",
#   n_ref       = "max",
#   make_plots  = TRUE,
#   out_prefix  = "cmp_global",
#   fit_scale_free = TRUE,
#   sft_target_R2 = 0.85,
#   sft_power_grid = c(1:10, seq(12, 30, 2))
# )
# > res <- compare_correlation_methods(  tissue_names = c("AC","MFBA9BA46","PCGBA23"),  tissue_expr_file_names = c(    "/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_AC.csv",    "/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_MFBA9BA46.csv",    "/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_PCGBA23.csv"),  perm_mode   = "analytic",  n_ref_mode  = "global",  n_ref       = "max",  make_plots  = TRUE,  out_prefix  = "cmp_global_5000",max_genes_per_tissue=5000,  fit_scale_free = TRUE,  sft_target_R2 = 0.8, sft_nBreaks = 12,  sft_power_grid = c(1:10, seq(12, 30, 2)))
# > res <- compare_correlation_methods(  tissue_names = c("AC","MFBA9BA46","PCGBA23"),  tissue_expr_file_names = c(    "/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_AC.csv",    "/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_MFBA9BA46.csv",    "/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_PCGBA23.csv"),  perm_mode   = "analytic",  n_ref_mode  = "global",  n_ref       = "max",  make_plots  = TRUE,  out_prefix  = "cmp_global_5000",max_genes_per_tissue=5000,  fit_scale_free = TRUE,  sft_target_R2 = 0.8, sft_nBreaks = 12,  sft_power_grid = c(1:10, seq(12, 30, 2)))

# Outputs:
#   output/cmp_global_edges_sample.csv
#   output/cmp_global_n_dependence_summary.csv
#   output/cmp_global_SFT_curves.csv
#   output/cmp_global_SFT_beta_summary.csv
#   plots/cmp_global_* (overlay R^2 vs β per pair, Δβ vs RAW, faceted |r|-vs-n, Δ|r|-vs-n, and scaleFreePlot PDFs)