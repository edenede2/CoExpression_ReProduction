# xwgcna_corshrink_addon.R
suppressPackageStartupMessages({
  requireNamespace("WGCNA", quietly = TRUE)
})

# ==== CorShrink availability (no auto-install mid-run) ====
.require_CorShrink <- function() {
  if (!requireNamespace("CorShrink", quietly = TRUE)) {
    stop("CorShrink is missing. Please install once:\n  install.packages('CorShrink')\n# ואם אין ב-CRAN שלך:\n# devtools::install_github('kkdey/CorShrink')")
  }
}

# ==== CorShrink on union of donors (NA-aware), extracting A×B block ====
.corshrink_cross_union <- function(
  Mi, Mj, nameA="A", nameB="B",
  block_B = NULL,           # if NULL and p small -> full run
  corshrink_max_p = 6000,   # threshold for "full" before blocking
  min_nonmiss = 4L,         # min non-NA per column
  min_sd = 1e-8             # min SD per column
){
  .require_CorShrink()
  genes_i <- colnames(Mi); genes_j <- colnames(Mj)
  donors_union <- sort(unique(c(rownames(Mi), rownames(Mj))))

  col_screen <- function(X) {
    if (!ncol(X)) return(X)
    nn <- colSums(is.finite(X))
    ss <- suppressWarnings(apply(X, 2, sd, na.rm = TRUE))
    keep <- (nn >= min_nonmiss) & is.finite(ss) & (ss > min_sd)
    X[, keep, drop = FALSE]
  }

  # build union-side matrices
  XA <- matrix(NA_real_, nrow = length(donors_union), ncol = length(genes_i),
               dimnames = list(donors_union, paste0(nameA, ":", genes_i)))
  if (length(rownames(Mi))) XA[rownames(Mi), ] <- as.matrix(Mi)
  XA <- col_screen(XA)
  genes_i <- sub(paste0("^", nameA, ":"), "", colnames(XA))

  # prepare result (rows = filtered A, cols = ORIGINAL B; B holes will be NA to be filled later)
  AB <- matrix(NA_real_, nrow = length(genes_i), ncol = length(genes_j),
               dimnames = list(genes_i, genes_j))

  run_corshrink <- function(X) {
    cs <- try(CorShrink::CorShrinkData(X), silent = TRUE)
    if (inherits(cs, "try-error")) return(NULL)
    S <- if (is.list(cs) && !is.null(cs$cor)) cs$cor else if (is.list(cs) && !is.null(cs$ash_cor)) cs$ash_cor else as.matrix(cs)
    if (!is.matrix(S)) return(NULL)
    if (anyNA(S)) return(NULL)
    S
  }

  do_full <- (length(genes_i) + length(genes_j)) <= corshrink_max_p && is.null(block_B)
  if (do_full) {
    XB <- matrix(NA_real_, nrow = length(donors_union), ncol = length(genes_j),
                 dimnames = list(donors_union, paste0(nameB, ":", genes_j)))
    if (length(rownames(Mj))) XB[rownames(Mj), ] <- as.matrix(Mj)
    XB <- col_screen(XB)

    Xfull <- cbind(XA, XB)
    message(sprintf("CorShrink (full): p=%d (|%s|=%d, |%s|=%d)", ncol(Xfull), nameA, ncol(XA), nameB, ncol(XB)))
    S <- run_corshrink(Xfull)
    if (is.null(S)) return(NULL)

    ixA <- grepl(paste0("^", nameA, ":"), colnames(S))
    ixB <- grepl(paste0("^", nameB, ":"), colnames(S))
    AB2 <- S[ixA, ixB, drop = FALSE]

    # tidy dimnames + clamp
    rownames(AB2) <- sub(paste0("^", nameA, ":"), "", rownames(AB2))
    colnames(AB2) <- sub(paste0("^", nameB, ":"), "", colnames(AB2))
    AB2 <- pmax(pmin(AB2, 1), -1)

    AB[rownames(AB2), colnames(AB2)] <- AB2
    return(AB)
  }

  # blockwise on B
  if (is.null(block_B)) block_B <- 800L
  nb <- if (length(genes_j)) ceiling(length(genes_j) / block_B) else 0L
  message(sprintf("CorShrink in blocks: |A|=%d, |B|=%d, block_B=%d", length(genes_i), length(genes_j), block_B))

  for (b in seq_len(nb)) {
    j1 <- (b - 1L) * block_B + 1L
    j2 <- min(b * block_B, length(genes_j))
    gj <- genes_j[j1:j2]

    XB <- matrix(NA_real_, nrow = length(donors_union), ncol = length(gj),
                 dimnames = list(donors_union, paste0(nameB, ":", gj)))
    if (length(rownames(Mj))) XB[rownames(Mj), ] <- as.matrix(Mj[, gj, drop = FALSE])
    XB <- col_screen(XB)
    if (!ncol(XB)) next

    X <- cbind(XA, XB)
    message(sprintf("  block %d/%d: p=%d", b, nb, ncol(X)))
    S <- run_corshrink(X)
    if (is.null(S)) next

    ixA <- grepl(paste0("^", nameA, ":"), colnames(S))
    ixB <- grepl(paste0("^", nameB, ":"), colnames(S))
    AB2 <- S[ixA, ixB, drop = FALSE]

    rownames(AB2) <- sub(paste0("^", nameA, ":"), "", rownames(AB2))
    colnames(AB2) <- sub(paste0("^", nameB, ":"), "", colnames(AB2))
    AB2 <- pmax(pmin(AB2, 1), -1)

    AB[rownames(AB2), colnames(AB2)] <- AB2
  }
  AB
}

# ==== naive cross-tissue correlation on common donors ====
.cor_cross_common <- function(Mi, Mj, cor_method="pearson", min_common=3L) {
  common <- intersect(rownames(Mi), rownames(Mj))
  if (length(common) < min_common) {
    return(matrix(0, nrow=ncol(Mi), ncol=ncol(Mj),
                  dimnames=list(colnames(Mi), colnames(Mj))))
  }
  S <- stats::cor(Mi[common, , drop=FALSE], Mj[common, , drop=FALSE],
                  use="pairwise.complete.obs", method=cor_method)
  S[!is.finite(S)] <- 0
  S <- pmax(pmin(S, 1), -1)
  dimnames(S) <- list(colnames(Mi), colnames(Mj))
  S
}

# ==== beta curve from correlation matrix ====
.CT_beta_curve_from_S <- function(S, powerVector=c(1:10, seq(12,20,2)),
                                  TOMType="unsigned", unsigned=TRUE,
                                  nBreaks=12, removeFirst=TRUE) {
  if (!is.matrix(S)) S <- as.matrix(S)
  S[!is.finite(S)] <- 0
  S <- pmax(pmin(S, 1), -1)
  if (unsigned) S <- abs(S)

  fit_df <- data.frame(Power=powerVector, SFT.R.sq=NA_real_, slope=NA_real_,
                       mean.k=NA_real_, median.k=NA_real_, max.k=NA_real_,
                       stringsAsFactors = FALSE)
  .adj_from_cor <- function(S, beta, TOMType) {
    tt <- tolower(TOMType)
    if (startsWith(tt, "signed")) ((S+1)/2)^beta else abs(S)^beta
  }
  for (ix in seq_along(powerVector)) {
    b <- powerVector[ix]
    A <- .adj_from_cor(S, b, TOMType)
    k <- c(rowSums(A, na.rm=TRUE), colSums(A, na.rm=TRUE))
    k <- k[is.finite(k) & k > 0]
    if (!length(k)) next
    sf <- WGCNA::scaleFreeFitIndex(k, nBreaks=nBreaks, removeFirst=removeFirst)
    fit_df$SFT.R.sq[ix] <- sf$Rsquared.SFT
    fit_df$slope[ix]    <- sf$slope.SFT
    fit_df$mean.k[ix]   <- mean(k); fit_df$median.k[ix] <- stats::median(k); fit_df$max.k[ix] <- max(k)
  }
  fit_df
}

# ==== before/after plot (robust to CorShrink failure) ====
plot_CT_beta_before_after <- function(
  Mi, Mj, pair_id,
  cor_method="pearson",
  powerVector=c(1:10, seq(12,20,2)),
  TOMType="unsigned", unsigned=TRUE,
  min_common=3L, out_dir="plots",
  nameA="A", nameB="B",
  nBreaks=12, removeFirst=TRUE
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  S_before <- .cor_cross_common(Mi, Mj, cor_method=cor_method, min_common=min_common)
  S_after  <- .corshrink_cross_union(Mi, Mj, nameA=nameA, nameB=nameB)
  if (is.null(S_after)) S_after <- S_before
  # fill any residual holes using S_before
  mask <- !is.finite(S_after)
  if (any(mask)) S_after[mask] <- S_before[mask]
  S_after <- pmax(pmin(S_after, 1), -1)

  df_before <- .CT_beta_curve_from_S(S_before, powerVector, TOMType, unsigned,
                                     nBreaks=nBreaks, removeFirst=removeFirst)
  df_after  <- .CT_beta_curve_from_S(S_after,  powerVector, TOMType, unsigned,
                                     nBreaks=nBreaks, removeFirst=removeFirst)

  fn <- file.path(out_dir, paste0("CT_beta_before_after__", gsub("[^A-Za-z0-9]+","_", pair_id), ".pdf"))
  pdf(fn, width=6, height=4); on.exit(dev.off(), add=TRUE)
  yr <- range(c(df_before$SFT.R.sq, df_after$SFT.R.sq), na.rm=TRUE); if (!all(is.finite(yr))) yr <- c(0,1)
  plot(df_before$Power, df_before$SFT.R.sq, type="o", pch=16, lwd=1.2, col="grey30",
       xlab=expression(beta), ylab=expression(R^2), main=paste0("CT beta-curves — ", pair_id),
       ylim=yr)
  lines(df_after$Power, df_after$SFT.R.sq, type="o", pch=17, lwd=1.4, col="steelblue")
  abline(h=0.8, lty=3, col="grey70")
  legend("bottomright", bty="n", pch=c(16,17), lwd=c(1.2,1.4),
         col=c("grey30","steelblue"),
         legend=c("Before (common donors, naive cor)", "After (CorShrink, all donors)"))
  invisible(fn)
}

# ==== AdjacencyFromExpr with CorShrink (no Fisher) ====
AdjacencyFromExpr <- function(
    tissue_names = NULL,
    tissue_expr_file_names = NULL,
    sd_quantile = 0.00,
    max_genes_per_tissue = 5000,
    cor_method = "pearson",
    TS_power_map = NULL,
    CT_power_map = NULL,
    default_TS = 6L,
    default_CT = 3L,
    ct_min_common = 3L,
    ct_too_few_action = c("stop","zeros"),
    signed = FALSE,
    # new flags:
    ct_use_all_donors = FALSE,
    ct_use_corshrink  = TRUE,
    plot_CT_beta_before_after = TRUE,
    powerVector_for_plots = c(1:10, seq(12,20,2)),
    TOMType_for_plots = "unsigned",
    beta_curves_dir = "plots",
    save_intermediates = FALSE,
    out_prefix = "xwgcna",
    corshrink_max_p = 6000,
    block_B = 800L
) {
  ct_too_few_action <- match.arg(ct_too_few_action)
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)

  expr_list   <- vector("list", T)
  donor_mat   <- vector("list", T)
  idx <- integer(T+1); rc_names <- character(0)

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

  A_final <- matrix(0, nrow = idx[T + 1], ncol = idx[T + 1], dimnames = list(rc_names, rc_names))
  A_pre  <- if (save_intermediates) matrix(0, nrow = nrow(A_final), ncol = ncol(A_final),
                                           dimnames = dimnames(A_final)) else NULL
  A_post <- if (save_intermediates) matrix(0, nrow = nrow(A_final), ncol = ncol(A_final),
                                           dimnames = dimnames(A_final)) else NULL

  # ---- TS (no CorShrink by default) ----
  for (i in 1:T) {
    rows <- (idx[i] + 1):idx[i + 1]
    pow_i <- if(!is.null(TS_power_map) && !is.na(TS_power_map[tissue_names[i]])) {
      as.integer(TS_power_map[tissue_names[i]])
    } else default_TS
    Sii <- stats::cor(expr_list[[i]], use = "pairwise.complete.obs", method = cor_method)
    if (!signed) Sii <- abs(Sii)
    Sii[!is.finite(Sii)] <- 0
    diag(Sii) <- 0
    Aii <- Sii ^ pow_i
    Aii <- pmin(pmax(Aii, 0), 1)
    A_final[rows, rows] <- Aii
    if (save_intermediates) { A_pre[rows, rows] <- Aii; A_post[rows, rows] <- Aii }
  }

  # ---- CT ----
  if (T >= 2) {
    for (i in 1:(T - 1)) {
      rows <- (idx[i] + 1):idx[i + 1]
      for (j in (i + 1):T) {
        cols <- (idx[j] + 1):idx[j + 1]
        pair_ij <- paste(tissue_names[i], tissue_names[j], sep = "||")
        pair_ji <- paste(tissue_names[j], tissue_names[i], sep = "||")
        pow_ij <- if(!is.null(CT_power_map)) {
          if(!is.na(CT_power_map[pair_ij])) as.integer(CT_power_map[pair_ij])
          else if(!is.na(CT_power_map[pair_ji])) as.integer(CT_power_map[pair_ji])
          else default_CT
        } else default_CT

        # naive (common donors) for "before" and as a safety net
        S_before <- .cor_cross_common(donor_mat[[i]], donor_mat[[j]],
                                      cor_method = cor_method, min_common = ct_min_common)

        if (isTRUE(ct_use_all_donors) && isTRUE(ct_use_corshrink)) {
          S_after <- .corshrink_cross_union(
            donor_mat[[i]], donor_mat[[j]],
            nameA = tissue_names[i], nameB = tissue_names[j],
            block_B = block_B, corshrink_max_p = corshrink_max_p,
            min_nonmiss = 4L, min_sd = 1e-8
          )

          if (is.null(S_after)) {
            message(sprintf("[CT %s] CorShrink returned NULL/NA → fallback to naive common-donors cor.", pair_ij))
            S_after <- S_before
          } else {
            # clamp & hole-filling
            S_after <- pmax(pmin(S_after, 1), -1)
            na_mask <- !is.finite(S_after)
            if (any(na_mask)) S_after[na_mask] <- S_before[na_mask]
          }

          # optional plots
          if (isTRUE(plot_CT_beta_before_after)) {
            try(
              plot_CT_beta_before_after(
                donor_mat[[i]], donor_mat[[j]], pair_id = pair_ij,
                cor_method = cor_method, powerVector = powerVector_for_plots,
                TOMType = TOMType_for_plots, unsigned = !signed,
                min_common = ct_min_common, out_dir = beta_curves_dir,
                nameA = tissue_names[i], nameB = tissue_names[j]
              ),
              silent = TRUE
            )
          }

          C_pre  <- (if (!signed) abs(S_before) else (S_before + 1)/2) ^ pow_ij
          C_post <- (if (!signed) abs(S_after)  else (S_after  + 1)/2) ^ pow_ij
          C_pre  <- pmin(pmax(C_pre,  0), 1)
          C_post <- pmin(pmax(C_post, 0), 1)

          A_final[rows, cols] <- C_post; A_final[cols, rows] <- t(C_post)
          if (save_intermediates) {
            A_pre [rows, cols] <- C_pre;  A_pre [cols, rows] <- t(C_pre)
            A_post[rows, cols] <- C_post; A_post[cols, rows] <- t(C_post)
          }
        } else {
          # default: naive
          C <- (if (!signed) abs(S_before) else (S_before + 1)/2) ^ pow_ij
          C <- pmin(pmax(C, 0), 1)
          A_final[rows, cols] <- C; A_final[cols, rows] <- t(C)
          if (save_intermediates) {
            A_pre [rows, cols] <- C; A_pre [cols, rows] <- t(C)
            A_post[rows, cols] <- C; A_post[cols, rows] <- t(C)
          }
        }
      }
    }
  }

  # ---- adjacency final sanitation ----
  A_final[!is.finite(A_final)] <- 0
  A_final <- pmin(pmax(A_final, 0), 1)
  diag(A_final) <- 0
  A_final <- (A_final + t(A_final)) / 2

  if (isTRUE(save_intermediates)) {
    dir.create("output", showWarnings = FALSE, recursive = TRUE)
    if (!is.null(A_pre))  saveRDS(A_pre,  file = file.path("output", paste0(out_prefix, "_adjacency_preNorm.rds")))
    if (!is.null(A_post)) saveRDS(A_post, file = file.path("output", paste0(out_prefix, "_adjacency_postNorm_CorShrink.rds")))
  }

  A_final
}
