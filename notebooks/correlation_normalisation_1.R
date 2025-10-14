# ==== Compare three correlation methods for cross-tissue edges ====
# Methods:
#   1) RAW Pearson (pairwise.complete.obs)
#   2) Permutation Z — exact as proposed: build a null by permuting one side and compute Z = (C - mean(null)) / sd(null).
#      We can also return an r-equivalent at a reference n for comparability; fast default: analytic Z (Fisher); optional Monte Carlo permutations.
#   3) CorShrink-style empirical Bayes shrinkage on Fisher z (via ashr), works for rectangular (A×B)
#
# Also computes the donor-intersection size n_ij for each edge, and summarizes
# dependence between n_ij and |r| (and also |Z| if requested) across methods.
#
# This file plugs into your existing helpers you shared (LoadExprData, .aggregate_by_donor, etc.).
# Make sure to source your helper script first so those functions are available.
#
# Example usage (after editing paths):
#   src <- "your_helpers_and_this_file.R"; source(src)
#   res <- compare_correlation_methods(
#            tissue_names = c("T1","T2","T3"),
#            tissue_expr_file_names = c("T1.csv","T2.csv","T3.csv"),
#            sd_quantile = 0.00, max_genes_per_tissue = 1500,
#            perm_mode = "analytic", n_ref = "median", edge_sample = 200000
#          )
#   # View per-pair/per-method Spearman rho(|r|, n)
#   res$summary
#   # Plots are written under plots/ by default


# ---------- small utils ----------
.require_packages <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = TRUE, quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse=", "),
                         "\nInstall: install.packages(c(",
                         paste(sprintf('"%s"', miss), collapse=", "), "))")
}

.null_coalesce <- function(a, b) if (!is.null(a)) a else b

# Compute pairwise sample sizes (n_ij) under pairwise.complete.obs for a bipartite block
# Mi: donors×genes_A, Mj: donors×genes_B
.compute_n_mat <- function(Mi, Mj) {
  Ri <- is.finite(as.matrix(Mi))
  Rj <- is.finite(as.matrix(Mj))
  t(Ri) %*% Rj
}

# Efficient Monte Carlo null for bipartite block: permute rows of Mj B times.
# Returns Z = (r - mean_perm)/sd_perm. (For big matrices this can be heavy—use small B or sample cols.)
.permZ_bipartite_mc <- function(Mi, Mj, B = 50L, cor_method = "pearson") {
  pA <- ncol(Mi); pB <- ncol(Mj)
  S  <- stats::cor(Mi, Mj, use = "pairwise.complete.obs", method = cor_method)
  S[!is.finite(S)] <- 0
  sumM   <- matrix(0, pA, pB)
  sumSqM <- matrix(0, pA, pB)
  for (b in seq_len(B)) {
    idx <- sample(nrow(Mj))
    Mjb <- Mj[idx, , drop = FALSE]
    Sb  <- stats::cor(Mi, Mjb, use = "pairwise.complete.obs", method = cor_method)
    Sb[!is.finite(Sb)] <- 0
    sumM   <- sumM   + Sb
    sumSqM <- sumSqM + Sb*Sb
  }
  mu  <- sumM / B
  var <- pmax(sumSqM / B - mu*mu, 1e-12)
  sd  <- sqrt(var)
  Z   <- (S - mu) / sd
  Z[!is.finite(Z)] <- 0
  dimnames(Z) <- dimnames(S)
  list(S = S, Z = Z)
}

# Analytic Z under null (Fisher):  Z ≈ atanh(r) * sqrt(n-3)
# Then map to an effect-size-equivalent at reference sample size n_ref: r_eff = tanh(Z/sqrt(n_ref-3))
.permZ_bipartite_analytic <- function(Mi, Mj, n_ref = c("median","max", 30L), cor_method = "pearson"){
  S <- stats::cor(Mi, Mj, use = "pairwise.complete.obs", method = cor_method)
  S[!is.finite(S)] <- 0
  n_mat <- .compute_n_mat(Mi, Mj)
  Z <- atanh(pmax(pmin(S, 1-1e-9), -1+1e-9)) * sqrt(pmax(n_mat - 3, 0))
  n_vec <- as.vector(n_mat[n_mat >= 4])
  if (!length(n_vec)) {
    n0 <- 4L
  } else if (is.character(n_ref[1])) {
    n0 <- if (n_ref[1] == "max") max(n_vec) else stats::median(n_vec)
  } else n0 <- as.integer(n_ref[1])
  n0 <- max(as.integer(round(n0)), 4L)
  R_eff <- tanh(Z / sqrt(n0 - 3))
  dimnames(R_eff) <- dimnames(S)
  list(S = S, Z = Z, R_eff = R_eff, n_mat = n_mat, n_ref = n0)
}

# CorShrink-like shrinkage for rectangular blocks using ashr on Fisher-z.
# If CorShrink is installed, can switch to it; otherwise this reproduces its core idea.
.corshrink_bipartite <- function(S, n_mat, ash_method = "fdr") {
  .require_packages(c("ashr"))
  Z   <- atanh(pmax(pmin(S, 1-1e-9), -1+1e-9))
  se  <- 1/sqrt(pmax(n_mat - 3, 1))
  # vectorize
  z_vec  <- as.vector(Z)
  se_vec <- as.vector(se)
  keep <- is.finite(z_vec) & is.finite(se_vec) & se_vec > 0
  fit <- ashr::ash(betahat = z_vec[keep], sebetahat = se_vec[keep],
                   method = ash_method, mixcompdist = "normal")
  post <- rep(NA_real_, length(z_vec))
  post[keep] <- ashr::get_pm(fit)   # posterior mean of z
  Z_shr <- matrix(post, nrow = nrow(S), ncol = ncol(S))
  R_shr <- tanh(Z_shr)
  dimnames(R_shr) <- dimnames(S)
  R_shr[!is.finite(R_shr)] <- 0
  R_shr
}

# Tidy edge sampler to keep memory reasonable
.sample_edges_df <- function(S_list, n_mat, tissueA, tissueB, edge_sample = 200000L) {
  stopifnot(length(S_list) >= 1)
  pA <- nrow(n_mat); pB <- ncol(n_mat)
  # sample indices
  tot <- pA * pB
  take <- min(edge_sample, tot)
  idx <- if (take < tot) sample.int(tot, take) else seq_len(tot)
  ia <- ((idx - 1L) %% pA) + 1L
  jb <- ((idx - 1L) %/% pA) + 1L
  geneA <- rownames(n_mat)[ia]; geneB <- colnames(n_mat)[jb]
  n_ij  <- n_mat[cbind(ia, jb)]
  out <- data.frame(tissueA = tissueA, tissueB = tissueB,
                    geneA = geneA, geneB = geneB, n = as.integer(n_ij))
  for (nm in names(S_list)) {
    out[[nm]] <- S_list[[nm]][cbind(ia, jb)]
  }
  out
}

# Summarize dependence between n and |r| for a tidy edges df with columns: method values
.summarize_n_dependence <- function(edges_df, methods = c("raw_r","perm_reff","corshrink_r","perm_Z")) {
  keep <- edges_df$n >= 3 & is.finite(edges_df$n)
  edges_df <- edges_df[keep, , drop = FALSE]
  res <- lapply(methods, function(m) {
    x <- edges_df[[m]]
    ok <- is.finite(x)
    x <- abs(x[ok]); n <- edges_df$n[ok]
    if (!length(n)) return(data.frame(method = m, spearman = NA_real_, lin_slope = NA_real_))
    rho <- suppressWarnings(stats::cor(n, x, method = "spearman"))
    fit <- tryCatch(coef(stats::lm(x ~ n)), error = function(e) c(NA, NA))
    data.frame(method = m, spearman = as.numeric(rho), lin_slope = as.numeric(fit[2]))
  })
  do.call(rbind, res)
}

# Optional: quick plots
.plot_n_vs_r <- function(edges_df, out_prefix, methods = c("raw_r","perm_reff","corshrink_r","perm_Z")) {
  .require_packages(c("ggplot2"))
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  for (m in methods) {
    if (!m %in% names(edges_df)) next
    df <- edges_df[is.finite(edges_df[[m]]) & edges_df$n >= 3, c("n", m)]
    names(df) <- c("n", "val")
    ylab <- if (m == "perm_Z") "|Z|" else "|r|"
    p <- ggplot2::ggplot(df, ggplot2::aes(x = n, y = abs(val))) +
      ggplot2::geom_point(alpha = 0.2, size = 0.6) +
      ggplot2::geom_smooth(method = "loess", se = FALSE) +
      ggplot2::labs(title = paste0(ylab, " vs n — ", m), x = "n (common donors)", y = ylab) +
      ggplot2::theme_bw()
    fn <- file.path("plots", paste0(out_prefix, "_nvs_", m, ".png"))
    ggplot2::ggsave(fn, p, width = 6.4, height = 4.5, dpi = 150)
  }
}

# ---------- Main entry ----------
compare_correlation_methods <- function(
  tissue_names,
  tissue_expr_file_names,
  sd_quantile = 0.00,
  max_genes_per_tissue = 2000,
  cor_method = "pearson",
  # RAW adjacency options (if you later want them):
  TS_beta_raw = 3L, CT_beta_raw = 2L,
  # Permutation-Z choices:
  perm_mode = c("analytic","montecarlo"),
  perm_B = 40L,
  n_ref = "max",
  # NEW: choose how to apply n_ref — per pair (default) or one global constant across all pairs
  n_ref_mode = c("pairwise","global"),
  # If n_ref_mode=="global": optionally force a specific numeric n_ref for all pairs
  # (if NULL, we'll derive a global n_ref from the distribution of |common donors| across pairs
  #  using the selector in `n_ref`, e.g., median/max)
  n_ref_global = NULL,
  # CorShrink (ashr) options:
  ash_method = "fdr",
  # Sampling to keep RAM reasonable in summaries/plots
  edge_sample = 200000L,
  make_plots = TRUE,
  out_prefix = "corr_compare",
  # NEW: also report Z itself (|Z| vs n)
  report_Z = TRUE
){
  perm_mode <- match.arg(perm_mode)
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)

  message(sprintf("[compare_correlation_methods] Start (%d tissues, perm_mode=%s, plots=%s)", T, perm_mode, make_plots))

  # Load and donor-aggregate per your helpers
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

    # Decide on global n_ref if requested
  n_ref_mode <- match.arg(n_ref_mode)
  global_n0 <- NA_integer_
  # We'll compute global_n0 after we know donor intersections (donor_list is ready now)
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
      message(sprintf("  - Global n_ref will be derived from %d pairwise donor intersections", length(n_sizes)))
      message(sprintf("    (min=%d, max=%d, median=%d, mean=%.1f)", min(n_sizes), max(n_sizes), stats::median(n_sizes), mean(n_sizes)))
      if (!length(n_sizes)) {
        global_n0 <- 4L
      } else if (is.character(n_ref[1])) {
        global_n0 <- if (n_ref[1] == "max") max(n_sizes) else stats::median(n_sizes)
      } else {
        global_n0 <- as.integer(n_ref[1])
      }
      global_n0 <- max(as.integer(round(global_n0)), 4L)
    }
    message(sprintf("  - Using global n_ref = %d (mode=%s)", global_n0, n_ref_mode))
  } else {
    message(sprintf("  - Using pairwise n_ref mode (%s)", n_ref_mode))
  }

  all_edges <- list()
  all_summ  <- list()

  # Iterate cross-tissue pairs
  if (T >= 2) {
    for (i in 1:(T - 1)) {
      for (j in (i + 1):T) {
        ti <- tissue_names[i]; tj <- tissue_names[j]
        di <- rownames(donor_list[[i]]); dj <- rownames(donor_list[[j]])
        common <- intersect(di, dj)
        message(sprintf("  - Pair %s || %s: %d common donors", ti, tj, length(common)))
        if (length(common) < 3L) {
          message(sprintf("Skip %s||%s: too few common donors (%d)", ti, tj, length(common)))
          next
        }
        Mi <- donor_list[[i]][common, , drop = FALSE]
        Mj <- donor_list[[j]][common, , drop = FALSE]

        # RAW
        S_raw <- stats::cor(Mi, Mj, use = "pairwise.complete.obs", method = cor_method)
        S_raw[!is.finite(S_raw)] <- 0
        n_mat <- .compute_n_mat(Mi, Mj)

        # perm-Z
        perm_Z_mat <- NULL
        if (perm_mode == "analytic") {
          message("    > Computing analytic perm-Z")
          az <- .permZ_bipartite_analytic(Mi, Mj, n_ref = n_ref[1], cor_method = cor_method)
          S_perm_eff <- az$R_eff
          n0 <- az$n_ref
          perm_Z_mat <- az$Z
        } else {
          # Monte Carlo (consider subsetting genes for speed before calling this)
          message(sprintf("    > Monte Carlo perm-Z with B=%d", perm_B))
          pm <- .permZ_bipartite_mc(Mi, Mj, B = perm_B, cor_method = cor_method)
          # Map Z to r at reference n
          n_vec <- as.vector(n_mat[n_mat >= 4])
          n0 <- if (!length(n_vec)) 4L else if (is.character(n_ref[1])) {
            if (n_ref[1] == "max") max(n_vec) else stats::median(n_vec)
          } else as.integer(n_ref[1])
          n0 <- max(as.integer(round(n0)), 4L)
          S_perm_eff <- tanh(pm$Z / sqrt(n0 - 3))
          dimnames(S_perm_eff) <- dimnames(pm$Z)
          perm_Z_mat <- pm$Z
        }        # Override mapping to r at a global reference n if requested
        if (n_ref_mode == "global" && is.finite(global_n0)) {
          S_perm_eff <- tanh(perm_Z_mat / sqrt(global_n0 - 3))
          n0 <- global_n0
          message(sprintf("    > Re-mapped perm results to global n_ref=%d", global_n0))
        }

        # CorShrink via ashr
        message("    > Applying CorShrink shrinkage")
        S_cs <- .corshrink_bipartite(S_raw, n_mat, ash_method = ash_method)

        # Collect tidy edges (sampled)
        S_list <- list(raw_r = S_raw, perm_reff = S_perm_eff, corshrink_r = S_cs)
        if (isTRUE(report_Z) && !is.null(perm_Z_mat)) S_list$perm_Z <- perm_Z_mat

        edf <- .sample_edges_df(
          S_list = S_list,
          n_mat = n_mat, tissueA = ti, tissueB = tj, edge_sample = edge_sample
        )
        edf$pair <- paste(ti, tj, sep = "||")
        all_edges[[length(all_edges) + 1L]] <- edf

        meths <- c("raw_r","perm_reff","corshrink_r")
        if (isTRUE(report_Z)) meths <- c(meths, "perm_Z")
        summ <- .summarize_n_dependence(edf, methods = meths)
        summ$tissueA <- ti; summ$tissueB <- tj; summ$pair <- paste(ti, tj, sep = "||")
        summ$n_ref   <- n0
        all_summ[[length(all_summ) + 1L]] <- summ

        if (isTRUE(make_plots)) {
          .plot_n_vs_r(edf, out_prefix = paste0(out_prefix, "_", ti, "__", tj), methods = meths)
        }
      }
    }
  }

  edges_df <- if (length(all_edges)) do.call(rbind, all_edges) else data.frame()
  summary_df <- if (length(all_summ)) do.call(rbind, all_summ) else data.frame()

  # Write outputs
  dir.create("output", showWarnings = FALSE, recursive = TRUE)
  if (nrow(edges_df)) {
    message(sprintf("[compare_correlation_methods] Writing %d sampled edges to output/%s_edges_sample.csv", nrow(edges_df), out_prefix))
    utils::write.csv(edges_df, file = file.path("output", paste0(out_prefix, "_edges_sample.csv")), row.names = FALSE)
  } else {
    message("[compare_correlation_methods] No sampled edges to write")
  }
  if (nrow(summary_df)) {
    message(sprintf("[compare_correlation_methods] Writing %d summary rows to output/%s_n_dependence_summary.csv", nrow(summary_df), out_prefix))
    utils::write.csv(summary_df, file = file.path("output", paste0(out_prefix, "_n_dependence_summary.csv")), row.names = FALSE)
  } else {
    message("[compare_correlation_methods] No summary rows to write")
  }

  message("[compare_correlation_methods] Done")
  list(edges = edges_df, summary = summary_df)
}


# ---------- Optional: build RAW adjacency matrices using your betas (TS=3, CT=2) ----------
# If you want to compare not just correlations but also adjacency under the RAW method.
# For the Z-perm and CorShrink methods, you can reuse the same mapping |r|^beta on their r_eff / r_shr outputs.

build_raw_adjacency_beta_fixed <- function(
  tissue_names, tissue_expr_file_names,
  sd_quantile = 0.00, max_genes_per_tissue = 2000,
  cor_method = "pearson", TS_beta = 3L, CT_beta = 2L
){
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)

  expr_list   <- vector("list", T)
  donor_list  <- vector("list", T)
  idx <- integer(T+1); rc_names <- character(0)
  for (i in seq_len(T)) {
    X <- LoadExprData(tissue_names[i], tissue_expr_file_names[i],
                      sd_quantile = sd_quantile, max_genes_per_tissue = max_genes_per_tissue)
    expr_list[[i]]  <- X
    donor_list[[i]] <- .aggregate_by_donor(X)
    idx[i+1] <- idx[i] + ncol(X)
    rc_names <- c(rc_names, colnames(X))
  }

  A <- matrix(0, nrow = idx[T+1], ncol = idx[T+1], dimnames = list(rc_names, rc_names))

  # TS blocks (within)
  for (i in seq_len(T)) {
    rows <- (idx[i]+1):idx[i+1]
    Sii <- abs(stats::cor(expr_list[[i]], use = "pairwise.complete.obs", method = cor_method))
    diag(Sii) <- 0
    A[rows, rows] <- Sii ^ TS_beta
  }

  # CT blocks (between)
  if (T >= 2) {
    for (i in 1:(T-1)) {
      di <- rownames(donor_list[[i]])
      rows <- (idx[i]+1):idx[i+1]
      for (j in (i+1):T) {
        cols <- (idx[j]+1):idx[j+1]
        dj <- rownames(donor_list[[j]])
        common <- intersect(di, dj)
        if (length(common) < 3L) next
        Mi <- donor_list[[i]][common, , drop = FALSE]
        Mj <- donor_list[[j]][common, , drop = FALSE]
        Sij <- abs(stats::cor(Mi, Mj, use = "pairwise.complete.obs", method = cor_method))
        A[rows, cols] <- Sij ^ CT_beta
        A[cols, rows] <- t(Sij) ^ CT_beta
      }
    }
  }
  A
}


# ======== Informative plots: enhanced visuals & diagnostics ========
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

.plot_n_vs_r <- function(edges_df, out_prefix, methods = c("raw_r","perm_reff","corshrink_r"), pair_label = NULL) {
  .require_or_stop(c("ggplot2", "scales"))
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  n_range <- range(edges_df$n[is.finite(edges_df$n)], na.rm = TRUE)
  null_df <- .r_null_curve_df(n_range, probs = c(0.95, 0.99))

  for (m in methods) {
    if (!m %in% names(edges_df)) next
    df <- edges_df[is.finite(edges_df[[m]]) & edges_df$n >= 3, c("n", m)]
    names(df) <- c("n", "r")
    if (!nrow(df)) next
    df$absr <- abs(df$r)
    rho <- suppressWarnings(stats::cor(df$n, df$absr, method = "spearman"))
    slope <- tryCatch(as.numeric(coef(stats::lm(absr ~ n, data = df))[2]), error = function(e) NA_real_)
    npts <- nrow(df)
    big  <- npts > 50000

    p <- ggplot2::ggplot(df, ggplot2::aes(x = n, y = absr)) +
      (if (big) ggplot2::geom_bin2d(bins = 60) else ggplot2::geom_point(alpha = 0.15, size = 0.6)) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.7) +
      ggplot2::geom_line(data = null_df, ggplot2::aes(x = n, y = r, linetype = q), color = "red") +
      ggplot2::scale_linetype_manual(values = c("95%" = "dashed", "99%" = "dotdash"), name = "Null |r| quantile") +
      ggplot2::scale_x_continuous(labels = scales::comma) +
      ggplot2::labs(title = paste0(ifelse(is.null(pair_label), "", paste0("[", pair_label, "] ")), "|r| vs n — ", m),
                    x = "n (common donors)", y = "|r|",
                    caption = sprintf("Spearman=%.3f, slope=%.4f, N=%s", rho, slope, scales::comma(npts))) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom",
                     plot.caption = ggplot2::element_text(size = 8))

    fn <- file.path("plots", paste0(out_prefix, "_nvsr_", ifelse(is.null(pair_label), "all", gsub("[^A-Za-z0-9]+","_", pair_label)), "_", m, ".png"))
    ggplot2::ggsave(fn, p, width = 7.2, height = 5.0, dpi = 150)
  }
}

plot_all_pairs_faceted <- function(edges_df, out_prefix, methods = c("raw_r","perm_reff","corshrink_r")) {
  .require_or_stop(c("ggplot2", "tidyr", "dplyr", "scales"))
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  keep <- c("pair", "n", intersect(methods, names(edges_df)))
  df <- edges_df[, keep, drop = FALSE]
  df <- tidyr::pivot_longer(df, cols = tidyselect::any_of(intersect(methods, names(edges_df))),
                            names_to = "method", values_to = "r")
  df <- df[is.finite(df$r) & df$n >= 3, ]
  if (!nrow(df)) return(invisible(NULL))
  df$absr <- abs(df$r)
  big <- nrow(df) > 200000
  p <- ggplot2::ggplot(df, ggplot2::aes(n, absr)) +
    (if (big) ggplot2::geom_bin2d(bins = 60) else ggplot2::geom_point(alpha = 0.10, size = 0.4)) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.5, color = "black") +
    ggplot2::facet_grid(pair ~ method) +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::labs(title = "|r| vs n across pairs and methods",
                  x = "n (common donors)", y = "|r|") +
    ggplot2::theme_bw() + ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
  fn <- file.path("plots", paste0(out_prefix, "_faceted_pairs_methods.png"))
  hh <- 0.9 + 2.0 * max(1L, length(unique(df$pair)))
  ggplot2::ggsave(fn, p, width = 12, height = hh, dpi = 140)
}

plot_delta_vs_n <- function(edges_df, out_prefix) {
  .require_or_stop(c("ggplot2", "tidyr", "dplyr", "scales"))
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  if (!all(c("raw_r","corshrink_r") %in% names(edges_df))) return(invisible(NULL))
  df <- edges_df[, c("pair", "n", intersect(c("raw_r","perm_reff","corshrink_r"), names(edges_df))), drop = FALSE]
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
    ggplot2::labs(title = "Method-induced change vs n",
                  x = "n (common donors)", y = "Δ|r| relative to RAW") +
    ggplot2::theme_bw()
  fn <- file.path("plots", paste0(out_prefix, "_delta_vs_n.png"))
  hh <- 0.9 + 2.0 * max(1L, length(unique(dfl$pair)))
  ggplot2::ggsave(fn, p, width = 12, height = hh, dpi = 140)
}


# # Convenience wrapper: run comparison + all enhanced plots automatically
# compare_correlation_methods_with_plots <- function(
#   ..., out_prefix = "corr_compare"
# ){
#   res <- compare_correlation_methods(..., out_prefix = out_prefix, make_plots = TRUE)
#   if (!is.null(res$edges) && nrow(res$edges)) {
#     plot_all_pairs_faceted(res$edges, out_prefix)
#     plot_delta_vs_n(res$edges, out_prefix)
#   }
#   invisible(res)
# }


# I used 
# res <- compare_correlation_methods(
#   tissue_names = c("AC","MFBA9BA46","PCGBA23"),
#   tissue_expr_file_names = c("/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_AC.csv","/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_MFBA9BA46.csv","/Users/edeneldar/CoExpression_ReProduction/rosmap_expr_files/ROSMAP_fixed_PCGBA23.csv"),
#   perm_mode   = "analytic",
#   n_ref_mode  = "global",   # <-- חדש: נרמול גלובלי
#   n_ref       = "max",   # נבחר n0 גלובלי = חציון |donor∩| בין כל זוגות הרקמות
#   make_plots  = TRUE,
#   out_prefix  = "cmp_global"
# )
# > source("/Users/edeneldar/CoExpression_ReProduction/notebooks/correlation_normalisation_1.R")
# > source("/Users/edeneldar/CoExpression_ReProduction/old_scripts/ClusteringBuilding.R")
# ============================
# XWGCNA_Clusters_autoBeta(): הרצה מלאה עם auto-β
# ============================
XWGCNA_Clusters_autoBeta <- function(
  tissue_names,
  tissue_expr_file_names,
  # סינון בסיסי
  sd_quantile = 0.00,
  max_genes_per_tissue = 5000,
  # סוג רשת / קורלציה
  TOMType = "unsigned",
  cor_method = "pearson",
  # Auto-β
  auto_beta = TRUE,
  beta_method = c("custom","wgcna","wgcna_new"),
  targetR2 = 0.80,
  require_neg_slope = TRUE,
  # רשת CT על Fisher אופציונלי
  ct_fisher = TRUE,
  ct_fisher_scheme   = c("to_ref","lambda"),
  ct_fisher_Nref     = c("median","max"),
  ct_fisher_cap_at_1 = TRUE,
  ct_fisher_lambda   = 10,
  ct_min_common      = 3L,
  # רשת/מודולים
  minClusterSize = 30,
  # plotting + IO
  out_prefix = "xwgcna_wgcnaNew",
  wgcna_powerVector = seq(0.5, 20, length.out = 20),
  beta_grid = c(1:10, seq(12, 20, 2)),   # fallback ל-"custom"
  make_plots = TRUE
){
  beta_method <- match.arg(beta_method)
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  dir.create("output", showWarnings = FALSE, recursive = TRUE)

  message(sprintf("[XWGCNA] auto_beta=%s, method=%s, TOMType=%s, cor=%s",
                  auto_beta, beta_method, TOMType, cor_method))

  # === שלב 1: בחירת β (TS/CT) ===
  if (auto_beta) {
    if (beta_method == "wgcna_new") {
      message("[XWGCNA] Picking betas via wgcna_new (TS by pickSoftThreshold; CT ע\"י scaleFreeFitIndex על דרגות bipartite)")
      beta_info <- wgcna_auto_pick_powers_new(
        tissue_names, tissue_expr_file_names,
        sd_quantile = sd_quantile,
        max_genes_per_tissue = max_genes_per_tissue,
        TOMType = TOMType, cor_method = cor_method,
        powerVector = wgcna_powerVector,
        targetR2 = targetR2, require_neg_slope = require_neg_slope,
        # Fisher ב־CT כדי לנרמל ל-n_ref
        ct_fisher = ct_fisher,
        ct_fisher_scheme = ct_fisher_scheme[1],
        ct_fisher_Nref   = ct_fisher_Nref[1],
        ct_fisher_cap_at_1 = ct_fisher_cap_at_1,
        ct_fisher_lambda   = ct_fisher_lambda
      )
    } else if (beta_method == "wgcna") {
      message("[XWGCNA] Picking betas via vanilla wgcna (pickSoftThreshold על expression/paired)")
      beta_info <- wgcna_auto_pick_powers(
        tissue_names, tissue_expr_file_names,
        sd_quantile = sd_quantile,
        max_genes_per_tissue = max_genes_per_tissue,
        TOMType = TOMType, cor_method = cor_method,
        powerVector = wgcna_powerVector,
        targetR2 = targetR2, require_neg_slope = require_neg_slope
      )
    } else {
      message("[XWGCNA] Picking betas via custom (log-binning degrees)")
      beta_info <- auto_pick_powers(
        tissue_names, tissue_expr_file_names,
        sd_quantile = sd_quantile,
        max_genes_per_tissue = max_genes_per_tissue,
        cor_method = cor_method,
        beta_grid = beta_grid,
        targetR2  = targetR2
      )
    }
    TS_map <- beta_info$TS_power_map
    CT_map <- beta_info$CT_power_map
  } else {
    stop("כרגע הפונקציה מיועדת לעבודה עם auto_beta=TRUE. אם תרצי, נוסיף גם מצב ידני.")
  }

  # === שלב 2: בניית adjacency לפי הבטאות שנבחרו (כולל Fisher ב-CT אם נדרש) ===
  A <- AdjacencyFromExpr(
    tissue_names = tissue_names,
    tissue_expr_file_names = tissue_expr_file_names,
    sd_quantile = sd_quantile,
    max_genes_per_tissue = max_genes_per_tissue,
    cor_method = cor_method,
    TS_power_map = TS_map,
    CT_power_map = CT_map,
    default_TS = ifelse(is.finite(beta_info$TS_power), beta_info$TS_power, 6L),
    default_CT = ifelse(is.finite(beta_info$CT_power), beta_info$CT_power, 3L),
    ct_fisher = ct_fisher,
    ct_fisher_scheme = ct_fisher_scheme[1],
    ct_fisher_Nref   = ct_fisher_Nref[1],
    ct_fisher_cap_at_1 = ct_fisher_cap_at_1,
    ct_fisher_lambda   = ct_fisher_lambda,
    ct_min_common      = ct_min_common,
    ct_too_few_action  = "stop"
  )

  # === שלב 3: TOM, קלאסטרים, טבלאות ופלטים ===
  TOM <- Cross_Tissue_TOM(A)
  clusters_tbl <- Clusters_Table(
    TOM_mat = TOM,
    minClusterSize = minClusterSize,
    plot_heatmap = TRUE,
    tissue_names = tissue_names,
    tissue_expr_file_names = tissue_expr_file_names,
    group = out_prefix
  )
  write.csv(clusters_tbl, file = file.path("output", paste0(out_prefix, "_clusters_table.csv")), row.names = FALSE)

  # פרטי קלאסטרים
  details <- Clusters_Details(clusters_tbl, cluster_type_thr = 0.95)
  write.csv(details, file = file.path("output", paste0(out_prefix, "_clusters_details.csv")), row.names = FALSE)

  # === Plots (בטות + סקייל-פרי) ===
  if (isTRUE(make_plots)) {
    # עקומות R^2~β + Heatmap של בטות מהחציות
    make_all_beta_plots(beta_info, tissue_names, out_prefix = out_prefix, thr = targetR2)

    # סקייל-פרי אמיתי על דרגות עבור הבטאות הנבחרות (TS+CT)
    plot_scaleFree_TS_CT(
      beta_info,
      tissue_names, tissue_expr_file_names,
      mode = "chosen", cor_method = cor_method,
      aggregate_by_donor_CT = FALSE,
      min_common_CT = ct_min_common,
      out_dir = file.path("plots", paste0(out_prefix, "_scaleFree"))
    )
  }

  message("✓ outputs:\n  - output/", out_prefix, "_clusters_table.csv",
          "\n  - output/", out_prefix, "_clusters_details.csv",
          "\n  - plots/", out_prefix, "_beta_series.pdf (וגם summaries/heatmap)",
          "\n  - plots/", out_prefix, "_scaleFree/*.pdf",
          "\n  - TOM heatmap נוצר תחת השם TOM_heatmap_WGCNA_", out_prefix, ".png")

  invisible(list(beta_info = beta_info, adj = A, TOM = TOM,
                 clusters_table = clusters_tbl, cluster_details = details))
}
