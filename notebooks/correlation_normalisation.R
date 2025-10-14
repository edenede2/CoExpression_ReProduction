# ==== Compare three correlation methods for cross-tissue edges ====
# Methods:
#   1) RAW Pearson (pairwise.complete.obs)
#   2) Permutation Z  → mapped back to an "effect-size-equivalent" r via Fisher equalization
#      (fast default: analytic Z using Fisher; optional Monte Carlo permutations for smaller runs)
#   3) CorShrink-style empirical Bayes shrinkage on Fisher z (via ashr), works for rectangular (A×B)
#
# Also computes the donor-intersection size n_ij for each edge, and summarizes
# dependence between n_ij and |r| across methods.
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
.summarize_n_dependence <- function(edges_df, methods = c("raw_r","perm_reff","corshrink_r")) {
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
.plot_n_vs_r <- function(edges_df, out_prefix, methods = c("raw_r","perm_reff","corshrink_r")) {
  .require_packages(c("ggplot2"))
  dir.create("plots", showWarnings = FALSE, recursive = TRUE)
  for (m in methods) {
    df <- edges_df[is.finite(edges_df[[m]]) & edges_df$n >= 3, c("n", m)]
    names(df) <- c("n", "r")
    p <- ggplot2::ggplot(df, ggplot2::aes(x = n, y = abs(r))) +
      ggplot2::geom_point(alpha = 0.2, size = 0.6) +
      ggplot2::geom_smooth(method = "loess", se = FALSE) +
      ggplot2::labs(title = paste0("|r| vs n — ", m), x = "n (common donors)", y = "|r|") +
      ggplot2::theme_bw()
    fn <- file.path("plots", paste0(out_prefix, "_nvsr_", m, ".png"))
    ggplot2::ggsave(fn, p, width = 6.4, height = 4.5, dpi = 150)
    message(sprintf("[.plot_n_vs_r] Saved plot %s", fn))
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
  n_ref = c("median","max", 30L),
  # CorShrink (ashr) options:
  ash_method = "fdr",
  # Sampling to keep RAM reasonable in summaries/plots
  edge_sample = 200000L,
  make_plots = TRUE,
  out_prefix = "corr_compare"
){
  perm_mode <- match.arg(perm_mode)
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)
  message(sprintf("[compare_correlation_methods] Starting run for %d tissues (perm_mode=%s)", T, perm_mode))

  # Load and donor-aggregate per your helpers
  expr_list   <- vector("list", T)
  donor_list  <- vector("list", T)
  for (i in seq_len(T)) {
    message(sprintf("[compare_correlation_methods] Loading tissue %s from %s", tissue_names[i], tissue_expr_file_names[i]))
    X <- LoadExprData(
      tissue_name = tissue_names[i],
      tissue_file_name = tissue_expr_file_names[i],
      sd_quantile = sd_quantile,
      max_genes_per_tissue = max_genes_per_tissue
    )
    expr_list[[i]]  <- X
    donor_list[[i]] <- .aggregate_by_donor(X)
    message(sprintf("[compare_correlation_methods] Completed donor aggregation for %s (genes=%d, donors=%d)", tissue_names[i], ncol(X), nrow(donor_list[[i]])))
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
        if (length(common) < 3L) {
          message(sprintf("Skip %s||%s: too few common donors (%d)", ti, tj, length(common)))
          next
        }
        message(sprintf("[compare_correlation_methods] Processing pair %s||%s with %d common donors", ti, tj, length(common)))
        Mi <- donor_list[[i]][common, , drop = FALSE]
        Mj <- donor_list[[j]][common, , drop = FALSE]

        # RAW
        S_raw <- stats::cor(Mi, Mj, use = "pairwise.complete.obs", method = cor_method)
        S_raw[!is.finite(S_raw)] <- 0
        n_mat <- .compute_n_mat(Mi, Mj)
        message(sprintf("[compare_correlation_methods] Computed raw correlations for %s||%s", ti, tj))

        # perm-Z
        if (perm_mode == "analytic") {
          az <- .permZ_bipartite_analytic(Mi, Mj, n_ref = n_ref[1], cor_method = cor_method)
          S_perm_eff <- az$R_eff
          n0 <- az$n_ref
          message(sprintf("[compare_correlation_methods] Analytic perm-Z complete for %s||%s (n_ref=%d)", ti, tj, n0))
        } else {
          # Monte Carlo (consider subsetting genes for speed before calling this)
          pm <- .permZ_bipartite_mc(Mi, Mj, B = perm_B, cor_method = cor_method)
          # Map Z to r at reference n
          n_vec <- as.vector(n_mat[n_mat >= 4])
          n0 <- if (!length(n_vec)) 4L else if (is.character(n_ref[1])) {
            if (n_ref[1] == "max") max(n_vec) else stats::median(n_vec)
          } else as.integer(n_ref[1])
          n0 <- max(as.integer(round(n0)), 4L)
          S_perm_eff <- tanh(pm$Z / sqrt(n0 - 3))
          dimnames(S_perm_eff) <- dimnames(pm$Z)
          message(sprintf("[compare_correlation_methods] Monte Carlo perm-Z complete for %s||%s (B=%d, n_ref=%d)", ti, tj, perm_B, n0))
        }

        # CorShrink via ashr
        S_cs <- .corshrink_bipartite(S_raw, n_mat, ash_method = ash_method)
        message(sprintf("[compare_correlation_methods] CorShrink shrinkage complete for %s||%s", ti, tj))

        # Collect tidy edges (sampled)
        edf <- .sample_edges_df(
          S_list = list(raw_r = S_raw, perm_reff = S_perm_eff, corshrink_r = S_cs),
          n_mat = n_mat, tissueA = ti, tissueB = tj, edge_sample = edge_sample
        )
        edf$pair <- paste(ti, tj, sep = "||")
        all_edges[[length(all_edges) + 1L]] <- edf
        message(sprintf("[compare_correlation_methods] Sampled %d edges for %s||%s", nrow(edf), ti, tj))

        summ <- .summarize_n_dependence(edf, methods = c("raw_r","perm_reff","corshrink_r"))
        summ$tissueA <- ti; summ$tissueB <- tj; summ$pair <- paste(ti, tj, sep = "||")
        summ$n_ref   <- n0
        all_summ[[length(all_summ) + 1L]] <- summ
        message(sprintf("[compare_correlation_methods] Summarized n dependence for %s||%s", ti, tj))

        if (isTRUE(make_plots)) {
          .plot_n_vs_r(edf, out_prefix = paste0(out_prefix, "_", ti, "__", tj))
          message(sprintf("[compare_correlation_methods] Plots written for %s||%s", ti, tj))
        }
      }
    }
  }

  edges_df <- if (length(all_edges)) do.call(rbind, all_edges) else data.frame()
  summary_df <- if (length(all_summ)) do.call(rbind, all_summ) else data.frame()

  # Write outputs
  dir.create("output", showWarnings = FALSE, recursive = TRUE)
  if (nrow(edges_df)) utils::write.csv(edges_df, file = file.path("output", paste0(out_prefix, "_edges_sample.csv")), row.names = FALSE)
  if (nrow(summary_df)) utils::write.csv(summary_df, file = file.path("output", paste0(out_prefix, "_n_dependence_summary.csv")), row.names = FALSE)
  message(sprintf("[compare_correlation_methods] Outputs written under output/ with prefix '%s'", out_prefix))

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
  message(sprintf("[build_raw_adjacency_beta_fixed] Starting adjacency build for %d tissues", T))

  expr_list   <- vector("list", T)
  donor_list  <- vector("list", T)
  idx <- integer(T+1); rc_names <- character(0)
  for (i in seq_len(T)) {
    message(sprintf("[build_raw_adjacency_beta_fixed] Loading tissue %s from %s", tissue_names[i], tissue_expr_file_names[i]))
    X <- LoadExprData(tissue_names[i], tissue_expr_file_names[i],
                      sd_quantile = sd_quantile, max_genes_per_tissue = max_genes_per_tissue)
    expr_list[[i]]  <- X
    donor_list[[i]] <- .aggregate_by_donor(X)
    idx[i+1] <- idx[i] + ncol(X)
    rc_names <- c(rc_names, colnames(X))
    message(sprintf("[build_raw_adjacency_beta_fixed] Aggregated donors for %s (genes=%d, donors=%d)", tissue_names[i], ncol(X), nrow(donor_list[[i]])))
  }

  A <- matrix(0, nrow = idx[T+1], ncol = idx[T+1], dimnames = list(rc_names, rc_names))
  message("[build_raw_adjacency_beta_fixed] Initialized adjacency matrix")

  # TS blocks (within)
  for (i in seq_len(T)) {
    rows <- (idx[i]+1):idx[i+1]
    Sii <- abs(stats::cor(expr_list[[i]], use = "pairwise.complete.obs", method = cor_method))
    diag(Sii) <- 0
    A[rows, rows] <- Sii ^ TS_beta
    message(sprintf("[build_raw_adjacency_beta_fixed] Filled TS block for %s", tissue_names[i]))
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
        message(sprintf("[build_raw_adjacency_beta_fixed] Filled CT block for %s||%s (common_donors=%d)", tissue_names[i], tissue_names[j], length(common)))
      }
    }
  }
  message("[build_raw_adjacency_beta_fixed] Completed adjacency matrix")
  A
}
