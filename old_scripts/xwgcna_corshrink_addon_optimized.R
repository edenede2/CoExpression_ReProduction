# xwgcna_corshrink_addon_optimized.R
# Optimized version with performance improvements
suppressPackageStartupMessages({
  requireNamespace("WGCNA", quietly = TRUE)
  requireNamespace("parallel", quietly = TRUE)
})

# ==== CorShrink availability (cached check) ====
.corshrink_available <- NULL
.require_CorShrink <- function() {
  if (is.null(.corshrink_available)) {
    .corshrink_available <<- requireNamespace("CorShrink", quietly = TRUE)
  }
  if (!.corshrink_available) {
    stop("CorShrink is missing. Please install once:\n  devtools::install_github('kkdey/CorShrink')")
  }
}

# ==== Optimized CorShrink on union of donors ====
.corshrink_cross_union <- function(
  Mi, Mj, nameA="A", nameB="B",
  block_B = NULL,
  corshrink_max_p = 6000,
  min_nonmiss = 4L,
  min_sd = 1e-8,
  use_parallel = FALSE,
  n_cores = NULL
){
  .require_CorShrink()
  
  # Early exit for empty matrices
  if (ncol(Mi) == 0 || ncol(Mj) == 0) {
    return(matrix(0, nrow = ncol(Mi), ncol = ncol(Mj),
                  dimnames = list(colnames(Mi), colnames(Mj))))
  }
  
  genes_i <- colnames(Mi)
  genes_j <- colnames(Mj)
  
  # More efficient union computation
  donors_i <- rownames(Mi)
  donors_j <- rownames(Mj)
  donors_union <- union(donors_i, donors_j)
  
  # Vectorized column screening function
  col_screen_fast <- function(X) {
    if (ncol(X) == 0) return(X)
    
    # Vectorized operations
    finite_mask <- is.finite(X)
    nn <- colSums(finite_mask)
    
    # Fast SD computation using colSums and colMeans
    col_means <- colSums(X, na.rm = TRUE) / nn
    col_vars <- (colSums(X^2, na.rm = TRUE) - nn * col_means^2) / pmax(nn - 1, 1)
    ss <- sqrt(col_vars)
    
    keep <- (nn >= min_nonmiss) & is.finite(ss) & (ss > min_sd)
    X[, keep, drop = FALSE]
  }

  # Pre-allocate matrices more efficiently
  n_union <- length(donors_union)
  XA <- matrix(NA_real_, nrow = n_union, ncol = length(genes_i))
  dimnames(XA) <- list(donors_union, paste0(nameA, ":", genes_i))
  
  # Efficient matrix assignment using match
  if (length(donors_i) > 0) {
    idx_i <- match(donors_i, donors_union)
    XA[idx_i, ] <- as.matrix(Mi)
  }
  
  XA <- col_screen_fast(XA)
  genes_i_filtered <- sub(paste0("^", nameA, ":"), "", colnames(XA))

  # Pre-allocate result matrix
  AB <- matrix(NA_real_, nrow = length(genes_i_filtered), ncol = length(genes_j),
               dimnames = list(genes_i_filtered, genes_j))

  # Optimized CorShrink wrapper
  run_corshrink_fast <- function(X) {
    if (ncol(X) < 2) return(NULL)
    
    cs <- try(CorShrink::CorShrinkData(X), silent = TRUE)
    if (inherits(cs, "try-error")) return(NULL)
    
    # More efficient result extraction
    if (is.list(cs)) {
      S <- cs$cor %||% cs$ash_cor %||% as.matrix(cs)
    } else {
      S <- as.matrix(cs)
    }
    
    if (!is.matrix(S) || anyNA(S)) return(NULL)
    S
  }

  # Decide on full vs blocked approach
  total_genes <- length(genes_i_filtered) + length(genes_j)
  do_full <- total_genes <= corshrink_max_p && is.null(block_B)
  
  if (do_full) {
    # Full matrix approach - optimized
    XB <- matrix(NA_real_, nrow = n_union, ncol = length(genes_j))
    dimnames(XB) <- list(donors_union, paste0(nameB, ":", genes_j))
    
    if (length(donors_j) > 0) {
      idx_j <- match(donors_j, donors_union)
      XB[idx_j, ] <- as.matrix(Mj)
    }
    
    XB <- col_screen_fast(XB)
    
    if (ncol(XB) == 0) return(AB)
    
    Xfull <- cbind(XA, XB)
    message(sprintf("CorShrink (full): p=%d (|%s|=%d, |%s|=%d)", 
                   ncol(Xfull), nameA, ncol(XA), nameB, ncol(XB)))
    
    S <- run_corshrink_fast(Xfull)
    if (is.null(S)) return(NULL)

    # Faster index computation
    ncol_A <- ncol(XA)
    AB2 <- S[1:ncol_A, (ncol_A + 1):ncol(S), drop = FALSE]
    
    # Clean dimnames and clamp
    rownames(AB2) <- sub(paste0("^", nameA, ":"), "", rownames(AB2))
    colnames(AB2) <- sub(paste0("^", nameB, ":"), "", colnames(AB2))
    AB2 <- pmax(pmin(AB2, 1), -1)

    AB[rownames(AB2), colnames(AB2)] <- AB2
    return(AB)
  }

  # Blocked approach with optional parallelization
  if (is.null(block_B)) block_B <- 800L
  nb <- ceiling(length(genes_j) / block_B)
  
  if (nb == 0) return(AB)
  
  message(sprintf("CorShrink in blocks: |A|=%d, |B|=%d, block_B=%d", 
                 length(genes_i_filtered), length(genes_j), block_B))

  # Process blocks (optionally in parallel)
  process_block <- function(b) {
    j1 <- (b - 1L) * block_B + 1L
    j2 <- min(b * block_B, length(genes_j))
    gj <- genes_j[j1:j2]

    XB_block <- matrix(NA_real_, nrow = n_union, ncol = length(gj))
    dimnames(XB_block) <- list(donors_union, paste0(nameB, ":", gj))
    
    if (length(donors_j) > 0) {
      idx_j <- match(donors_j, donors_union)
      XB_block[idx_j, ] <- as.matrix(Mj[, gj, drop = FALSE])
    }
    
    XB_block <- col_screen_fast(XB_block)
    if (ncol(XB_block) == 0) return(NULL)

    X <- cbind(XA, XB_block)
    message(sprintf("  block %d/%d: p=%d", b, nb, ncol(X)))
    
    S <- run_corshrink_fast(X)
    if (is.null(S)) return(NULL)

    ncol_A <- ncol(XA)
    AB_block <- S[1:ncol_A, (ncol_A + 1):ncol(S), drop = FALSE]

    rownames(AB_block) <- sub(paste0("^", nameA, ":"), "", rownames(AB_block))
    colnames(AB_block) <- sub(paste0("^", nameB, ":"), "", colnames(AB_block))
    AB_block <- pmax(pmin(AB_block, 1), -1)
    
    list(block = b, result = AB_block, genes = gj)
  }

  # Execute blocks
  if (use_parallel && nb > 1) {
    if (is.null(n_cores)) n_cores <- min(parallel::detectCores() - 1, nb)
    block_results <- parallel::mclapply(seq_len(nb), process_block, mc.cores = n_cores)
  } else {
    block_results <- lapply(seq_len(nb), process_block)
  }

  # Combine results
  for (result in block_results) {
    if (!is.null(result) && !is.null(result$result)) {
      AB[rownames(result$result), colnames(result$result)] <- result$result
    }
  }
  
  AB
}

# ==== Optimized naive cross-tissue correlation ====
.cor_cross_common <- function(Mi, Mj, cor_method="pearson", min_common=3L) {
  # Early exit
  if (ncol(Mi) == 0 || ncol(Mj) == 0) {
    return(matrix(0, nrow = ncol(Mi), ncol = ncol(Mj),
                  dimnames = list(colnames(Mi), colnames(Mj))))
  }
  
  # Faster intersection using match
  donors_i <- rownames(Mi)
  donors_j <- rownames(Mj)
  common <- intersect(donors_i, donors_j)
  
  if (length(common) < min_common) {
    return(matrix(0, nrow = ncol(Mi), ncol = ncol(Mj),
                  dimnames = list(colnames(Mi), colnames(Mj))))
  }
  
  # More efficient subsetting
  Mi_common <- Mi[common, , drop = FALSE]
  Mj_common <- Mj[common, , drop = FALSE]
  
  S <- stats::cor(Mi_common, Mj_common, use = "pairwise.complete.obs", method = cor_method)
  S[!is.finite(S)] <- 0
  S <- pmax(pmin(S, 1), -1)
  dimnames(S) <- list(colnames(Mi), colnames(Mj))
  S
}

# ==== Vectorized beta curve computation ====
.CT_beta_curve_from_S <- function(S, powerVector = c(1:10, seq(12,20,2)),
                                  TOMType = "unsigned", unsigned = TRUE,
                                  nBreaks = 12, removeFirst = TRUE) {
  if (!is.matrix(S)) S <- as.matrix(S)
  S[!is.finite(S)] <- 0
  S <- pmax(pmin(S, 1), -1)
  if (unsigned) S <- abs(S)

  # Pre-allocate result
  n_powers <- length(powerVector)
  fit_df <- data.frame(
    Power = powerVector,
    SFT.R.sq = numeric(n_powers),
    slope = numeric(n_powers),
    mean.k = numeric(n_powers),
    median.k = numeric(n_powers),
    max.k = numeric(n_powers),
    stringsAsFactors = FALSE
  )
  
  # Optimized adjacency function
  is_signed <- startsWith(tolower(TOMType), "signed")
  
  # Vectorized computation where possible
  for (ix in seq_len(n_powers)) {
    b <- powerVector[ix]
    
    if (is_signed) {
      A <- ((S + 1) / 2) ^ b
    } else {
      A <- S ^ b
    }
    
    # Efficient connectivity calculation
    k_row <- rowSums(A, na.rm = TRUE)
    k_col <- colSums(A, na.rm = TRUE)
    k <- c(k_row, k_col)
    k <- k[is.finite(k) & k > 0]
    
    if (length(k) == 0) next
    
    sf <- WGCNA::scaleFreeFitIndex(k, nBreaks = nBreaks, removeFirst = removeFirst)
    fit_df$SFT.R.sq[ix] <- sf$Rsquared.SFT
    fit_df$slope[ix] <- sf$slope.SFT
    fit_df$mean.k[ix] <- mean(k)
    fit_df$median.k[ix] <- stats::median(k)
    fit_df$max.k[ix] <- max(k)
  }
  
  fit_df
}

# ==== Optimized before/after plot ====
plot_CT_beta_before_after <- function(
  Mi, Mj, pair_id,
  cor_method = "pearson",
  powerVector = c(1:10, seq(12,20,2)),
  TOMType = "unsigned", unsigned = TRUE,
  min_common = 3L, out_dir = "plots",
  nameA = "A", nameB = "B",
  nBreaks = 12, removeFirst = TRUE,
  use_parallel = FALSE
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Compute correlations (reuse optimized functions)
  S_before <- .cor_cross_common(Mi, Mj, cor_method = cor_method, min_common = min_common)
  S_after <- .corshrink_cross_union(Mi, Mj, nameA = nameA, nameB = nameB, 
                                   use_parallel = use_parallel)
  
  if (is.null(S_after)) S_after <- S_before
  
  # Fill holes
  mask <- !is.finite(S_after)
  if (any(mask)) S_after[mask] <- S_before[mask]
  S_after <- pmax(pmin(S_after, 1), -1)

  # Compute beta curves
  df_before <- .CT_beta_curve_from_S(S_before, powerVector, TOMType, unsigned,
                                    nBreaks = nBreaks, removeFirst = removeFirst)
  df_after <- .CT_beta_curve_from_S(S_after, powerVector, TOMType, unsigned,
                                   nBreaks = nBreaks, removeFirst = removeFirst)

  # Generate plot
  fn <- file.path(out_dir, paste0("CT_beta_before_after__", 
                                 gsub("[^A-Za-z0-9]+", "_", pair_id), ".pdf"))
  
  pdf(fn, width = 6, height = 4)
  on.exit(dev.off(), add = TRUE)
  
  yr <- range(c(df_before$SFT.R.sq, df_after$SFT.R.sq), na.rm = TRUE)
  if (!all(is.finite(yr))) yr <- c(0, 1)
  
  plot(df_before$Power, df_before$SFT.R.sq, type = "o", pch = 16, lwd = 1.2, 
       col = "grey30", xlab = expression(beta), ylab = expression(R^2), 
       main = paste0("CT beta-curves — ", pair_id), ylim = yr)
  
  lines(df_after$Power, df_after$SFT.R.sq, type = "o", pch = 17, lwd = 1.4, 
        col = "steelblue")
  
  abline(h = 0.8, lty = 3, col = "grey70")
  legend("bottomright", bty = "n", pch = c(16, 17), lwd = c(1.2, 1.4),
         col = c("grey30", "steelblue"),
         legend = c("Before (common donors)", "After (CorShrink, all donors)"))
  
  invisible(fn)
}

# ==== Optimized main adjacency function ====
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
    # optimization flags:
    ct_use_all_donors = FALSE,
    ct_use_corshrink = TRUE,
    plot_CT_beta_before_after = TRUE,
    powerVector_for_plots = c(1:10, seq(12,20,2)),
    TOMType_for_plots = "unsigned",
    beta_curves_dir = "plots",
    save_intermediates = FALSE,
    out_prefix = "xwgcna",
    corshrink_max_p = 6000,
    block_B = 800L,
    use_parallel = FALSE,
    n_cores = NULL,
    # memory optimization
    preallocate_matrices = TRUE
) {
  ct_too_few_action <- match.arg(ct_too_few_action)
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)

  # Optimized data loading with progress
  message("Loading expression data...")
  expr_list <- vector("list", T)
  donor_mat <- vector("list", T)
  
  # Pre-compute sizes for efficient matrix allocation
  gene_counts <- integer(T)
  for (i in seq_len(T)) {
    X <- LoadExprData(
      tissue_names[i], tissue_expr_file_names[i],
      sd_quantile = sd_quantile,
      max_genes_per_tissue = max_genes_per_tissue
    )
    expr_list[[i]] <- X
    donor_mat[[i]] <- .aggregate_by_donor(X)
    gene_counts[i] <- ncol(X)
  }
  
  # Efficient index computation
  total_genes <- sum(gene_counts)
  idx <- c(0, cumsum(gene_counts))
  rc_names <- unlist(lapply(expr_list, colnames), use.names = FALSE)

  # Pre-allocate matrices efficiently
  if (preallocate_matrices) {
    A_final <- matrix(0, nrow = total_genes, ncol = total_genes)
    dimnames(A_final) <- list(rc_names, rc_names)
    
    A_pre <- if (save_intermediates) {
      matrix(0, nrow = total_genes, ncol = total_genes, dimnames = dimnames(A_final))
    } else NULL
    
    A_post <- if (save_intermediates) {
      matrix(0, nrow = total_genes, ncol = total_genes, dimnames = dimnames(A_final))
    } else NULL
  } else {
    # Fallback for memory-constrained systems
    A_final <- matrix(0, nrow = total_genes, ncol = total_genes, 
                     dimnames = list(rc_names, rc_names))
    A_pre <- A_post <- NULL
  }

  # ---- Optimized TS computation ----
  message("Computing tissue-specific adjacencies...")
  for (i in seq_len(T)) {
    rows <- (idx[i] + 1):idx[i + 1]
    
    pow_i <- if (!is.null(TS_power_map) && !is.na(TS_power_map[tissue_names[i]])) {
      as.integer(TS_power_map[tissue_names[i]])
    } else default_TS
    
    # More efficient correlation computation
    Sii <- stats::cor(expr_list[[i]], use = "pairwise.complete.obs", method = cor_method)
    
    if (!signed) Sii <- abs(Sii)
    Sii[!is.finite(Sii)] <- 0
    diag(Sii) <- 0
    
    # Vectorized power operation
    Aii <- Sii ^ pow_i
    Aii <- pmax(pmin(Aii, 1), 0)
    
    A_final[rows, rows] <- Aii
    
    if (save_intermediates) {
      if (!is.null(A_pre)) A_pre[rows, rows] <- Aii
      if (!is.null(A_post)) A_post[rows, rows] <- Aii
    }
  }

  # ---- Optimized CT computation ----
  if (T >= 2) {
    message("Computing cross-tissue adjacencies...")
    
    # Pre-compute all tissue pairs
    pairs <- combn(T, 2, simplify = FALSE)
    n_pairs <- length(pairs)
    
    for (p in seq_len(n_pairs)) {
      i <- pairs[[p]][1]
      j <- pairs[[p]][2]
      
      rows <- (idx[i] + 1):idx[i + 1]
      cols <- (idx[j] + 1):idx[j + 1]
      
      pair_ij <- paste(tissue_names[i], tissue_names[j], sep = "||")
      pair_ji <- paste(tissue_names[j], tissue_names[i], sep = "||")
      
      pow_ij <- if (!is.null(CT_power_map)) {
        if (!is.na(CT_power_map[pair_ij])) as.integer(CT_power_map[pair_ij])
        else if (!is.na(CT_power_map[pair_ji])) as.integer(CT_power_map[pair_ji])
        else default_CT
      } else default_CT

      # Always compute baseline for comparison/fallback
      S_before <- .cor_cross_common(donor_mat[[i]], donor_mat[[j]],
                                   cor_method = cor_method, min_common = ct_min_common)

      if (isTRUE(ct_use_all_donors) && isTRUE(ct_use_corshrink)) {
        S_after <- .corshrink_cross_union(
          donor_mat[[i]], donor_mat[[j]],
          nameA = tissue_names[i], nameB = tissue_names[j],
          block_B = block_B, corshrink_max_p = corshrink_max_p,
          min_nonmiss = 4L, min_sd = 1e-8,
          use_parallel = use_parallel, n_cores = n_cores
        )

        if (is.null(S_after)) {
          message(sprintf("[CT %s] CorShrink failed → fallback to common donors", pair_ij))
          S_after <- S_before
        } else {
          # Efficient clamping and hole-filling
          S_after <- pmax(pmin(S_after, 1), -1)
          na_mask <- !is.finite(S_after)
          if (any(na_mask)) S_after[na_mask] <- S_before[na_mask]
        }

        # Optional plotting (skip if too many pairs)
        if (isTRUE(plot_CT_beta_before_after) && n_pairs <= 20) {
          try(
            plot_CT_beta_before_after(
              donor_mat[[i]], donor_mat[[j]], pair_id = pair_ij,
              cor_method = cor_method, powerVector = powerVector_for_plots,
              TOMType = TOMType_for_plots, unsigned = !signed,
              min_common = ct_min_common, out_dir = beta_curves_dir,
              nameA = tissue_names[i], nameB = tissue_names[j],
              use_parallel = use_parallel
            ),
            silent = TRUE
          )
        }

        # Efficient adjacency computation
        if (signed) {
          C_pre <- ((S_before + 1) / 2) ^ pow_ij
          C_post <- ((S_after + 1) / 2) ^ pow_ij
        } else {
          C_pre <- abs(S_before) ^ pow_ij
          C_post <- abs(S_after) ^ pow_ij
        }
        
        C_pre <- pmax(pmin(C_pre, 1), 0)
        C_post <- pmax(pmin(C_post, 1), 0)

        A_final[rows, cols] <- C_post
        A_final[cols, rows] <- t(C_post)
        
        if (save_intermediates) {
          if (!is.null(A_pre)) {
            A_pre[rows, cols] <- C_pre
            A_pre[cols, rows] <- t(C_pre)
          }
          if (!is.null(A_post)) {
            A_post[rows, cols] <- C_post
            A_post[cols, rows] <- t(C_post)
          }
        }
      } else {
        # Standard approach only
        if (signed) {
          C <- ((S_before + 1) / 2) ^ pow_ij
        } else {
          C <- abs(S_before) ^ pow_ij
        }
        
        C <- pmax(pmin(C, 1), 0)
        A_final[rows, cols] <- C
        A_final[cols, rows] <- t(C)
        
        if (save_intermediates) {
          if (!is.null(A_pre)) {
            A_pre[rows, cols] <- C
            A_pre[cols, rows] <- t(C)
          }
          if (!is.null(A_post)) {
            A_post[rows, cols] <- C
            A_post[cols, rows] <- t(C)
          }
        }
      }
    }
  }

  # ---- Final adjacency cleanup (vectorized) ----
  message("Finalizing adjacency matrix...")
  A_final[!is.finite(A_final)] <- 0
  A_final <- pmax(pmin(A_final, 1), 0)
  diag(A_final) <- 0
  A_final <- (A_final + t(A_final)) / 2

  # Save intermediates if requested
  if (isTRUE(save_intermediates)) {
    dir.create("output", showWarnings = FALSE, recursive = TRUE)
    if (!is.null(A_pre)) {
      saveRDS(A_pre, file = file.path("output", paste0(out_prefix, "_adjacency_preNorm.rds")))
    }
    if (!is.null(A_post)) {
      saveRDS(A_post, file = file.path("output", paste0(out_prefix, "_adjacency_postNorm_CorShrink.rds")))
    }
  }

  message("Adjacency matrix computation complete.")
  A_final
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x