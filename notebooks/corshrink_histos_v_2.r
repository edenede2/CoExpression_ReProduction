#!/usr/bin/env Rscript

# ---- Auto-install and load packages (first run may take a minute) ----
cran <- getOption("repos")["CRAN"]
if (is.null(cran) || identical(cran, "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
need_pkgs <- c(
  "optparse", "data.table", "matrixStats", "ggplot2", "stringr", "glue",
  "WGCNA", "ashr", "fastcluster", "dynamicTreeCut", "mixsqp"
)
inst <- need_pkgs[!vapply(need_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(inst)) {
  install.packages(inst, dependencies = TRUE)
}

suppressPackageStartupMessages({
  for (p in need_pkgs) {
    library(p, character.only = TRUE)
  }
  library(parallel)
})

options(stringsAsFactors = FALSE)
disableWGCNAThreads()

# ----------------------- CLI -----------------------
opt_list <- list(
  make_option(c("--input_dir"), type="character", default=NULL,
              help="Directory with CSV/CSV.GZ files (each = one tissue)."),
  make_option(c("--files"), type="character", default=NULL,
              help="Comma-separated list of CSV/CSV.GZ files (each = one tissue)."),
  make_option(c("--out_dir"), type="character", default="out_corshrink",
              help="Output directory [default: %default]"),
  make_option(c("--gene_col_name"), type="character", default="auto",
              help="Gene column name, or 'auto' to detect [default: %default]"),
  make_option(c("--gene_id_clean"), type="character", default="none",
              help="Standardize gene IDs across files: none | ensembl_strip_version | upper | lower [default: %default]"),
  make_option(c("--id_map_mode"), type="character", default="gtex",
              help="How to map sample columns to donor IDs for cross-tissue overlap: gtex | regex | none [default: %default]"),
  make_option(c("--donor_regex"), type="character", default="^([^-]+-[^-]+)",
              help="Regex with one capturing group that extracts donor ID from sample column names (used if --id_map_mode=regex) [default: %default]"),
  make_option(c("--collapse_fun"), type="character", default="mean",
              help="How to collapse multiple samples from the same donor within a tissue: mean | median | first [default: %default]"),
  make_option(c("--within_max_genes"), type="integer", default=2000,
              help="Max genes per tissue for within-tissue analysis [default: %default]"),
  make_option(c("--within_subset_mode"), type="character", default="topvar",
              help="Subset genes for within-tissue: topvar | random [default: %default]"),
  make_option(c("--within_max_pairs_for_plots"), type="integer", default=1e6,
              help="Max pairs sampled for within-tissue plots/stats [default: %default]"),
  make_option(c("--between_max_genes"), type="integer", default=0,
              help="Optional: limit common genes per between-tissue pair (0 = no limit) [default: %default]"),
  make_option(c("--between_subset_mode"), type="character", default="topvar",
              help="If limiting between genes: topvar | random [default: %default]"),
  make_option(c("--min_overlap_donors"), type="integer", default=4,
              help="Minimum overlapping donors required per tissue-pair [default: %default]"),
  make_option(c("--min_common_genes"), type="integer", default=50,
              help="Minimum common genes required per tissue-pair [default: %default]"),
  make_option(c("--threads"), type="integer", default=4,
              help="Threads for WGCNA::cor [default: %default]"),
  make_option(c("--seed"), type="integer", default=123,
              help="Random seed [default: %default]"),
  make_option(c("--skip_within"), type="logical", default=FALSE,
              help="Skip within-tissue analysis [default: %default]"),
  make_option(c("--skip_between"), type="logical", default=FALSE,
              help="Skip between-tissue (same-gene) analysis [default: %default]"),
  make_option(c("--cor_mode"), type="character", default="signed",
              help="Use 'signed' (keep sign) or 'unsigned' (use |r| for plots/stats) [default: %default]"),
  make_option(c("--cv_repeats"), type="integer", default=3,
              help="Number of CV repeats for between-tissue validation [default: %default]"),
  make_option(c("--cv_frac_train"), type="double", default=0.7,
              help="Fraction of donors for training in CV [default: %default]"),
  make_option(c("--reliab_k"), type="integer", default=50,
              help="Number of donors per split for reliability bootstrap [default: %default]"),
  make_option(c("--reliab_repeats"), type="integer", default=2,
              help="Number of bootstrap repeats for reliability [default: %default]"),
  make_option(c("--lfsr_threshold"), type="double", default=0.05,
              help="LFSR threshold for significance [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
set.seed(opt$seed)

if (is.null(opt$input_dir) && is.null(opt$files)) {
  stop("Provide either --input_dir or --files")
}

# ----------------------- IO helpers -----------------------
list_input_files <- function(input_dir, files_arg) {
  paths <- character(0)
  if (!is.null(input_dir)) {
    stopifnot(dir.exists(input_dir))
    a <- list.files(input_dir, pattern="\\.(csv|CSV)(\\.gz)?$", full.names = TRUE)
    paths <- c(paths, a)
  }
  if (!is.null(files_arg)) {
    b <- str_split(files_arg, ",", simplify = TRUE)
    b <- trimws(as.character(b[b != ""]))
    if (length(b)) paths <- c(paths, b)
  }
  paths <- unique(paths)
  if (!length(paths)) stop("No CSV files found.")
  paths
}

extract_tissue_name <- function(path) {
  base <- basename(path)
  base <- sub("\\.csv(\\.gz)?$", "", base, ignore.case = TRUE)
  base <- sub("^GTEX_(full|young|old)_", "", base, ignore.case = TRUE)
  base <- sub("^GTEX_", "", base, ignore.case = TRUE)
  base <- gsub("[^A-Za-z0-9]+", "", base)  # compact readable key
  if (nchar(base) == 0) base <- tools::file_path_sans_ext(basename(path))
  base
}

# Robust gene-column detection
detect_gene_col <- function(dt, gene_col_name) {
  if (gene_col_name != "auto") {
    if (!(gene_col_name %in% names(dt))) stop(glue("Gene column '{gene_col_name}' not found"))
    return(gene_col_name)
  }
  cand <- intersect(names(dt), c("gene","gene_id","geneid","symbol","ensembl","ensembl_id","Gene","GeneID"))
  if (length(cand) >= 1) return(cand[1])
  is_num <- vapply(dt, is.numeric, logical(1))
  nn <- which(!is_num)
  if (length(nn) >= 1) return(names(dt)[nn[1]])
  names(dt)[1]
}

# Gene ID standardization across files
clean_gene_ids <- function(x, mode) {
  x <- trimws(as.character(x))
  mode <- tolower(mode)
  if (mode == "ensembl_strip_version") {
    x <- sub("\\.\\d+$", "", x)
  } else if (mode == "upper") {
    x <- toupper(x)
  } else if (mode == "lower") {
    x <- tolower(x)
  }
  x
}

# Map sample columns to donors and collapse samples per donor (GTEx-friendly)
map_and_collapse_donors <- function(X, mode="gtex", donor_regex="^([^-]+-[^-]+)", fun="mean") {
  cn <- trimws(colnames(X))
  mode <- tolower(mode); fun <- tolower(fun)
  donor_id <- switch(mode,
    gtex = sub(paste0(donor_regex, ".*"), "\\1", cn),   # Fixed: \\1 instead of \\\\1
    regex = sub(paste0(donor_regex, ".*"), "\\1", cn),  # Fixed: \\1 instead of \\\\1
    none  = cn,
    cn
  )
  donor_id[donor_id == "" | is.na(donor_id)] <- cn[donor_id == "" | is.na(donor_id)]
  
  # Debug: show how many unique donors we found
  u <- unique(donor_id)
  if (mode != "none" && length(cn) > 10) {
    first_few <- paste(head(u, 3), collapse=", ")
    message(glue("  Donor extraction: {length(cn)} samples -> {length(u)} unique donors (first few: {first_few})"))
  }
  agg <- switch(fun,
    mean = function(m) rowMeans(m, na.rm=TRUE),
    median = function(m) matrixStats::rowMedians(m, na.rm=TRUE),
    first = function(m) m[,1,drop=FALSE],
    function(m) rowMeans(m, na.rm=TRUE)
  )
  # Use approach similar to ClusteringBuilding.R
  split_idx <- split(seq_len(ncol(X)), donor_id)
  res <- matrix(NA_real_, nrow=nrow(X), ncol=length(split_idx))
  
  for (i in seq_along(split_idx)) {
    idx <- split_idx[[i]]
    if (length(idx) == 1L) {
      res[, i] <- X[, idx]
    } else {
      res[, i] <- agg(X[, idx, drop=FALSE])
    }
  }
  
  colnames(res) <- names(split_idx)
  rownames(res) <- rownames(X)
  storage.mode(res) <- "double"
  res
}



load_tissue_csv <- function(path, gene_col_name="auto") {
  tissue <- extract_tissue_name(path)
  dt <- suppressWarnings(data.table::fread(path, nThread=1, showProgress=FALSE))
  if (ncol(dt) < 2) stop(glue("[{tissue}] CSV must have >=2 columns"))

  setnames(dt, old=names(dt), new=trimws(names(dt)))
  is_num <- vapply(dt, is.numeric, logical(1))
  non_num_cols <- names(dt)[!is_num]
  num_cols     <- names(dt)[ is_num]

  # האם זה samples×genes? (עמודה/ות לא-מספריות לדגימות, המון עמודות מספריות שנראות כמו ENSEMBL)
  ensg_pat <- "^ENSG[0-9]+(\\.[0-9]+)?$"
  num_is_ensg <- sum(grepl(ensg_pat, num_cols))
  
  # Additional patterns for gene IDs
  other_gene_patterns <- c(
    "^ENSG",           # Any ENSEMBL gene
    "^[A-Z0-9]+$",     # Gene symbols (all caps)
    "^[A-Za-z0-9_-]+\\.[0-9]+$"  # Versioned IDs (fixed character class order)
  )
  
  # Count how many numeric columns look like genes
  total_gene_like <- num_is_ensg
  for (pat in other_gene_patterns) {
    total_gene_like <- total_gene_like + sum(grepl(pat, num_cols))
  }
  
  # Check if we have GTEx-like sample IDs in the first non-numeric column
  gtex_sample_pattern <- "^GTEX-[A-Z0-9]+-[0-9]+-SM-[A-Z0-9]+$"
  has_gtex_samples <- FALSE
  if (length(non_num_cols) >= 1) {
    first_col_vals <- dt[[non_num_cols[1]]]
    has_gtex_samples <- sum(grepl(gtex_sample_pattern, first_col_vals, ignore.case = TRUE)) >= 10
  }
  
  # More flexible detection criteria
  min_genes_threshold <- max(100, floor(0.3 * length(num_cols)))
  samples_x_genes <- (length(non_num_cols) >= 1 && 
                     (total_gene_like >= min_genes_threshold || has_gtex_samples) &&
                     length(num_cols) >= 100)  # Must have many numeric columns
  
  # Debug message
  if (samples_x_genes) {
    message(glue("[{tissue}] Detected samples×genes format: {nrow(dt)} samples, {length(num_cols)} gene columns"))
  } else {
    message(glue("[{tissue}] Using genes×samples format: {nrow(dt)} genes, {length(num_cols)} sample columns"))
  }

  if (samples_x_genes) {
    sample_col <- non_num_cols[1]                    # לדוגמה: 'Unnamed: 0'
    samples <- dt[[sample_col]]
    gene_cols <- num_cols
    clean_genes <- clean_gene_ids(gene_cols, opt$gene_id_clean)
    X <- as.matrix(dt[, ..gene_cols])
    storage.mode(X) <- "double"
    rownames(X) <- make.unique(as.character(samples))   # שורות = דגימות
    colnames(X) <- make.unique(clean_genes)             # עמודות = גנים
    m <- t(X)                                           # הופך ל- genes × samples
    return(list(tissue=tissue, expr=m,
                donors_before=ncol(m), genes_before=nrow(m)))
  }

  # אחרת: הפורמט הישן genes×samples
  gcol <- detect_gene_col(dt, gene_col_name)
  genes <- clean_gene_ids(dt[[gcol]], opt$gene_id_clean)
  donor_cols <- setdiff(names(dt)[is_num], gcol)
  if (!length(donor_cols)) {
    maybe_num <- setdiff(names(dt), gcol)
    for (nm in maybe_num) {
      suppressWarnings({
        x <- as.numeric(dt[[nm]])
        if (sum(is.finite(x)) >= 0.5*length(x)) dt[[nm]] <- x
      })
    }
    is_num <- vapply(dt, is.numeric, logical(1))
    donor_cols <- setdiff(names(dt)[is_num], gcol)
  }
  if (!length(donor_cols)) stop(glue("[{tissue}] No numeric donor columns detected."))
  expr_dt <- dt[, ..donor_cols]
  m <- as.matrix(expr_dt)
  storage.mode(m) <- "double"
  rownames(m) <- make.unique(as.character(genes))
  colnames(m) <- trimws(colnames(expr_dt))
  list(tissue=tissue, expr=m,
       donors_before=ncol(m), genes_before=nrow(m))
}



# ----------------------- math helpers -----------------------
upper_tri_vec <- function(M) {
  M[upper.tri(M, diag=FALSE)]
}

ash_shrink <- function(z, se) {
  keep <- is.finite(z) & is.finite(se) & se > 0
  zhat <- rep(NA_real_, length(z))
  if (sum(keep) > 0) {
    fit <- ashr::ash(z[keep], se[keep], mixcompdist="normal", method="fdr")
    zhat[keep] <- ashr::get_pm(fit)
  }
  zhat
}

apply_cor_mode <- function(r, mode) {
  if (identical(tolower(mode), "unsigned")) return(abs(r))
  r
}

spearman_absr_vs_n <- function(r, n) {
  ok <- is.finite(r) & is.finite(n) & n >= 4
  if (sum(ok) < 10) return(list(rho=NA_real_, p=NA_real_, N=sum(ok)))
  
  # Check for variation in n (sample sizes)
  abs_r <- abs(r[ok])
  n_vals <- n[ok]
  
  # If all sample sizes are the same (common in GTEx), skip Spearman correlation
  if (length(unique(n_vals)) < 2) {
    # Return NA for Spearman but include sample info for debugging
    return(list(rho=NA_real_, p=NA_real_, N=sum(ok), constant_n=unique(n_vals)[1]))
  }
  
  # If not enough variation in correlations, also skip
  if (length(unique(abs_r)) < 2) {
    return(list(rho=NA_real_, p=NA_real_, N=sum(ok)))
  }
  
  ct <- suppressWarnings(cor.test(abs_r, n_vals, method="spearman", exact=FALSE))
  list(rho=unname(ct$estimate), p=ct$p.value, N=sum(ok))
}

# Efficient row-wise correlation over per-row overlap (same genes across two tissues)
rowwise_cor_overlap <- function(X, Y) {
  stopifnot(identical(dim(X), dim(Y)))
  M <- (!is.na(X)) & (!is.na(Y))
  n <- rowSums(M)
  X2 <- X; X2[!M] <- 0
  Y2 <- Y; Y2[!M] <- 0
  mX <- rowSums(X2) / pmax(n, 1)
  mY <- rowSums(Y2) / pmax(n, 1)
  Xc <- sweep(X2, 1, mX, "-")
  Yc <- sweep(Y2, 1, mY, "-")
  Xc[!M] <- 0; Yc[!M] <- 0
  num <- rowSums(Xc * Yc)
  den <- sqrt(rowSums(Xc^2) * rowSums(Yc^2))
  r <- num / den
  r[!is.finite(r) | n < 4 | den == 0] <- NA_real_
  list(r=r, n=n)
}

plot_hist_pre_post <- function(r_pre, r_post, title, outfile) {
  NMAX <- 2e6
  if (length(r_pre) > NMAX) {
    ix <- sample.int(length(r_pre), NMAX)
    r_pre <- r_pre[ix]
    r_post <- r_post[ix]
  }
  df <- rbind(
    data.frame(r=r_pre, state="pre"),
    data.frame(r=r_post, state="post")
  )
  p <- ggplot(df, aes(x=r, fill=state)) +
    geom_histogram(aes(y=after_stat(density)), bins=80, alpha=0.5, position="identity") +
    geom_vline(xintercept=0, linetype="dashed") +
    theme_minimal(base_size = 12) +
    labs(title=title, x="correlation", y="density")
  ggsave(outfile, p, width=7, height=4.2, dpi=160)
}

# Cross-validation metrics for between-tissue analysis
cv_metrics_between <- function(X, Y, B, frac_train) {
  stopifnot(identical(dim(X), dim(Y)))
  donors <- colnames(X)
  n_train <- max(4, floor(length(donors) * frac_train))
  
  mse_pre <- mse_post <- numeric(B)
  spear_pre <- spear_post <- numeric(B)
  signacc_pre <- signacc_post <- numeric(B)
  
  for (b in seq_len(B)) {
    train_idx <- sample.int(length(donors), n_train)
    test_idx <- setdiff(seq_along(donors), train_idx)
    if (length(test_idx) < 4) next
    
    # Train
    X_train <- X[, train_idx, drop=FALSE]
    Y_train <- Y[, train_idx, drop=FALSE]
    rw_train <- rowwise_cor_overlap(X_train, Y_train)
    r_train <- rw_train$r
    n_train_vec <- rw_train$n
    ok_train <- is.finite(r_train) & is.finite(n_train_vec) & n_train_vec >= 4
    if (sum(ok_train) < 10) next
    
    z_train <- atanh(r_train[ok_train])
    se_train <- sqrt(1 / (n_train_vec[ok_train] - 3))
    zhat_train <- ash_shrink(z_train, se_train)
    
    # Test
    X_test <- X[, test_idx, drop=FALSE]
    Y_test <- Y[, test_idx, drop=FALSE]
    rw_test <- rowwise_cor_overlap(X_test, Y_test)
    r_test <- rw_test$r[ok_train]
    n_test_vec <- rw_test$n[ok_train]
    ok_test <- is.finite(r_test) & is.finite(n_test_vec) & n_test_vec >= 4
    if (sum(ok_test) < 10) next
    
    z_test <- atanh(r_test[ok_test])
    z_train_pred <- z_train[ok_test]
    zhat_train_pred <- zhat_train[ok_test]
    
    # Metrics
    mse_pre[b] <- mean((z_train_pred - z_test)^2, na.rm=TRUE)
    mse_post[b] <- mean((zhat_train_pred - z_test)^2, na.rm=TRUE)
    spear_pre[b] <- suppressWarnings(cor(z_train_pred, z_test, method="spearman", use="complete.obs"))
    spear_post[b] <- suppressWarnings(cor(zhat_train_pred, z_test, method="spearman", use="complete.obs"))
    signacc_pre[b] <- mean(sign(z_train_pred) == sign(z_test), na.rm=TRUE)
    signacc_post[b] <- mean(sign(zhat_train_pred) == sign(z_test), na.rm=TRUE)
  }
  
  list(
    cv_mse_z_pre = mean(mse_pre[is.finite(mse_pre)], na.rm=TRUE),
    cv_mse_z_post = mean(mse_post[is.finite(mse_post)], na.rm=TRUE),
    cv_spear_pre = mean(spear_pre[is.finite(spear_pre)], na.rm=TRUE),
    cv_spear_post = mean(spear_post[is.finite(spear_post)], na.rm=TRUE),
    cv_signacc_pre = mean(signacc_pre[is.finite(signacc_pre)], na.rm=TRUE),
    cv_signacc_post = mean(signacc_post[is.finite(signacc_post)], na.rm=TRUE)
  )
}

# Bootstrap reliability for between-tissue analysis
bootstrap_reliability_between <- function(X, Y, K, B) {
  stopifnot(identical(dim(X), dim(Y)))
  donors <- colnames(X)
  if (length(donors) < 2*K) {
    return(list(reliab_pre=NA_real_, reliab_post=NA_real_))
  }
  
  rel_pre <- rel_post <- numeric(B)
  
  for (b in seq_len(B)) {
    idx <- sample.int(length(donors), 2*K)
    idx1 <- idx[1:K]
    idx2 <- idx[(K+1):(2*K)]
    
    # Split 1
    X1 <- X[, idx1, drop=FALSE]
    Y1 <- Y[, idx1, drop=FALSE]
    rw1 <- rowwise_cor_overlap(X1, Y1)
    r1_pre <- rw1$r
    n1 <- rw1$n
    ok1 <- is.finite(r1_pre) & is.finite(n1) & n1 >= 4
    if (sum(ok1) < 10) next
    
    z1 <- atanh(r1_pre[ok1])
    se1 <- sqrt(1 / (n1[ok1] - 3))
    zhat1 <- ash_shrink(z1, se1)
    r1_post <- tanh(zhat1)
    
    # Split 2
    X2 <- X[, idx2, drop=FALSE]
    Y2 <- Y[, idx2, drop=FALSE]
    rw2 <- rowwise_cor_overlap(X2, Y2)
    r2_pre <- rw2$r[ok1]
    n2 <- rw2$n[ok1]
    ok2 <- is.finite(r2_pre) & is.finite(n2) & n2 >= 4
    if (sum(ok2) < 10) next
    
    z2 <- atanh(r2_pre[ok2])
    se2 <- sqrt(1 / (n2[ok2] - 3))
    zhat2 <- ash_shrink(z2, se2)
    r2_post <- tanh(zhat2)
    
    # Reliability
    r1_pre_ok <- r1_pre[ok1][ok2]
    r1_post_ok <- r1_post[ok2]
    rel_pre[b] <- suppressWarnings(cor(r1_pre_ok, r2_pre[ok2], method="spearman", use="complete.obs"))
    rel_post[b] <- suppressWarnings(cor(r1_post_ok, r2_post, method="spearman", use="complete.obs"))
  }
  
  list(
    reliab_pre = mean(rel_pre[is.finite(rel_pre)], na.rm=TRUE),
    reliab_post = mean(rel_post[is.finite(rel_post)], na.rm=TRUE)
  )
}

# ----------------------- load inputs -----------------------
paths <- list_input_files(opt$input_dir, opt$files)
dir.create(opt$out_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$out_dir, "hists_within"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opt$out_dir, "hists_between"), showWarnings=FALSE, recursive=TRUE)
summary_path <- file.path(opt$out_dir, "summary_metrics.csv")
debug_inventory_path <- file.path(opt$out_dir, "debug_inventory.csv")
debug_pairs_path <- file.path(opt$out_dir, "debug_pairs.csv")
shrinkage_alpha_path <- file.path(opt$out_dir, "shrinkage_alpha_points.csv")
summary_dt <- data.table()
debug_inv <- data.table()
debug_pairs <- data.table()
alpha_points <- data.table()

message(glue("Found {length(paths)} file(s)."))

# Read all tissues
message("Loading tissues…")
tissues <- lapply(paths, load_tissue_csv, gene_col_name = opt$gene_col_name)
names(tissues) <- vapply(tissues, `[[`, "", "tissue")

# Clean gene IDs and collapse samples -> donors
for (i in seq_along(tissues)) {
  t <- tissues[[i]]
  # collapse to donors (GTEx-aware)
  X <- t$expr
  X <- map_and_collapse_donors(X, mode=opt$id_map_mode, donor_regex=opt$donor_regex, fun=opt$collapse_fun)
  # apply gene-id cleaning (already applied to row names when loading)
  rownames(X) <- clean_gene_ids(rownames(X), opt$gene_id_clean)
  tissues[[i]]$expr <- X
  tissues[[i]]$donors_after <- ncol(X)
  tissues[[i]]$genes_after <- nrow(X)
  debug_inv <- rbind(debug_inv, data.table(
    tissue = tissues[[i]]$tissue,
    donors_before = tissues[[i]]$donors_before,
    donors_after  = tissues[[i]]$donors_after,
    genes_before  = tissues[[i]]$genes_before,
    genes_after   = tissues[[i]]$genes_after
  ))
}

# Write inventory early
fwrite(debug_inv, debug_inventory_path)

# ----------------------- WITHIN-TISSUE -----------------------
if (!isTRUE(opt$skip_within)) {
  allowWGCNAThreads(nThreads = opt$threads)
  for (tiss in tissues) {
    tissue_name <- tiss$tissue
    X <- tiss$expr   # genes x donors
    if (ncol(X) < 4) {
      message(glue("[{tissue_name}] < 4 donors — skipping within-tissue"))
      next
    }
    g <- nrow(X)
    sel_idx <- seq_len(g)
    if (g > opt$within_max_genes) {
      if (opt$within_subset_mode == "topvar") {
        v <- rowVars(X, na.rm=TRUE)
        sel_idx <- order(v, decreasing = TRUE)[seq_len(opt$within_max_genes)]
      } else {
        sel_idx <- sample.int(g, opt$within_max_genes)
      }
    }
    Xs <- X[sel_idx, , drop=FALSE]
    p <- nrow(Xs)

    # Function to perform correlation analysis for a specific method
    analyze_correlation <- function(method_name) {
      message(glue("    Computing {method_name} correlations..."))
      
      # pairwise-complete correlation
      R <- suppressWarnings(WGCNA::cor(t(Xs), use="pairwise.complete.obs", method=tolower(method_name)))
      M <- t(!is.na(Xs))            # donors x genes (logical)
      Nmat <- crossprod(M)          # genes x genes (integer)
      
      # Debug: check the missing data pattern (only for first method to avoid repetition)
      if (method_name == "Pearson") {
        missing_per_gene <- rowSums(is.na(Xs))
        missing_per_sample <- colSums(is.na(Xs))
        message(glue("    Debug: Missing data - genes with missing: {sum(missing_per_gene > 0)}, samples with missing: {sum(missing_per_sample > 0)}"))
        message(glue("    Debug: Nmat range: [{min(Nmat)}, {max(Nmat)}], unique values: {length(unique(as.vector(Nmat)))}"))
        if (length(unique(as.vector(Nmat))) <= 5) {
          message(glue("    Debug: Nmat unique values: {paste(unique(as.vector(Nmat)), collapse=', ')}"))
        }
        message(glue("    Debug: R matrix dimensions: {nrow(R)}x{ncol(R)}, R range: [{round(min(R, na.rm=TRUE), 3)}, {round(max(R, na.rm=TRUE), 3)}]"))
      }

      r_pre <- upper_tri_vec(R)
      n_vec <- upper_tri_vec(Nmat)

      ok <- is.finite(r_pre) & is.finite(n_vec) & n_vec >= 4
      z <- atanh(r_pre[ok])
      se <- sqrt(1 / (n_vec[ok] - 3))
      
      # ASH shrinkage with diagnostics
      keep_ash <- is.finite(z) & is.finite(se) & se > 0
      zhat <- rep(NA_real_, length(z))
      lfsr <- rep(NA_real_, length(z))
      alpha <- rep(NA_real_, length(z))
      
      if (sum(keep_ash) > 0) {
        fit <- ashr::ash(z[keep_ash], se[keep_ash], mixcompdist="normal", method="fdr")
        zhat[keep_ash] <- ashr::get_pm(fit)
        lfsr[keep_ash] <- ashr::get_lfsr(fit)
        alpha[keep_ash] <- zhat[keep_ash] / z[keep_ash]
      }
      
      r_post <- r_pre
      r_post[ok] <- tanh(zhat)

      # apply mode for plots/stats
      r_pre_eff  <- apply_cor_mode(r_pre[ok],  opt$cor_mode)
      r_post_eff <- apply_cor_mode(r_post[ok], opt$cor_mode)
      n_eff <- n_vec[ok]

      # sample for plotting if huge
      if (length(r_pre_eff) > opt$within_max_pairs_for_plots) {
        ix <- sample.int(length(r_pre_eff), opt$within_max_pairs_for_plots)
        r_pre_eff  <- r_pre_eff[ix]
        r_post_eff <- r_post_eff[ix]
        n_eff      <- n_eff[ix]
      }

      plot_hist_pre_post(
        r_pre_eff, r_post_eff,
        title = glue("{tissue_name}: within-tissue {method_name} correlations ({opt$cor_mode})"),
        outfile = file.path(opt$out_dir, "hists_within", glue("within_{tissue_name}_{tolower(method_name)}.png"))
      )

      # Debug: check the correlation values
      if (length(r_pre_eff) > 0) {
        message(glue("    Debug {method_name}: r_pre_eff range: [{round(min(r_pre_eff, na.rm=TRUE), 3)}, {round(max(r_pre_eff, na.rm=TRUE), 3)}], n_eff range: [{min(n_eff, na.rm=TRUE)}, {max(n_eff, na.rm=TRUE)}]"))
        message(glue("    Debug {method_name}: Total pairs analyzed: {length(r_pre_eff)}"))
      }
      
      # LFSR diagnostics
      prop_lfsr <- as.numeric(mean(lfsr < opt$lfsr_threshold, na.rm=TRUE))
      median_lfsr <- as.numeric(median(lfsr, na.rm=TRUE))
      q10_lfsr <- as.numeric(quantile(lfsr, 0.1, na.rm=TRUE, names=FALSE))
      
      # Alpha vs n Spearman
      alpha_spear <- as.numeric(suppressWarnings(cor(alpha, n_vec[ok], method="spearman", use="complete.obs")))
      
      # Store alpha points for global plot
      if (method_name == "Pearson") {
        alpha_points <<- rbind(alpha_points, data.table(
          tissue = tissue_name,
          level = "within",
          n = n_vec[ok],
          alpha = alpha
        ))
      }
      
      return(list(
        r_pre = r_pre_eff,
        r_post = r_post_eff,
        n_eff = n_eff,
        method = method_name,
        prop_lfsr = prop_lfsr,
        median_lfsr = median_lfsr,
        q10_lfsr = q10_lfsr,
        alpha_spear = alpha_spear
      ))
    }
    
    # Analyze both correlation methods
    pearson_results <- analyze_correlation("Pearson")
    spearman_results <- analyze_correlation("Spearman")
    
    # Use Pearson results for summary statistics (consistent with original behavior)
    r_pre_eff <- pearson_results$r_pre
    r_post_eff <- pearson_results$r_post  
    n_eff <- pearson_results$n_eff
    
    sp_pre  <- spearman_absr_vs_n(r_pre_eff,  n_eff)
    sp_post <- spearman_absr_vs_n(r_post_eff, n_eff)

    # Calculate metrics for n vs |r| analysis (within-tissue)
    abs_r_pre <- abs(r_pre_eff)
    abs_r_post <- abs(r_post_eff)
    var_abs_r_pre <- var(abs_r_pre, na.rm = TRUE)
    var_abs_r_post <- var(abs_r_post, na.rm = TRUE)
    mean_abs_r_pre <- mean(abs_r_pre, na.rm = TRUE)
    mean_abs_r_post <- mean(abs_r_post, na.rm = TRUE)
    
    # For within-tissue, n_donors is the same for all gene pairs (constant per tissue)
    unique_n_donors <- unique(n_eff)[1]
    
    summary_dt <- rbindlist(
      list(summary_dt,
      data.table(
        level = "within",
        unit  = tissue_name,
        mode  = opt$cor_mode,
        n_points = as.integer(length(r_pre_eff)),
        spearman_absr_n_pre  = as.numeric(sp_pre$rho),
        p_pre = as.numeric(sp_pre$p),
        spearman_absr_n_post = as.numeric(sp_post$rho),
        p_post = as.numeric(sp_post$p),
        mean_val_pre  = as.numeric(mean(r_pre_eff,  na.rm=TRUE)),
        mean_val_post = as.numeric(mean(r_post_eff, na.rm=TRUE)),
        # New metrics for n vs |r| analysis
        n_donors = as.integer(unique_n_donors),
        mean_abs_r_pre = as.numeric(mean_abs_r_pre),
        mean_abs_r_post = as.numeric(mean_abs_r_post),
        var_abs_r_pre = as.numeric(var_abs_r_pre),
        var_abs_r_post = as.numeric(var_abs_r_post),
        # CV metrics (NA for within-tissue)
        cv_mse_z_pre = as.numeric(NA),
        cv_mse_z_post = as.numeric(NA),
        cv_spear_pre = as.numeric(NA),
        cv_spear_post = as.numeric(NA),
        cv_signacc_pre = as.numeric(NA),
        cv_signacc_post = as.numeric(NA),
        # Reliability (NA for within-tissue)
        reliab_pre = as.numeric(NA),
        reliab_post = as.numeric(NA),
        # ASH diagnostics from Pearson
        prop_lfsr = as.numeric(pearson_results$prop_lfsr),
        median_lfsr = as.numeric(pearson_results$median_lfsr),
        q10_lfsr = as.numeric(pearson_results$q10_lfsr),
        alpha_spear = as.numeric(pearson_results$alpha_spear)
      )),
      use.names=TRUE, fill=TRUE
    )
    # Create informative message based on whether Spearman correlation is available
    if (is.na(sp_pre$rho) && !is.null(sp_pre$constant_n)) {
      message(glue("[{tissue_name}] within: N={length(r_pre_eff)}, constant_n={sp_pre$constant_n} (Spearman rho not applicable - no sample size variation)"))
    } else {
      message(glue("[{tissue_name}] within: N={length(r_pre_eff)}, rho_pre={round(sp_pre$rho,3)} -> rho_post={round(sp_post$rho,3)}"))
    }
  }
}

# ----------------------- BETWEEN-TISSUE (same gene) -----------------------
if (!isTRUE(opt$skip_between)) {
  tn <- names(tissues)
  if (length(tn) >= 2) {
    pairs <- t(combn(tn, 2))
    for (k in seq_len(nrow(pairs))) {
      t1 <- tissues[[ pairs[k,1] ]]
      t2 <- tissues[[ pairs[k,2] ]]
      nm1 <- t1$tissue; nm2 <- t2$tissue

      donors <- intersect(colnames(t1$expr), colnames(t2$expr))
      genes  <- intersect(rownames(t1$expr), rownames(t2$expr))
      debug_pairs <- rbind(debug_pairs, data.table(
        pair = glue("{nm1}__{nm2}"),
        donors_overlap = length(donors),
        genes_overlap  = length(genes)
      ))

      if (length(donors) < opt$min_overlap_donors) {
        message(glue("[{nm1} vs {nm2}] < {opt$min_overlap_donors} overlapping donors — skipping"))
        next
      }
      if (length(genes) < opt$min_common_genes) {
        message(glue("[{nm1} vs {nm2}] too few common genes (<{opt$min_common_genes}) — skipping"))
        next
      }

      if (opt$between_max_genes > 0 && length(genes) > opt$between_max_genes) {
        if (opt$between_subset_mode == "topvar") {
          v <- rowVars(t1$expr[genes, donors, drop=FALSE], na.rm=TRUE)
          ord <- order(v, decreasing = TRUE)
          genes <- genes[ord[seq_len(opt$between_max_genes)]]
        } else {
          genes <- sample(genes, opt$between_max_genes)
        }
      }

      X <- t1$expr[genes, donors, drop=FALSE]
      Y <- t2$expr[genes, donors, drop=FALSE]

      rw <- rowwise_cor_overlap(X, Y)
      r_pre <- rw$r
      n_vec <- rw$n
      ok <- is.finite(r_pre) & is.finite(n_vec) & n_vec >= 4
      if (sum(ok) < 50) {
        message(glue("[{nm1} vs {nm2}] too few valid pairs — skipping"))
        next
      }

      z <- atanh(r_pre[ok])
      se <- sqrt(1 / (n_vec[ok] - 3))
      
      # ASH shrinkage with diagnostics
      keep_ash <- is.finite(z) & is.finite(se) & se > 0
      zhat <- rep(NA_real_, length(z))
      lfsr <- rep(NA_real_, length(z))
      alpha <- rep(NA_real_, length(z))
      
      if (sum(keep_ash) > 0) {
        fit <- ashr::ash(z[keep_ash], se[keep_ash], mixcompdist="normal", method="fdr")
        zhat[keep_ash] <- ashr::get_pm(fit)
        lfsr[keep_ash] <- ashr::get_lfsr(fit)
        alpha[keep_ash] <- zhat[keep_ash] / z[keep_ash]
      }
      
      r_post <- r_pre
      r_post[ok] <- tanh(zhat)
      
      # LFSR diagnostics
      prop_lfsr <- as.numeric(mean(lfsr < opt$lfsr_threshold, na.rm=TRUE))
      median_lfsr <- as.numeric(median(lfsr, na.rm=TRUE))
      q10_lfsr <- as.numeric(quantile(lfsr, 0.1, na.rm=TRUE, names=FALSE))
      alpha_spear <- as.numeric(suppressWarnings(cor(alpha, n_vec[ok], method="spearman", use="complete.obs")))
      
      # Store alpha points
      alpha_points <- rbind(alpha_points, data.table(
        tissue = glue("{nm1}__{nm2}"),
        level = "between",
        n = n_vec[ok],
        alpha = alpha
      ))
      
      # Cross-validation and reliability
      cv <- cv_metrics_between(X, Y, B = opt$cv_repeats, frac_train = opt$cv_frac_train)
      rel <- bootstrap_reliability_between(X, Y, K = opt$reliab_k, B = opt$reliab_repeats)

      r_pre_eff  <- apply_cor_mode(r_pre[ok],  opt$cor_mode)
      r_post_eff <- apply_cor_mode(r_post[ok], opt$cor_mode)
      n_eff <- n_vec[ok]

      plot_hist_pre_post(
        r_pre_eff, r_post_eff,
        title = glue("{nm1} vs {nm2}: cross-tissue same-gene ({opt$cor_mode})"),
        outfile = file.path(opt$out_dir, "hists_between", glue("between_{nm1}__{nm2}.png"))
      )

      sp_pre  <- spearman_absr_vs_n(r_pre_eff,  n_eff)
      sp_post <- spearman_absr_vs_n(r_post_eff, n_eff)

      # Calculate metrics for n vs |r| analysis
      abs_r_pre <- abs(r_pre_eff)
      abs_r_post <- abs(r_post_eff)
      var_abs_r_pre <- var(abs_r_pre, na.rm = TRUE)
      var_abs_r_post <- var(abs_r_post, na.rm = TRUE)
      mean_abs_r_pre <- mean(abs_r_pre, na.rm = TRUE)
      mean_abs_r_post <- mean(abs_r_post, na.rm = TRUE)
      
      # Correlation between n_donors and abs(r) - use unique n value since all genes have same donor count
      unique_n_donors <- unique(n_eff)[1]  # All genes have same number of donors in between-tissue analysis
      
      summary_dt <- rbindlist(
        list(summary_dt,
        data.table(
          level = "between",
          unit  = glue("{nm1}__{nm2}"),
          mode  = opt$cor_mode,
          n_points = as.integer(length(r_pre_eff)),
          spearman_absr_n_pre  = as.numeric(sp_pre$rho),
          p_pre = as.numeric(sp_pre$p),
          spearman_absr_n_post = as.numeric(sp_post$rho),
          p_post = as.numeric(sp_post$p),
          mean_val_pre  = as.numeric(mean(r_pre_eff,  na.rm=TRUE)),
          mean_val_post = as.numeric(mean(r_post_eff, na.rm=TRUE)),
          # New metrics for n vs |r| analysis
          n_donors = as.integer(unique_n_donors),
          mean_abs_r_pre = as.numeric(mean_abs_r_pre),
          mean_abs_r_post = as.numeric(mean_abs_r_post),
          var_abs_r_pre = as.numeric(var_abs_r_pre),
          var_abs_r_post = as.numeric(var_abs_r_post),
          # CV metrics
          cv_mse_z_pre = as.numeric(cv$cv_mse_z_pre),
          cv_mse_z_post = as.numeric(cv$cv_mse_z_post),
          cv_spear_pre = as.numeric(cv$cv_spear_pre),
          cv_spear_post = as.numeric(cv$cv_spear_post),
          cv_signacc_pre = as.numeric(cv$cv_signacc_pre),
          cv_signacc_post = as.numeric(cv$cv_signacc_post),
          # Reliability
          reliab_pre = as.numeric(rel$reliab_pre),
          reliab_post = as.numeric(rel$reliab_post),
          # ASH diagnostics
          prop_lfsr = as.numeric(prop_lfsr),
          median_lfsr = as.numeric(median_lfsr),
          q10_lfsr = as.numeric(q10_lfsr),
          alpha_spear = as.numeric(alpha_spear)
        )),
        use.names=TRUE, fill=TRUE
      )
      message(glue("[{nm1} vs {nm2}] between: donors={length(donors)}, genes={length(genes)}, N={length(r_pre_eff)}, rho_pre={round(sp_pre$rho,3)} -> rho_post={round(sp_post$rho,3)}"))
    }
  }
}

# ----------------------- Fairness diagnostic -----------------------
if (nrow(summary_dt[level == "between"]) > 0) {
  between_dt <- summary_dt[level == "between"]
  q <- quantile(between_dt$n_donors, c(0.25, 0.75), na.rm=TRUE)
  q1_rows <- between_dt[n_donors <= q[1]]
  q3_rows <- between_dt[n_donors >= q[2]]
  
  fairness_gap <- data.table(
    fairness_gap_pre = median(q3_rows$mean_abs_r_pre, na.rm=TRUE) - median(q1_rows$mean_abs_r_pre, na.rm=TRUE),
    fairness_gap_post = median(q3_rows$mean_abs_r_post, na.rm=TRUE) - median(q1_rows$mean_abs_r_post, na.rm=TRUE)
  )
  
  fairness_gap_path <- file.path(opt$out_dir, "fairness_gap.csv")
  fwrite(fairness_gap, fairness_gap_path)
  message(glue("Fairness gap: {fairness_gap_path}"))
}

# ----------------------- Save summary and alpha points -----------------------
dir.create(dirname(summary_path), showWarnings = FALSE, recursive = TRUE)
fwrite(summary_dt, summary_path)
fwrite(debug_inv, debug_inventory_path)
fwrite(debug_pairs, debug_pairs_path)
fwrite(alpha_points, shrinkage_alpha_path)

# Plot shrinkage vs n
if (nrow(alpha_points) > 0) {
  alpha_plot_path <- file.path(opt$out_dir, "shrinkage_vs_n.png")
  alpha_clean <- alpha_points[is.finite(alpha) & is.finite(n)]
  if (nrow(alpha_clean) > 0) {
    p <- ggplot(alpha_clean, aes(x=n, y=alpha)) +
      geom_point(alpha=0.05, size=0.5) +
      geom_smooth(method="loess", color="red") +
      facet_wrap(~level, scales="free") +
      theme_minimal(base_size=12) +
      labs(title="Shrinkage factor (alpha = zhat/z) vs sample size",
           x="Sample size (n)", y="Alpha")
    ggsave(alpha_plot_path, p, width=10, height=5, dpi=160)
    message(glue("Shrinkage plot: {alpha_plot_path}"))
  }
}

message(glue("Done. Summary: {summary_path}"))
message(glue("Debug: {debug_inventory_path}, {debug_pairs_path}"))
message(glue("Alpha points: {shrinkage_alpha_path}"))
