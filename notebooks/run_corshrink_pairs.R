
# run_corshrink_pairs.R
# Matthew Stephens-style Empirical Bayes shrinkage of cross-tissue correlations using CorShrink
# This script:
#  1) Loads three tissue expression CSVs (first column = sample IDs; other columns = genes).
#  2) For each tissue pair, builds an NA-aware Donors×(A:genes + B:genes) matrix on COMMON genes.
#  3) Runs CorShrinkData(), extracts the A×B cross-tissue correlation block, writes:
#       - *_CorShrunk_R.csv  (shrunken cross-tissue correlation)
#       - *_Adjacency_beta<beta>.csv  (|r|^beta for TOM/WGCNA)
#
# Edit the three PATHS below to your files as needed.

# -------- EDIT THESE PATHS IF NEEDED --------
file_EsophagusMuscularis <- "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_EsophagusMuscularis.csv"
file_EsophagusMucosa     <- "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_EsophagusMucosa.csv"
file_Liver               <- "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_Liver.csv"

beta <- 2  # soft-threshold power for adjacency (you can change)

# -------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("ashr", quietly = TRUE)) devtools::install_github("stephenslab/ashr")
if (!requireNamespace("CorShrink", quietly = TRUE)) devtools::install_github("kkdey/CorShrink")

suppressPackageStartupMessages({
  library(CorShrink)
})

# Utility: robust CSV loader (first column = sample IDs; others = genes)
load_tissue <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  if (ncol(df) < 2) stop(paste("File has <2 columns:", path))
  colnames(df)[1] <- "sample"
  # drop duplicate samples, keep first
  if (anyDuplicated(df$sample)) {
    dup <- df$sample[duplicated(df$sample)]
    message("Dropping duplicated samples in ", path, ": ", paste(unique(dup), collapse=", "))
    df <- df[!duplicated(df$sample), ]
  }
  rownames(df) <- df$sample
  df$sample <- NULL
  # coerce to numeric where possible
  for (j in seq_len(ncol(df))) {
    df[[j]] <- suppressWarnings(as.numeric(df[[j]]))
  }
  as.matrix(df)
}

# Build CorShrink input matrix Donors × (A:genes + B:genes) on COMMON genes only
build_pair_matrix <- function(A, B, nameA="A", nameB="B") {
  genes <- intersect(colnames(A), colnames(B))
  if (length(genes) < 2) stop("Need >=2 common genes between tissues")
  A2 <- A[, genes, drop=FALSE]
  B2 <- B[, genes, drop=FALSE]
  donors_union <- sort(unique(c(rownames(A2), rownames(B2))))
  X <- matrix(NA_real_, nrow=length(donors_union), ncol=2*length(genes),
              dimnames=list(donors_union, c(paste0(nameA, ":", genes), paste0(nameB, ":", genes))))
  # fill blocks
  X[rownames(A2), paste0(nameA, ":", genes)] <- A2
  X[rownames(B2), paste0(nameB, ":", genes)] <- B2
  X
}

# Process a pair: builds X, runs CorShrinkData, writes outputs
process_pair <- function(A, B, nameA, nameB, out_prefix, beta=2) {
  X <- build_pair_matrix(A, B, nameA, nameB)
  cs <- CorShrinkData(X)
  cor_mat <- if (is.list(cs)) {
    if (!is.null(cs$cor)) cs$cor else if (!is.null(cs$ash_cor)) cs$ash_cor else stop("Unexpected CorShrinkData output")
  } else as.matrix(cs)

  ixA <- grepl(paste0("^", nameA, ":"), colnames(cor_mat))
  ixB <- grepl(paste0("^", nameB, ":"), colnames(cor_mat))
  R_ab <- cor_mat[ixA, ixB, drop = FALSE]
  A_ab <- abs(R_ab)^beta

  write.csv(R_ab, paste0(out_prefix, "_CorShrunk_R.csv"), row.names = TRUE)
  write.csv(A_ab, paste0(out_prefix, "_Adjacency_beta", beta, ".csv"), row.names = TRUE)

  invisible(list(R_cross=R_ab, A_cross=A_ab))
}

# Load all three tissues
A <- load_tissue(file_EsophagusMuscularis)
B <- load_tissue(file_EsophagusMucosa)
C <- load_tissue(file_Liver)

# Choose an output directory (same as first file)
out_dir <- dirname(file_EsophagusMuscularis)

# Process all three pairs
pairs <- list(
  list(A, B, "A", "B", file.path(out_dir, "EsophagusMuscularis_vs_EsophagusMucosa")),
  list(A, C, "A", "B", file.path(out_dir, "EsophagusMuscularis_vs_Liver")),
  list(B, C, "A", "B", file.path(out_dir, "EsophagusMucosa_vs_Liver"))
)

for (p in pairs) {
  message("Processing ", basename(p[[5]]), " ...")
  process_pair(A=p[[1]], B=p[[2]], nameA=p[[3]], nameB=p[[4]], out_prefix=p[[5]], beta=beta)
  message("Done.")
}

message("All pairs processed. Outputs written to: ", out_dir)

# (Optional) Within-tissue blocks:
# To get within-tissue shrunken correlations (e.g., A×A), you can run:
# cs_within_A <- CorShrinkData(A); R_AA <- if (is.list(cs_within_A)) { if (!is.null(cs_within_A$cor)) cs_within_A$cor else cs_within_A$ash_cor } else as.matrix(cs_within_A)
# write.csv(R_AA, file.path(out_dir, "EsophagusMuscularis_within_CorShrunk_R.csv"))
# and similarly for B, C.
