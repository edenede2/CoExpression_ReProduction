#!/usr/bin/env Rscript

# =====================================================================
# mdc_run_cli.R — One-shot MDC with sample & gene permutations (CLI)
# =====================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Package 'optparse' is required. Install with install.packages('optparse')")
  }
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required. Install with install.packages('WGCNA') or BiocManager::install('WGCNA')")
  }
  library(optparse)
  library(WGCNA)
})

# ---------- Helpers ----------
.extract_donor <- function(x) sub("^([^-]+-[^-]+).*", "\\1", x)

read_tissue_expr <- function(path, tissue, check.names=FALSE) {
  X <- read.csv(path, check.names = check.names)
  if (ncol(X) < 2) stop("Bad CSV: expect first column=sample IDs, others=genes: ", path)
  rownames(X) <- X[,1]
  X <- as.matrix(X[,-1, drop=FALSE])
  colnames(X) <- paste0(tissue, "_", colnames(X))  # Tissue_Gene
  storage.mode(X) <- "double"
  X
}

aggregate_by_donor <- function(M) {
  d <- .extract_donor(rownames(M))
  idx <- split(seq_len(nrow(M)), d)
  A <- do.call(rbind, lapply(idx, function(ix) colMeans(M[ix, , drop=FALSE], na.rm=TRUE)))
  rownames(A) <- names(idx)
  A
}

assemble_cohort_matrix <- function(files_by_tissue, check.names=FALSE) {
  mats <- lapply(names(files_by_tissue), function(t) {
    M <- read_tissue_expr(files_by_tissue[[t]], t, check.names = check.names)
    aggregate_by_donor(M)
  })
  names(mats) <- names(files_by_tissue)
  all_donors <- sort(unique(unlist(lapply(mats, rownames))))
  all_genes  <- sort(unique(unlist(lapply(mats, colnames))))
  Z <- matrix(NA_real_, nrow=length(all_donors), ncol=length(all_genes),
              dimnames=list(all_donors, all_genes))
  for (t in names(mats)) {
    Z[rownames(mats[[t]]), colnames(mats[[t]])] <- mats[[t]]
  }
  Z
}

permuteVect <- function(v) sample(v, length(v), replace=FALSE)
lower_mean  <- function(M) { lt <- lower.tri(M); mean(abs(M[lt]), na.rm=TRUE) }

compute_MDC <- function(exprB, exprA, gene_idx, power=6) {
  if (length(gene_idx) < 2) return(NA_real_)
  XB <- exprB[, gene_idx, drop=FALSE]
  XA <- exprA[, gene_idx, drop=FALSE]
  corB <- abs(cor(XB, use="pairwise.complete.obs"))^power; diag(corB) <- 0
  corA <- abs(cor(XA, use="pairwise.complete.obs"))^power; diag(corA) <- 0
  lower_mean(corB) / lower_mean(corA)
}

empirical_p <- function(obs, null_vec) {
  (1 + sum(null_vec >= obs, na.rm=TRUE)) / (1 + sum(!is.na(null_vec)))
}

parse_named_pairs <- function(vec) {
  # vec like c("Adipose=/path/a.csv", "Muscle=/path/b.csv")
  if (is.null(vec) || length(vec) == 0) return(character())
  kv <- strsplit(vec, "=", fixed = TRUE)
  ok <- vapply(kv, function(x) length(x) == 2 && nzchar(x[1]) && nzchar(x[2]), logical(1))
  if (!all(ok)) stop("Bad --old/--young entry. Use Tissue=/full/path.csv")
  vals <- vapply(kv, function(x) x[2], character(1))
  nms  <- vapply(kv, function(x) x[1], character(1))
  names(vals) <- nms
  vals
}

# ---------- CLI options ----------
option_list <- list(
  make_option(c("-c", "--clusters"), type = "character", help = "Clusters table path (TS/CT output)."),
  make_option(c("-o", "--old"), type = "character", action = "append",
              help = "Old cohort tissue mapping: Tissue=/path/to/old.csv (repeat per tissue)"),
  make_option(c("-y", "--young"), type = "character", action = "append",
              help = "Young cohort tissue mapping: Tissue=/path/to/young.csv (repeat per tissue)"),
  make_option(c("--out_dir"), type = "character", default = "MDC_out", help = "Output directory [default %default]"),
  make_option(c("--run_name"), type = "character", default = "young_vs_old", help = "Run name prefix [default %default]"),
  make_option(c("--power"), type = "integer", default = 6, help = "Soft-threshold power for correlations [default %default]"),
  make_option(c("--perms"), type = "integer", default = 100, help = "Number of permutations per null [default %default]"),
  make_option(c("--alpha"), type = "double", default = 0.10, help = "Alpha for marking significance [default %default]"),
  make_option(c("--topK"), type = "integer", default = 12, help = "Number of top modules for heatmaps [default %default]"),
  make_option(c("--moduleMode"), type = "character", default = "numeric", help = "Module mode: numeric|color [default %default]"),
  make_option(c("--check_names"), action = "store_true", default = FALSE, help = "Use check.names=TRUE when reading CSVs")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------- Validate inputs ----------
if (is.null(opt$clusters)) stop("--clusters is required")
files_old   <- parse_named_pairs(opt$old)
files_young <- parse_named_pairs(opt$young)
if (length(files_old) == 0 || length(files_young) == 0) stop("Provide at least one --old and --young mapping")
if (!identical(sort(names(files_old)), sort(names(files_young)))) {
  stop("Tissue names in --old and --young must match.")
}

out_dir  <- opt$out_dir
run_name <- opt$run_name
corrlpower <- opt$power
no.perms <- opt$perms
alpha <- opt$alpha
moduleMode <- match.arg(tolower(opt$moduleMode), c("numeric","color"))

# ---------- IO setup ----------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)

# ---------- 1) Build cohort matrices ----------
message("Assembling cohort matrices ...")
expr_old   <- assemble_cohort_matrix(files_old,   check.names = opt$check_names)
expr_young <- assemble_cohort_matrix(files_young, check.names = opt$check_names)

common <- intersect(colnames(expr_old), colnames(expr_young))
expr_old   <- expr_old[,   common, drop=FALSE]
expr_young <- expr_young[, common, drop=FALSE]
message(sprintf("OLD: %d donors × %d genes | YOUNG: %d donors × %d genes | Common genes: %d",
                nrow(expr_old), ncol(expr_old), nrow(expr_young), ncol(expr_young), length(common)))

# ---------- 2) Map modules from clusters table ----------
clu <- read.delim(opt$clusters, header=TRUE, check.names=FALSE)
required_cols <- c("Cluster ID", "Tissue", "Gene Symbol")
miss <- setdiff(required_cols, colnames(clu))
if (length(miss) > 0) stop("clusters table missing columns: ", paste(miss, collapse=", "))
clu$GeneID <- paste0(clu$Tissue, "_", clu$`Gene Symbol`)

if (moduleMode == "numeric") {
  gene2mod <- setNames(as.integer(clu$`Cluster ID`), clu$GeneID)
  modules  <- gene2mod[colnames(expr_young)]
  modtb    <- table(modules, useNA = "no")
  mod_ids  <- as.integer(names(modtb))
  if (length(mod_ids) == 0) stop("No modules overlap the expression columns.")
  genes_idx_of_module <- function(mid) which(modules == mid)
} else {
  # color mode
  clu$ModuleColor <- labels2colors(as.integer(clu$`Cluster ID`))
  module_of_gene  <- setNames(clu$ModuleColor, clu$GeneID)
  modules <- module_of_gene[colnames(expr_young)]
  modules[is.na(modules)] <- "grey"
  modtb <- table(modules)
  modtb <- modtb[names(modtb) != "grey"]
  mod_ids <- names(modtb)
  if (length(mod_ids) == 0) stop("No non-grey modules overlap the expression columns.")
  genes_idx_of_module <- function(mid) which(modules == mid)
}

# ---------- 3) True MDC per module ----------
message("Computing true MDC per module ...")
yfold <- sapply(mod_ids, function(mid) {
  idx <- genes_idx_of_module(mid)
  compute_MDC(expr_young, expr_old, idx, power = corrlpower)
})
names(yfold) <- as.character(mod_ids)

# ---------- 4) Sample permutations ----------
message(sprintf("Sample permutations × %d ...", no.perms))
meanFoldChange <- replicate(no.perms, {
  XBY <- apply(expr_young, 2, permuteVect)
  XAO <- apply(expr_old,   2, permuteVect)
  sapply(mod_ids, function(mid) {
    idx <- genes_idx_of_module(mid)
    compute_MDC(XBY, XAO, idx, power = corrlpower)
  })
})
if (is.numeric(mod_ids)) rownames(meanFoldChange) <- as.character(mod_ids)
colnames(meanFoldChange) <- sprintf("MDC_random_samples_%02d", seq_len(no.perms))

# ---------- 5) Gene permutations ----------
message(sprintf("Gene permutations × %d ...", no.perms))
all_genes <- seq_len(ncol(expr_young))
GlobalyRandomMDC <- replicate(no.perms, {
  sapply(mod_ids, function(mid) {
    msize <- as.integer(modtb[as.character(mid)])
    if (is.na(msize) || msize < 2) return(NA_real_)
    idx <- sample(all_genes, msize, replace = FALSE)
    compute_MDC(expr_young, expr_old, idx, power = corrlpower)
  })
})
if (is.numeric(mod_ids)) rownames(GlobalyRandomMDC) <- as.character(mod_ids)
colnames(GlobalyRandomMDC) <- sprintf("MDC_random_genes_%02d", seq_len(no.perms))

# ---------- 6) p-values & FDR ----------
p_samp <- vapply(seq_along(yfold), function(i) {
  empirical_p(yfold[i], meanFoldChange[ names(yfold)[i], ])
}, numeric(1))
names(p_samp) <- names(yfold)

p_gene <- vapply(seq_along(yfold), function(i) {
  empirical_p(yfold[i], GlobalyRandomMDC[ names(yfold)[i], ])
}, numeric(1))
names(p_gene) <- names(yfold)

q_samp <- p.adjust(p_samp, method = "BH")
q_gene <- p.adjust(p_gene, method = "BH")
q_max  <- pmax(q_samp, q_gene)

# ---------- 7) Save tables ----------
mod_ids_char <- names(yfold)
final <- data.frame(
  module = if (is.numeric(mod_ids)) as.integer(mod_ids_char) else mod_ids_char,
  size   = as.integer(modtb[mod_ids_char]),
  MDC    = round(unname(yfold), 4),
  p_random_samples = signif(unname(p_samp), 4),
  q_random_samples = signif(unname(q_samp), 4),
  p_random_genes   = signif(unname(p_gene), 4),
  q_random_genes   = signif(unname(q_gene), 4),
  q_combined_max   = signif(unname(q_max), 4),
  check.names = FALSE
)
final <- final[order(-final$MDC), ]

stopifnot(all(rownames(meanFoldChange)   %in% mod_ids_char) || !is.numeric(mod_ids))
stopifnot(all(rownames(GlobalyRandomMDC) %in% mod_ids_char) || !is.numeric(mod_ids))

randoms <- data.frame(
  module = if (is.numeric(mod_ids)) as.integer(mod_ids_char) else mod_ids_char,
  as.data.frame(meanFoldChange[mod_ids_char, , drop=FALSE], check.names = FALSE),
  as.data.frame(GlobalyRandomMDC[mod_ids_char, , drop=FALSE], check.names = FALSE),
  check.names = FALSE
)
colnames(randoms) <- c(
  "module",
  sprintf("MDC_random_samples_%02d", seq_len(ncol(meanFoldChange))),
  sprintf("MDC_random_genes_%02d",  seq_len(ncol(GlobalyRandomMDC)))
)

write.table(final,
            file = file.path(out_dir, paste0(run_name, "_MDC_wFDR.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(randoms,
            file = file.path(out_dir, paste0(run_name, "_MDC_randoms.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: ", file.path(out_dir, paste0(run_name, "_MDC_wFDR.tsv")))
message("Saved: ", file.path(out_dir, paste0(run_name, "_MDC_randoms.tsv")))

# ---------- 8) Simple barplot ----------
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
png(file.path(out_dir, "plots", paste0(run_name, "_MDC_barplot.png")),
    width = 2000, height = 900, res = 180)
op <- par(mar = c(9, 5, 3, 1))
lab <- final$module
cols <- rep_len(setdiff(colors(), c("white","grey","gray","lightgray")), length(lab))
stars <- ifelse(final$q_combined_max <= alpha, "*", "")
bp <- barplot(final$MDC, names.arg = lab, las = 2,
              ylab = sprintf("MDC (k_%s / k_%s)", "young", "old"),
              col = cols[seq_along(lab)], border = NA)
text(bp, final$MDC, labels = stars, pos = 3, cex = 1.1)
abline(h = 1, lty = 2)
par(op)
dev.off()
message("Saved: ", file.path(out_dir, "plots", paste0(run_name, "_MDC_barplot.png")))

# ---------- 9) Optional: per-module heatmaps (top K) ----------
dir.create(file.path(out_dir, "plots", "heatmaps"), showWarnings = FALSE)
topK <- min(as.integer(opt$topK), nrow(final))
for (z in seq_len(topK)) {
  mid <- final$module[z]
  idx <- genes_idx_of_module(mid)
  if (length(idx) < 2) next
  corB <- abs(cor(expr_young[, idx, drop=FALSE], use="pairwise.complete.obs"))^corrlpower
  corA <- abs(cor(expr_old  [, idx, drop=FALSE], use="pairwise.complete.obs"))^corrlpower
  diag(corB) <- 0; diag(corA) <- 0
  ord <- order(rowSums(corB, na.rm=TRUE), decreasing = TRUE)
  CB <- corB[ord, ord]; CA <- corA[ord, ord]
  U  <- CB; L <- CA
  U[lower.tri(U, diag=TRUE)] <- NA
  L[upper.tri(L, diag=TRUE)] <- NA
  M <- U; M[is.na(M)] <- L[is.na(M)]
  png(file.path(out_dir, "plots", "heatmaps",
                sprintf("%s_heatmap_module%s.png", run_name, as.character(mid))),
      width = 800, height = 800, res = 140)
  par(mar=c(2,2,2,2))
  image(t(apply(M, 2, rev)), axes=FALSE, col = heat.colors(64),
        main = sprintf("Module %s (top %d/%d genes)", as.character(mid), length(idx), as.integer(modtb[as.character(mid)])))
  box()
  dev.off()
}
message("Done MDC pipeline.")

