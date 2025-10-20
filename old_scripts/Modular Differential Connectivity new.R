# =====================================================================
# MDC_run.R  —  One-shot MDC with sample & gene permutations (standalone)
# =====================================================================

# ---------- USER PARAMETERS (EDIT THESE) ----------
files_old <- c(
  Adipose="/Users/edeneldar/CoExpression_ReProduction/old_scripts/outputs/Adipose - Subcutaneous/old_matrix.csv",
  Muscle ="/Users/edeneldar/CoExpression_ReProduction/old_scripts/outputs/Muscle - Skeletal/old_matrix.csv",
  Brain  ="/Users/edeneldar/CoExpression_ReProduction/old_scripts/outputs/Brain - Cortex/old_matrix.csv"
)

files_old_shaked <- c(
  Adipose="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Adipose - Subcutaneous_old.csv",
  Muscle ="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Muscle - Skeletal_old.csv",
  Brain  ="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Brain - Cortex_old.csv"
)

files_young <- c(
  Adipose="/Users/edeneldar/CoExpression_ReProduction/old_scripts/outputs/Adipose - Subcutaneous/young_matrix.csv",
  Muscle ="/Users/edeneldar/CoExpression_ReProduction/old_scripts/outputs/Muscle - Skeletal/young_matrix.csv",
  Brain  ="/Users/edeneldar/CoExpression_ReProduction/old_scripts/outputs/Brain - Cortex/young_matrix.csv"
)

files_young_shaked <- c(
  Adipose="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Adipose - Subcutaneous_young.csv",
  Muscle ="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Muscle - Skeletal_young.csv",
  Brain  ="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Brain - Cortex_young.csv"
)
clusters_table_path_shaked <- "/Users/edeneldar/CoExpression_ReProduction/xwgcna_young_original_run9_Clusters_table.txt"  # columns: Cluster ID,Tissue,Gene Symbol 
clusters_table_path <- "xwgcna_young_auto_run1_Clusters_table.txt"  # columns: Cluster ID,Tissue,Gene Symbol
# clu <- read.delim("xwgcna_young_auto_run1_Clusters_table.txt", check.names = FALSE)  

out_dir     <- "MDC_out"
run_name    <- "young_vs_old"
set.seed(42)

# Analysis knobs
corrlpower  <- 6          # like WGCNA soft-threshold (β) used for within-cohort correlation
no.perms    <- 100         # increase to 100/200 for stronger nulls
alpha       <- 0.10       # significance threshold for bars

# ---------- Packages ----------
suppressPackageStartupMessages({
  library(WGCNA)   # for labels2colors
})

# ---------- Helpers ----------
.extract_donor <- function(x) sub("^([^-]+-[^-]+).*", "\\1", x)

read_tissue_expr <- function(path, tissue, check.names=FALSE) {
  X <- read.csv(path, check.names = check.names)
  if (ncol(X) < 2) stop("Bad CSV: expect first column=sample IDs, others=genes")
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

assemble_cohort_matrix <- function(files_by_tissue) {
  mats <- lapply(names(files_by_tissue), function(t) {
    M <- read_tissue_expr(files_by_tissue[[t]], t)
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
  # one-sided: large MDC supports difference (≥)
  (1 + sum(null_vec >= obs, na.rm=TRUE)) / (1 + sum(!is.na(null_vec)))
}

# ---------- IO setup ----------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)

# ---------- 1) Build cohort matrices ----------
message("Assembling cohort matrices ...")
expr_old   <- assemble_cohort_matrix(files_old_shaked)
expr_young <- assemble_cohort_matrix(files_young_shaked)
message(sprintf("OLD: %d donors × %d genes | YOUNG: %d donors × %d genes",
                nrow(expr_old), ncol(expr_old), nrow(expr_young), ncol(expr_young)))

# Align columns to common genes (Tissue_Gene)
common <- intersect(colnames(expr_old), colnames(expr_young))
expr_old   <- expr_old[,   common, drop=FALSE]
expr_young <- expr_young[, common, drop=FALSE]
message(sprintf("OLD: %d donors × %d genes | YOUNG: %d donors × %d genes | Common genes: %d",
                nrow(expr_old), ncol(expr_old), nrow(expr_young), ncol(expr_young), length(common)))

# ---------- 2) Map modules by original Cluster ID (numeric) ----------
message("Mapping modules from clusters table ...")
clu <- read.delim(clusters_table_path_shaked, header=TRUE, check.names=FALSE)
required_cols <- c("Cluster ID", "Tissue", "Gene Symbol")
miss <- setdiff(required_cols, colnames(clu))
if (length(miss) > 0) stop("clusters table missing columns: ", paste(miss, collapse=", "))

message(sprintf("Clusters table: %d rows × %d columns", nrow(clu), ncol(clu)))
clu$GeneID <- paste0(clu$Tissue, "_", clu$`Gene Symbol`)
message(sprintf("Found %d unique genes in clusters table", length(unique(clu$GeneID))))
gene2modID <- setNames(as.integer(clu$`Cluster ID`), clu$GeneID)
message(sprintf("Found %d unique modules in clusters table", length(unique(gene2modID))))


modulesID_young <- gene2modID[colnames(expr_young)]
message(sprintf("Found %d modules in young cohort: %s",
                length(unique(modulesID_young)), paste(unique(modulesID_young), collapse=", ")))

modtb <- table(modulesID_young, useNA = "no")
mod_ids <- as.integer(names(modtb))       
if (length(mod_ids) == 0) stop("No modules overlap the expression columns.")
message(sprintf("Found %d modules in young cohort: %s",
                length(mod_ids), paste(mod_ids, collapse=", ")))
genes_idx_of_module <- function(mid) which(modulesID_young == mid)
message(sprintf("Modules sizes: %s", paste(modtb, collapse=", ")))
gene_tissue <- sub("_.*$", "", colnames(expr_young))

message("Modules by tissue:")
for (t in sort(unique(gene_tissue))) {
  idx <- which(gene_tissue == t)
  msize <- sum(modulesID_young[idx] %in% mod_ids)
  message(sprintf("  %s: %d genes in %d modules", t, length(idx), msize))
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
rownames(meanFoldChange) <- as.character(mod_ids)
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
rownames(GlobalyRandomMDC) <- as.character(mod_ids)
colnames(GlobalyRandomMDC) <- sprintf("MDC_random_genes_%02d", seq_len(no.perms))

# ---------- 6) p-values & FDR ----------
# empirical p per module from each permutation family
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
  module = as.integer(mod_ids_char),
  size   = as.integer(modtb[mod_ids_char]),
  MDC    = round(unname(yfold), 10),
  p_random_samples = signif(unname(p_samp), 10),
  q_random_samples = signif(unname(q_samp), 10),
  p_random_genes   = signif(unname(p_gene), 10),
  q_random_genes   = signif(unname(q_gene), 10),
  q_combined_max   = signif(unname(q_max), 10),
  check.names = FALSE
)
final <- final[order(-final$MDC), ]

stopifnot(all(rownames(meanFoldChange)   %in% mod_ids_char))
stopifnot(all(rownames(GlobalyRandomMDC) %in% mod_ids_char))

# Rebuild randoms properly (avoid transposition confusion):
randoms <- data.frame(
  module = as.integer(mod_ids_char),
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
png(file.path(out_dir, "plots", paste0(run_name, "_MDC_barplot.png")),
    width = 2000, height = 900, res = 180)
op <- par(mar = c(9, 5, 3, 1))
lab <- final$module
cols <- rep_len(setdiff(colors(), c("white","grey","gray","lightgray")), length(lab))
stars <- ifelse(final$q_random_genes <= alpha, "*", "")
bp <- barplot(final$MDC, names.arg = lab, las = 2,
              ylab = sprintf("MDC (k_%s / k_%s)", "young", "old"),
              col = cols[seq_along(lab)], border = NA)
text(bp, final$MDC, labels = stars, pos = 3, cex = 1.1)
abline(h = 1, lty = 2)
par(op)
dev.off()
message("Saved: ", file.path(out_dir, "plots", paste0(run_name, "_MDC_barplot.png")))

# ---------- 9) Optional: per-module heatmaps (top K) ----------
topK <- min(12L, nrow(final))
dir.create(file.path(out_dir, "plots", "heatmaps"), showWarnings = FALSE)
# Order modules by MDC already in 'final'
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
                sprintf("%s_heatmap_module%d.png", run_name, mid)),
      width = 800, height = 800, res = 140)
  par(mar=c(2,2,2,2))
  image(t(apply(M, 2, rev)), axes=FALSE, col = heat.colors(64),
        main = sprintf("Module %d (top %d/%d genes)", mid, length(idx), as.integer(modtb[as.character(mid)])))
  box()
  dev.off()
}
message("Done numeric module MDC pipeline.")
# =====================================================================



# ---------- 2) Map modules from clusters table (colors) ----------
clu <- read.delim(clusters_table_path_shaked, header=TRUE, check.names=FALSE)
# Expect columns: "Cluster ID", "Tissue", "Gene Symbol"
required_cols <- c("Cluster ID", "Tissue", "Gene Symbol")
miss <- setdiff(required_cols, colnames(clu))
if (length(miss) > 0) stop("clusters table missing columns: ", paste(miss, collapse=", "))

clu$GeneID <- paste0(clu$Tissue, "_", clu$`Gene Symbol`)
# map Cluster ID -> color (like WGCNA)
clu$ModuleColor <- labels2colors(as.integer(clu$`Cluster ID`))
module_of_gene  <- setNames(clu$ModuleColor, clu$GeneID)

modulescolorB <- module_of_gene[colnames(expr_young)]
modulescolorB[is.na(modulescolorB)] <- "grey"
modtb <- table(modulescolorB)
modtb <- modtb[names(modtb) != "grey"]
modulenames <- names(modtb)
if (length(modulenames) == 0) stop("No non-grey modules overlap the expression columns.")

# ---------- 3) True MDC per module ----------
message("Computing true MDC per module ...")
yfold <- sapply(modulenames, function(mcol) {
  idx <- which(modulescolorB == mcol)
  compute_MDC(expr_young, expr_old, idx, power = corrlpower)
})

# ---------- 4) Sample permutations ----------
message(sprintf("Sample permutations × %d ...", no.perms))
meanFoldChange <- replicate(no.perms, {
  XBY <- apply(expr_young, 2, permuteVect)  # permute donors independently per gene
  XAO <- apply(expr_old,   2, permuteVect)
  sapply(modulenames, function(mcol) {
    idx <- which(modulescolorB == mcol)
    compute_MDC(XBY, XAO, idx, power = corrlpower)
  })
})
colnames(meanFoldChange) <- sprintf("MDC_random_samples_%02d", seq_len(no.perms))

# ---------- 5) Gene permutations ----------
message(sprintf("Gene permutations × %d ...", no.perms))
all_genes <- seq_len(ncol(expr_young))
GlobalyRandomMDC <- replicate(no.perms, {
  sapply(modulenames, function(mcol) {
    msize <- as.integer(modtb[mcol])
    if (msize < 2) return(NA_real_)
    idx <- sample(all_genes, msize, replace = FALSE)
    compute_MDC(expr_young, expr_old, idx, power = corrlpower)
  })
})
colnames(GlobalyRandomMDC) <- sprintf("MDC_random_genes_%02d", seq_len(no.perms))

# ---------- 6) p-values & FDR ----------
# empirical p per module from each permutation family
p_samp <- mapply(function(obs, i) empirical_p(obs, meanFoldChange[i, ]), obs = yfold, i = seq_along(yfold))
p_gene <- mapply(function(obs, i) empirical_p(obs, GlobalyRandomMDC[i, ]), obs = yfold, i = seq_along(yfold))

# BH FDR
q_samp <- p.adjust(p_samp, method = "BH")
q_gene <- p.adjust(p_gene, method = "BH")
q_max  <- pmax(q_samp, q_gene)   # conservative combination (like taking max of two FDRs)

# ---------- 7) Save tables ----------
final <- data.frame(
  module = modulenames,
  size   = as.integer(modtb[modulenames]),
  MDC    = round(yfold, 4),
  p_random_samples = signif(p_samp, 4),
  q_random_samples = signif(q_samp, 4),
  p_random_genes   = signif(p_gene, 4),
  q_random_genes   = signif(q_gene, 4),
  q_combined_max   = signif(q_max, 4),
  check.names = FALSE
)
final <- final[order(-final$MDC), ]

stopifnot(nrow(meanFoldChange)    == length(modulenames))
stopifnot(nrow(GlobalyRandomMDC)  == length(modulenames))

randoms <- cbind(
    module = modulenames,
    setNames(as.data.frame(meanFoldChange,    check.names = FALSE),
                   sprintf("MDC_random_samples_%02d", seq_len(ncol(meanFoldChange)))),
    setNames(as.data.frame(GlobalyRandomMDC,  check.names = FALSE),
                   sprintf("MDC_random_genes_%02d",  seq_len(ncol(GlobalyRandomMDC))))
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
png(file.path(out_dir, "plots", paste0(run_name, "_MDC_barplot.png")),
    width = 2000, height = 900, res = 180)
op <- par(mar = c(9, 5, 3, 1))
lab <- final$module
cols <- rep_len(setdiff(colors(), c("white","grey","gray","lightgray")), length(lab))
# mark significance by combined FDR
stars <- ifelse(final$q_combined_max <= alpha, "*", "")
bp <- barplot(final$MDC, names.arg = lab, las = 2,
              ylab = sprintf("MDC (k_%s / k_%s)", "young", "old"),
              col = cols[seq_along(lab)], border = NA)
text(bp, final$MDC, labels = stars, pos = 3, cex = 1.2)
abline(h = 1, lty = 2)
par(op)
dev.off()
message("Saved: ", file.path(out_dir, "plots", paste0(run_name, "_MDC_barplot.png")))

# ---------- 9) Optional: per-module heatmaps (top K) ----------
topK <- min(12L, nrow(final))
dir.create(file.path(out_dir, "plots", "heatmaps"), showWarnings = FALSE)
for (z in seq_len(topK)) {
  mcol <- final$module[z]
  idx  <- which(modulescolorB == mcol)
  if (length(idx) < 2) next
  corB <- abs(cor(expr_young[, idx, drop=FALSE], use="pairwise.complete.obs"))^corrlpower
  corA <- abs(cor(expr_old  [, idx, drop=FALSE], use="pairwise.complete.obs"))^corrlpower
  diag(corB) <- 0; diag(corA) <- 0
  
  # combine into one tile: upper = B, lower = A (like the original idea)
  od <- order(rowSums(corB, na.rm=TRUE), decreasing = TRUE)
  CB <- corB[od, od]; CA <- corA[od, od]
  U  <- CB; L <- CA
  U[lower.tri(U, diag=TRUE)] <- NA
  L[upper.tri(L, diag=TRUE)] <- NA
  M <- U; M[is.na(M)] <- L[is.na(M)]
  
  png(file.path(out_dir, "plots", "heatmaps",
                sprintf("%s_heatmap_%s.png", run_name, mcol)),
      width = 800, height = 800, res = 140)
  par(mar=c(2,2,2,2))
  image(t(apply(M, 2, rev)), axes=FALSE, col = heat.colors(64),
        main = sprintf("Module %s (top %d/%d genes)", mcol, length(idx), as.integer(modtb[mcol])))
  box()
  dev.off()
}
message("Done.")
# =====================================================================
