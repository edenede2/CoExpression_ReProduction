# Modular Differential Connectivity (MDC): 
# Compute ratio of the mean intra-module connectivity in one network to that of the same set of genes in another network
# This code uses functions from: https://github.com/songw01/melanoma_network_codes
# Author(s): Bin Zhang Cell 153(3):707-720. PMID: 23622250 
#
#
#    Input:
#         1) inputfnameB       -- mRNA data (tab-delimited text file) for one disease state B
#         2) inputfnameBModule -- module assignment file as output from WINA or WGCNA
#         3) inputfname        -- mRNA data (tab-delimited text file) for another disease state A
#         4) shortnames        -- shortnames for B and A, respectively. Order is important
#         5) headCol           -- # of gene information columns in inputfname
#         6) headCol2          -- # of gene information columns in inputfnameB
#         7) corrlpower        -- the exponent of the power function, determined by WINA or WGCNA
#
#   Output:
#         1) "DiffConnectvty\*_wFDR.xls"             -- MDC with false discovery rates
#         2) "DiffConnectvty\*.png"                  -- MDC plot
#         3) "DiffConnectvty\Heatmaps\*_iHMP0_*.png" -- MDC heatmaps of mulitple modules 
#         3) "DiffConnectvty\Heatmaps\*_iHMP-*.png" -- MDC heatmaps of individual modules
#
#
################################################################################################
#memory.limit(size = 35000)
rm(list = ls())
library(lattice) # require is design for use inside functions 
library(plotrix)
#library(sma) # this is needed for plot.mat below
source("/Users/edeneldar/CoExpression_ReProduction/old_scripts/R-functions_MDC.R")
source("/Users/edeneldar/CoExpression_ReProduction/old_scripts/colPalette_200.txt")


##### 
inputfnameB       = "/Users/edeneldar/CoExpression_ReProduction/xwgcna_young_original_run9_adjacency.rds"; headCol2 =2;
inputfnameBModule = "/Users/edeneldar/CoExpression_ReProduction/xwgcna_young_original_run9_Clusters_table.txt"
inputfname        = "/Users/edeneldar/CoExpression_ReProduction/xwgcna_old_original_run9_adjacency.rds"; headCol =2;
shortnames        = c("young", "old")
corrlpower        = 6

# MDC powers for expression-based MDC computation (sample permutations only)
# Apply TS (within-tissue) pairs with power 3 and CT (cross-tissue) pairs with power 6
if (!exists("MDC_TS_POWER")) MDC_TS_POWER <- 3
if (!exists("MDC_CT_POWER")) MDC_CT_POWER <- 6


files_old_shaked <- c(
  Adipose="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Adipose - Subcutaneous_old.csv",
  Muscle ="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Muscle - Skeletal_old.csv",
  Brain  ="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Brain - Cortex_old.csv"
)


files_young_shaked <- c(
  Adipose="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Adipose - Subcutaneous_young.csv",
  Muscle ="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Muscle - Skeletal_young.csv",
  Brain  ="/Users/edeneldar/CoExpression_ReProduction/old_outputs/Brain - Cortex_young.csv"
)


# ------------------- helpers to assemble expression matrices -------------------
.extract_donor <- function(x) sub("^([^-]+-[^-]+).*", "\\1", x)

read_tissue_expr <- function(path, tissue, check.names=FALSE) {
  X <- read.csv(path, check.names = check.names)
  if (ncol(X) < 2) stop(paste0("Bad CSV for ", tissue, ": expect first column=sample IDs, others=genes"))
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

# ------------------- build expr matrices if file paths are provided -----------
if (exists("files_old_shaked") && exists("files_young_shaked")) {
  message("Assembling cohort matrices for sample permutations ...")
  expr_old   <- assemble_cohort_matrix(files_old_shaked)
  expr_young <- assemble_cohort_matrix(files_young_shaked)
  common <- intersect(colnames(expr_old), colnames(expr_young))
  if (length(common) >= 2) {
    expr_old   <- expr_old[,   common, drop=FALSE]
    expr_young <- expr_young[, common, drop=FALSE]
    message(sprintf("OLD: %d donors × %d genes | YOUNG: %d donors × %d genes | Common genes: %d",
                    nrow(expr_old), ncol(expr_old), nrow(expr_young), ncol(expr_young), length(common)))
  } else {
    warning("No common genes between old/young expression inputs; sample permutations will be skipped.")
  }
}


# User-tunable parameters (can be overridden before sourcing this file)
if (!exists("NO_PERMS")) NO_PERMS <- 50
if (!exists("PERMUTE_GENES")) PERMUTE_GENES <- TRUE
if (!exists("PERMUTE_SAMPLES")) PERMUTE_SAMPLES <- FALSE  
if (!exists("COMPUTE_PVALUES")) COMPUTE_PVALUES <- TRUE

# Finish and generalize permutation-based validation
validateMDC <- function(no_perms = NO_PERMS,
                        permute_genes = PERMUTE_GENES,
                        permute_samples = PERMUTE_SAMPLES,
                        get_pVal = NULL) {
  # Dependencies: uses datExpr, datExpr2, modtb, no.genes, no.modules, yfold from parent env
  # Returns: list(random_genes, random_samples, p_values)

  # Sanity checks
  if (!exists("datExpr") || !exists("datExpr2")) stop("datExpr/datExpr2 not found. Load adjacency matrices first.")
  if (!exists("modtb")) stop("modtb (module sizes) not found. Build modules first.")
  if (!exists("yfold")) warning("Observed MDC (yfold) not found yet; p-values will be NA unless computed later.")

  # Helper: robust lower-tri mean of absolute values
  lower_mean_local <- function(M) {
    if (is.null(M) || any(dim(M) == 0)) return(NA_real_)
    diag(M) <- 0
    m <- abs(M)
    mean(m[lower.tri(m)], na.rm = TRUE)
  }

  random_genes <- if (isTRUE(permute_genes)) matrix(NA_real_, nrow = length(modtb), ncol = no_perms) else NULL
  random_samples <- if (isTRUE(permute_samples)) matrix(NA_real_, nrow = length(modtb), ncol = no_perms) else NULL

  # Gene-label permutations: sample random gene sets matching each module size
  if (isTRUE(permute_genes)) {
    for (i in seq_len(no_perms)) {
      if (i %% 5 == 0) {
        print(paste("********* random genes", i, " ........"))
      }
      for (x in seq_along(modtb)) {
        k <- as.integer(modtb[x])
        if (is.na(k) || k < 2) {
          random_genes[x, i] <- NA_real_
          next
        }
        xsel <- sample.int(no.genes, k, replace = FALSE)
        A <- datExpr [xsel, xsel, drop = FALSE]  # old
        B <- datExpr2[xsel, xsel, drop = FALSE]  # young
        # Ratio consistent with yfold definition below (young/old)
        random_genes[x, i] <- lower_mean_local(B) / lower_mean_local(A)
      }
      collect_garbage()
    }
  }

  # Sample-label permutations using raw expression matrices (if available)
  if (isTRUE(permute_samples)) {
    # Expect global expr_young and expr_old matrices with donor rows and gene columns matching Tissue_Gene IDs
    if (!exists("expr_young") || !exists("expr_old")) {
      warning("permute_samples=TRUE, but expr_young/expr_old not found; skipping sample permutations.")
      random_samples <- NULL
    } else {
      # Align expression matrices to a common gene set and to module gene IDs
      genes_adj <- colnames(datExpr2)
      common <- intersect(intersect(colnames(expr_young), colnames(expr_old)), genes_adj)
      if (length(common) < 2) {
        warning("No overlapping genes between expression matrices and adjacency/module labels; skipping sample permutations.")
        random_samples <- NULL
      } else {
        EY <- expr_young[, common, drop = FALSE]
        EO <- expr_old  [, common, drop = FALSE]
        # Map module indices in expression space
        mod_idx_in_expr <- lapply(modulenames, function(mcol) {
          gnames <- genes_adj[modulescolorB == mcol]
          which(colnames(EY) %in% gnames)
        })

        # helpers
        permuteVect <- function(v) sample(v, length(v), replace = FALSE)
        compute_MDC_expr <- function(exprB, exprA, idx, ts_power = MDC_TS_POWER, ct_power = MDC_CT_POWER) {
          if (length(idx) < 2) return(NA_real_)
          XB <- exprB[, idx, drop = FALSE]
          XA <- exprA[, idx, drop = FALSE]
          # Determine tissue for each column (expects Tissue_Gene)
          gnames <- colnames(XB)
          tissues <- sub("_.*$", "", gnames)
          same_tissue <- outer(tissues, tissues, FUN = "==")
          # Correlations for young (B) and old (A)
          RB <- abs(stats::cor(XB, use = "pairwise.complete.obs")); diag(RB) <- 0
          RA <- abs(stats::cor(XA, use = "pairwise.complete.obs")); diag(RA) <- 0
          # Apply pair-specific powers: TS -> ts_power, CT -> ct_power
          CB <- RB
          CA <- RA
          CB[same_tissue] <- CB[same_tissue]^ts_power
          CB[!same_tissue] <- CB[!same_tissue]^ct_power
          CA[same_tissue] <- CA[same_tissue]^ts_power
          CA[!same_tissue] <- CA[!same_tissue]^ct_power
          lower_mean_local(CB) / lower_mean_local(CA)
        }

        for (i in seq_len(no_perms)) {
          if (i %% 5 == 0) print(paste("********* random samples", i, " ........"))
          XBY <- apply(EY, 2, permuteVect)
          XAO <- apply(EO, 2, permuteVect)
          for (x in seq_along(mod_idx_in_expr)) {
            idx <- mod_idx_in_expr[[x]]
            random_samples[x, i] <- compute_MDC_expr(XBY, XAO, idx, ts_power = MDC_TS_POWER, ct_power = MDC_CT_POWER)
          }
          collect_garbage()
        }
      }
    }
  }

  # Empirical p-values (one-sided; more connected in young → MDC > 1) unless a custom function is provided
  p_values <- NULL
  p_values_samples <- NULL
  if (isTRUE(COMPUTE_PVALUES) && exists("yfold") && !is.null(random_genes)) {
    if (is.null(get_pVal)) {
      # Default empirical p-values per module using gene permutations (upper-tail)
      p_values <- rep(NA_real_, length(yfold))
      for (r in seq_along(yfold)) {
        perms <- random_genes[r, ]
        perms <- perms[is.finite(perms)]
        if (length(perms) < 3 || !is.finite(yfold[r])) {
          p_values[r] <- NA_real_
        } else {
          # One-sided: probability permutation MDC >= observed MDC
          p_values[r] <- (1 + sum(perms >= yfold[r])) / (length(perms) + 1)
        }
      }
    } else {
      # Allow user to supply a custom p-value function
      # Expected signature: get_pVal(observed_vector, perm_matrix) -> vector p-values
      p_values <- tryCatch(get_pVal(yfold, random_genes), error = function(e) {
        warning(paste("get_pVal failed:", e$message)); NULL
      })
    }
  }
  if (isTRUE(COMPUTE_PVALUES) && exists("yfold") && !is.null(random_samples)) {
    if (is.null(get_pVal)) {
      p_values_samples <- rep(NA_real_, length(yfold))
      for (r in seq_along(yfold)) {
        perms <- random_samples[r, ]
        perms <- perms[is.finite(perms)]
        if (length(perms) < 3 || !is.finite(yfold[r])) {
          p_values_samples[r] <- NA_real_
        } else {
          # One-sided: probability permutation MDC >= observed MDC
          p_values_samples[r] <- (1 + sum(perms >= yfold[r])) / (length(perms) + 1)
        }
      }
    } else {
      p_values_samples <- tryCatch(get_pVal(yfold, random_samples), error = function(e) {
        warning(paste("get_pVal failed (samples):", e$message)); NULL
      })
    }
  }

  list(random_genes = random_genes, random_samples = random_samples, p_values = p_values, p_values_samples = p_values_samples)
}
#
# -----------------------------End of Parameters to be changed --------------------------------------

# Number of permutations used by this run
no.perms = NO_PERMS
imgwid=600; imghei=400

#heatmapColorRG = rev( rgcolors.func(50) )
heatmapColor = rev(heat.colors(100))

# specify the directory for holding analysis results, you need make such a sub-directory 
#    under your working directory if it doesn't exist
outputDir1  ="DiffConnectvty_young_to_old_new/"
dir.create(outputDir1)

outputDir  ="DiffConnectvty_young_to_old_new/Random/"
dir.create(outputDir)

outputDirTom = "DiffConnectvty_young_to_old_new/Heatmaps/"
dir.create(outputDirTom)

# image type
imagetype="png" #"pdf"


fname      =getFileNameNopath(inputfname)
fnameB     =getFileNameNopath(inputfnameBModule)
extname    =getFileExtension(inputfname)

fname = paste(fnameB, "_VS_", shortnames[2], "_MDC", sep="")
flog  = paste(outputDir1, fname, "_Randoms", ".xls",       sep='')
flog2 = paste(outputDir1, fname, "_wFDR", ".xls",       sep='')
fimgRt= paste(outputDir1, fname, ".png",sep='')  #png only

#------------------------- second file --------------------------------------

datExpr2 = readRDS(inputfnameB)
dim(datExpr2)

no.genes <- dim(datExpr2)[2]
genesInfor2 <- colnames(datExpr2)

# ------ modules from 2nd file --------------------------------
allMatrix2modinfo <- read.delim(inputfnameBModule,sep="\t", header=T)
allMatrix2modinfo$Gene.Symbol <- paste(allMatrix2modinfo$Tissue,"_",allMatrix2modinfo$Gene.Symbol,sep="")
dim(allMatrix2modinfo)
ncols2m = dim(allMatrix2modinfo)[2]

allMatrix2modinfo <- as.matrix(allMatrix2modinfo)[,c(1,ncols2m)]
genes2Idx = cbind(as.matrix(genesInfor2)[,1], c(1:no.genes) )
merged = merge(genes2Idx, allMatrix2modinfo, by.x=1, by.y=2, all.x=T)
merged = as.matrix(merged)
dim(merged)
morder = order(as.integer(merged[,2]))
allMatrix2modinfoXX = merged[morder,] # aligned to the other gene infos
modulescolorB = as.character(merged[morder,3])
modulescolorB = ifelse(is.na(modulescolorB), "grey", modulescolorB)

modtb = table(modulescolorB)
modulenames = names(modtb)
modtb = modtb[modulenames!="grey"]

if (length(modulenames) > length(col.names))
{
  col.names <- colors();
  col.names <- col.names[-grep("grey|gray|white",col.names)]
  if (length(modulenames) > length(col.names))
  {
    col.names <- rep(col.names,ceiling(length(modulenames)/length(col.names)))
  }
}


modulenames = names(modtb)

umodulesizes = union(sort(as.integer(modtb)), NULL)

no.modules  = length(modulenames)
no.nets     = 2

#------------------------- first file --------------------------------------

datExpr <- readRDS(inputfname)
dim(datExpr)

#corhelp<- abs(corhelp)

#*-------------------------------------------------------------------------------------
#* initilization

uno = length(umodulesizes)

meanPermodule   = matrix(0, no.modules, no.nets)
sdPermodule     = matrix(0, no.modules, no.nets)
ks.pvalues      = rep(1, no.modules)
linkCor         = matrix(-1, no.modules, no.nets)

meanPermoduleTop   = matrix(0, no.modules, no.nets)
sdPermoduleTop     = matrix(0, no.modules, no.nets)
ks.pvaluesTop      = rep(1, no.modules)

heatmaps.list      = as.list(rep(NA, no.modules))
modsizes           = rep(NA, no.modules)

# ---- NEW: compute true MDC per module (yfold) and build heatmap tiles ----
lower_mean <- function(M) {
  if (is.null(M) || any(dim(M) == 0)) return(NA_real_)
  diag(M) <- 0
  m <- abs(M)
  mean(m[lower.tri(m)], na.rm = TRUE)
}

for (x in seq_len(no.modules)) {
  idx <- which(modulescolorB == modulenames[x])
  modsizes[x] <- length(idx)
  if (length(idx) < 2) {
    meanPermodule[x, ] <- c(NA_real_, NA_real_)
    next
  }
  # A = old, B = young (keep the order consistent with shortnames)
  A <- datExpr [idx, idx, drop = FALSE]
  B <- datExpr2[idx, idx, drop = FALSE]
  meanPermodule[x, 2] <- lower_mean(A)  # old
  meanPermodule[x, 1] <- lower_mean(B)  # young

  # heatmap tile: upper=young (B), lower=old (A)
  AU <- abs(B); diag(AU) <- 0
  AL <- abs(A); diag(AL) <- 0
  ord <- order(rowSums(AU, na.rm = TRUE), decreasing = TRUE)
  CB <- AU[ord, ord, drop = FALSE]
  CA <- AL[ord, ord, drop = FALSE]
  U <- CB; L <- CA
  U[lower.tri(U, diag = TRUE)] <- NA
  L[upper.tri(L, diag = TRUE)] <- NA
  M <- U; M[is.na(M)] <- L[is.na(M)]
  heatmaps.list[[x]] <- M
}

yfold <- meanPermodule[, 1] / meanPermodule[, 2]

USE_SAMPLE_PERMS <- PERMUTE_SAMPLES


#############################################################################
#-------------------- run permutations via validateMDC ---------------------

# perm_res_Genes50 <- validateMDC(no_perms = 50,
#                         permute_genes = PERMUTE_GENES,
#                         permute_samples = FALSE)

# perm_res_Genes100 <- validateMDC(no_perms = 100,
#                         permute_genes = PERMUTE_GENES,
#                         permute_samples = FALSE)

perm_res50 <- validateMDC(no_perms = 50,
                        permute_genes = PERMUTE_GENES,
                        permute_samples = TRUE)

perm_res100 <- validateMDC(no_perms = 100,
                        permute_genes = PERMUTE_GENES,
                        permute_samples = TRUE)

# GlobalyRandomMDC_genes50 <- perm_res_Genes50$random_genes
# GlobalyRandomMDC_genes100 <- perm_res_Genes100$random_genes

GlobalyRandomMDC_genes50 <- perm_res50$random_genes
GlobalyRandomMDC_genes100 <- perm_res100$random_genes
SampleRandomMDC_50 <- perm_res50$random_samples
SampleRandomMDC_100 <- perm_res100$random_samples

# Helper to write outputs for a single validation
write_validation_outputs <- function(run_tag, perm_res_obj) {
  RG <- perm_res_obj$random_genes
  RS <- perm_res_obj$random_samples

  # FDRs per family
  if (!is.null(RG)) {
    GR_FDR <- apply(cbind(yfold, RG), 1, MDC_FDR)
  } else {
    GR_FDR <- rep(NA_real_, length(yfold))
  }
  if (!is.null(RS)) {
    mdcFDR <- apply(cbind(yfold, RS), 1, MDC_FDR)
  } else {
    mdcFDR <- rep(NA_real_, length(yfold))
  }
  comFDR <- pmax(GR_FDR, mdcFDR, na.rm = TRUE)
  comFDR[!is.finite(comFDR)] <- NA_real_

  # Save randoms
  if (!is.null(RG)) {
    flog_genes  <- paste(outputDir1, fname, "_", run_tag, "_Randoms_genes", ".xls", sep='')
    xfinal <- cbind(module = modulenames, RG)
    colnames(xfinal) <- c("module", sprintf("MDC_random_genes_%d", seq_len(ncol(RG))))
    write.table(xfinal, flog_genes, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  }
  if (!is.null(RS)) {
    flog_samples  <- paste(outputDir1, fname, "_", run_tag, "_Randoms_samples", ".xls", sep='')
    xsamp <- cbind(module = modulenames, RS)
    colnames(xsamp) <- c("module", sprintf("MDC_random_samples_%d", seq_len(ncol(RS))))
    write.table(xsamp, flog_samples, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  }

  # Final table
  final <- cbind(modulenames, round(yfold, 4),
                 signif(comFDR, 4),
                 signif(mdcFDR, 4),
                 signif(GR_FDR, 4),
                 signif(if (!is.null(perm_res_obj$p_values)) perm_res_obj$p_values else NA_real_, 4),
                 signif(if (!is.null(perm_res_obj$p_values_samples)) perm_res_obj$p_values_samples else NA_real_, 4))
  colnames(final) <- c("module", "MDC", "FDR", "FDR_random_samples", "FDR_random_genes", "p_empirical_genes", "p_empirical_samples")

  flog2_run <- paste(outputDir1, fname, "_", run_tag, "_wFDR", ".xls", sep='')
  write.table(final, flog2_run, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

  # Per-run barplot
  rylab <- paste("MDC(", shortnames[1], ", ", shortnames[2], ")", sep="")
  od <- order(-yfold)
  dispDat <- cbind(yfold[od]); dispNames <- modulenames[od]
  rownames(dispDat) <- dispNames
  dispSign <- ifelse(comFDR[od] <= 0.10, "*", "")
  yymax <- closest_integer_bound(max(yfold))
  fimgRt_run <- paste(outputDir1, fname, "_", run_tag, ".png", sep='')
  barplot_by_oneGroup(datMatrix=dispDat, maskMatrix=cbind(dispSign), 
                      rowcol=col.names[1:length(dispNames)], vlines=c(2,1.5,1),
                      title="", ylab=rylab, imgfile=fimgRt_run, imgwid=2000, imghei=800, 
                      mcex=0.5, xlab_cex=0.65, xaxisLabel_rotate=45, legendcols=1, legendxDelta=0, legendCex=1, 
                      ipointsize=12, iunits="px", ires=300, icompression="lzw",
                      show_value=FALSE, showLegend=FALSE, show_xaxisLabel=TRUE, show_mask=TRUE,
                      xmargin=3, ymargin=5, zmargin=0.5, xlabel_orientation=0, ValueXdelta=0, 
                      ValueYdelta=2*yymax/100, valuecex=1, yMAX=yymax, yMIN=0)
}

# Four validations
# write_validation_outputs("genes_p50",  perm_res_Genes50)
write_validation_outputs("both_p50",   perm_res50)
# write_validation_outputs("genes_p100", perm_res_Genes100)
write_validation_outputs("both_p100",  perm_res100)


#################################################################################
# Per-validation barplots are generated above; removed single-run barplot.

#************************************************************************************
#
#  Plot Heatmaps
#
#************************************************************************************

orderedColor = modulescolorB

moduleHeatmaps = as.list(rep(NA, no.modules))

yhei=1000; ydelta = 0.2
if(no.modules<25) {
  rowfigs = as.integer(no.modules^0.5); subfigs = as.integer(no.modules/rowfigs+0.5)*rowfigs;
  if(subfigs < no.modules) {subfigs =subfigs + rowfigs}
  #subfigs = 12; rowfigs = 4
  print(paste(subfigs, rowfigs)) 
} else {
  subfigs = 24; rowfigs = 4
}

ywid = as.integer(yhei*(rowfigs^2)/subfigs);
tcexbase = 2.2 #1.5-unix #3.5, windows

#zrange=c(0, 1-max(max(dist1), max(dist1)) )
zrange=c(0,1)

module_order_by_size = FALSE

# ---------------- combined plot heatmap ----------------------
#
for (z in c(1:no.modules) ) {
  
  #image(heatmaps.list[[z]], xlim=c(0,1), ylim=c(0,1+ydelta), zlim=zrange, axes=F, col=heatmapColor)
  
  if(z%%subfigs==1) {
    end = subfigs+z-1; if(end>no.modules){end=no.modules}
    if(module_order_by_size) {
      imgHeatMap  =paste(outputDirTom, fname, "_iHMP0_", z, "-", end, ".png",sep='')  #png only
    } else{
      imgHeatMap  =paste(outputDirTom, fname, "_iHMP0_", modulenames[z], "-", modulenames[end], ".png",sep='')  #png only
    }
    
    openImgDev(imgHeatMap, iwidth =ywid, iheight = yhei)
    par(mfrow=c(subfigs/rowfigs,rowfigs), mar=c(0, 0, 0, 0),cex=1)
    #par(mfrow=c(2,2), mar=c(0, 0, 0, 0), cex=1)
    print(imgHeatMap)
  }
  
  par(mar=c(0.2, 0, 0.2, 0)+0.4)
  # image((heatmaps.list[[z]])^2, xlim=c(0,1), ylim=c(0,1+ydelta), zlim=zrange, axes=F, col=heatmapColor)
  # Compute z-limits from data
  Mz <- heatmaps.list[[z]]
  zr <- tryCatch(quantile(Mz[is.finite(Mz)], c(0.02, 0.98), na.rm=TRUE), error=function(e) c(0,1))
  if (!is.finite(diff(zr)) || diff(zr) <= 0) zr <- range(Mz[is.finite(Mz)], na.rm=TRUE)
  if (!all(is.finite(zr)) || diff(zr) <= 0) zr <- c(0, 1e-6)
  image(Mz, xlim=c(0,1), ylim=c(0,1+ydelta), zlim=zr, axes=F, col=heatmapColor)
  gradient.rect(0,1.05,1,1+ydelta,col=col.names[modulenames[z]],border=F)
  
  # font size
  colorchars = unlist(strsplit(modulenames[z], ""))
  if (length(colorchars)>12) {
    tcex = tcexbase#*(10/length(colorchars))^0.4
  } else{tcex=tcexbase}
  
  tcolor = SetTextContrastColor(col.names[z])
  text(0.5,1+ydelta*0.7,labels=modulenames[z],cex=tcex, col=tcolor)
  
  if(z%%subfigs==0) {
    par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    dev.off();
  }       
}

if(z%%subfigs>0) {
  
  scales = rev(seq(0,1, 0.01)); nscales = length(scales)
  scalesM= repmat(scales, 1, 100)
  textsel= (scales*100)%%25 ==0;
  scalesLabel = scales[textsel];
  #scalesLabel = as.integer(100*(zrange[1]+(zrange[2]-zrange[1])*scales[textsel])/100;
  ydel = 1.0
  
  par(mar=c(0, 0, 0, 0),cex=1)
  image(scalesM, xlim=c(-0.05,1.05), ylim=c(-ydel,1+1), zlim=zrange, axes=F, col=heatmapColor)

  # draw axis and ticks
  linecol = "brown"; ypos = -ydel/8; tickhei=ydel/4
  lines(x=c(0,1),y=c(ypos, ypos), col=linecol, lwd=1)
  for(iv in scalesLabel){
    lines(x=c(iv, iv), y=c(ypos,ypos-tickhei), col=linecol, lwd=1)
  }
  text(x=scalesLabel,y=rep(-ydel*0.7), labels=as.character(scalesLabel),cex=tcexbase*0.7, col="brown")
  
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
  dev.off();
  
}

#------------ plot individual heatmap ---------------------------------
#
yywid = 200
for (z in c(1:no.modules) ) {
  
  print(paste("Plot iHMP: ", modulenames[z]))
  
  imgHeatMap  =paste(outputDirTom, "/", fname, "_iHMP-", modulenames[z],".png",sep='')  #png only
  openImgDev(imgHeatMap, iwidth =yywid, iheight =yywid)
  par(mfrow=c(1,1), mar=c(0, 0, 0, 0),cex=1)
  
  # image((heatmaps.list[[z]])^2, xlim=c(0,1), ylim=c(0,1), zlim=zrange, axes=F, col=heatmapColor)
  Mz <- heatmaps.list[[z]]
  zr <- tryCatch(quantile(Mz[is.finite(Mz)], c(0.02, 0.98), na.rm=TRUE), error=function(e) c(0,1))
  if (!is.finite(diff(zr)) || diff(zr) <= 0) zr <- range(Mz[is.finite(Mz)], na.rm=TRUE)
  if (!all(is.finite(zr)) || diff(zr) <= 0) zr <- c(0, 1e-6)
  image(Mz, xlim=c(0,1), ylim=c(0,1), zlim=zr, axes=F, col=heatmapColor)

  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
  dev.off();
}

#quit(save = "no",status = 0)
rm(heatmaps.list)

