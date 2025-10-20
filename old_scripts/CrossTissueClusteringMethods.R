# This code uses functions from hustal/X-WGCNA Repository
# Repository URL: https://github.com/hustal/X-WGCNA
# Author(s): Husain Ahammad Talukdar
# License: GNU General Public License v3.0 (GPLv3)
# The original code has been modified for this project.
# The modified code follows the same GPLv3 license.


library(WGCNA)
options(stringsAsFactors = FALSE)
suppressMessages(WGCNA::allowWGCNAThreads())  # multi-thread WGCNA when available
.extract_donor <- function(x) sub("^([^-]+-[^-]+).*", "\\1", x)

# Average duplicate samples per donor (rows) within a tissue
.aggregate_by_donor <- function(mat) {
  d <- .extract_donor(rownames(mat))
  split_idx <- split(seq_len(nrow(mat)), d)
  res <- do.call(rbind, lapply(split_idx, function(ix) {
    colMeans(mat[ix, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(res) <- names(split_idx)
  storage.mode(res) <- "double"
  res
}

# -------------------Load gene expression data  --------------
LoadExprData<-function(tissue_name, tissue_file_name, 
                       # MV_sd_thr
                       sd_quantile = 0.90,
                       max_genes_per_tissue = 4000
                       ) {
  
  if (grepl("\\.RData$", tissue_file_name)){
    datExpr <- get(load(tissue_file_name))
  } else if (grepl("\\.csv$", tissue_file_name)) {
    datExpr <- read.csv(tissue_file_name, check.names = FALSE)
    rownames(datExpr) <- datExpr[,1]
    datExpr <- as.matrix(datExpr[,-1,drop=FALSE])
  } else stop("Unsupported input file !!!")
  
  # present <-colMeans(!is.na(datExpr) & datExpr > min_expr) >= min_prop
  # datExpr <- datExpr[, present, drop = FALSE]
  
  sds <- apply(datExpr, 2, sd, na.rm = TRUE)
  thr <- stats::quantile(sds, probs = sd_quantile, na.rm = TRUE)
  keep_sd <- sds >= thr
  datExpr <- datExpr[, keep_sd, drop = FALSE]
  
  
  if (ncol(datExpr) > max_genes_per_tissue) {
    o <- order(sds[colnames(datExpr)], decreasing = TRUE)
    datExpr <- datExpr[, o[seq_len(max_genes_per_tissue)], drop = FALSE]
  }
  
  colnames(datExpr) <- paste0(tissue_name, "_", colnames(datExpr))
  storage.mode(datExpr) <- "double"
  datExpr
}
  
  #rows indicate samples and columns indicate gene symbols
  # if(substr(tissue_file_name, start = nchar(tissue_file_name)-3, stop = nchar(tissue_file_name)) == ".RData") {
    
    # datExpr <- get(load(tissue_file_name))
  # } else if(substr(tissue_file_name, start = nchar(tissue_file_name)-3, stop = nchar(tissue_file_name)) == ".csv") {
    
    # datExpr <-read.csv(tissue_file_name)
    # datExpr <- as.matrix(datExpr)
    # rownames(datExpr) <- datExpr[ ,1]
    # datExpr <- datExpr[ ,-1]
    # datExpr <- apply(datExpr, c(1,2), as.numeric)
#   } else stop("Unsupported input file !!!")
  
  # colnames(datExpr) <- paste(tissue_name,"_",colnames(datExpr),sep="")
 
#   return(datExpr)
  
# }


AdjacencyFromExpr <- function(
    tissue_names = NULL,
    tissue_expr_file_names = NULL,
    sd_quantile = 0.90,
    max_genes_per_tissue = 4000,
    min_expr = 0,
    min_prop = 0.10,
    TS_power = 6, CT_power = 3,
    cor_method = "pearson",
    zero_fill_if_no_overlap = TRUE   # <- new
){
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  total_tissues <- length(tissue_names)
  
  vector_expr <- vector("list", total_tissues)
  donors_mat  <- vector("list", total_tissues)  # donor-aggregated per tissue (for CT)
  tissue_index_adj <- integer(total_tissues + 1)
  rc_names <- character()
  
  # Load + filter each tissue; also precompute donor-aggregated matrices
  for (i in seq_len(total_tissues)) {
    mat <- LoadExprData(
      tissue_names[i], tissue_expr_file_names[i],
      sd_quantile = sd_quantile,
      max_genes_per_tissue = max_genes_per_tissue
    )
    
    if (ncol(mat) == 0) {
      stop(sprintf("[%s] 0 genes after filtering; lower sd_quantile or max_genes_per_tissue, or lower min_prop.",
                   tissue_names[i]))
    }
    if (nrow(mat) < 3) {
      warning(sprintf("[%s] fewer than 3 samples; correlations may be unstable.", tissue_names[i]))
    }
    
    vector_expr[[i]] <- mat
    donors_mat[[i]]  <- .aggregate_by_donor(mat)
    
    tissue_index_adj[i + 1] <- tissue_index_adj[i] + ncol(mat)
    rc_names <- c(rc_names, colnames(mat))
    
    message(sprintf("[%s] kept %d genes across %d samples (%d donors)",
                    tissue_names[i], ncol(mat), nrow(mat), nrow(donors_mat[[i]])))
  }
  
  # Allocate adjacency
  adj_mat <- matrix(0, nrow = tissue_index_adj[total_tissues + 1],
                    ncol = tissue_index_adj[total_tissues + 1],
                    dimnames = list(rc_names, rc_names))
  
  # Within-tissue blocks (TS)
  for (i in 1:total_tissues) {
    idx_i <- (tissue_index_adj[i] + 1):tissue_index_adj[i + 1]
    adj_mat[idx_i, idx_i] <- abs(cor(
      vector_expr[[i]], use = "pairwise.complete.obs", method = cor_method
    ))^TS_power
  }
  
  # Cross-tissue blocks (CT) via donor intersection
  if (total_tissues >= 2) {
    for (i in 1:(total_tissues - 1)) {
      idx_i <- (tissue_index_adj[i] + 1):tissue_index_adj[i + 1]
      for (j in (i + 1):total_tissues) {
        idx_j <- (tissue_index_adj[j] + 1):tissue_index_adj[j + 1]
        
        di <- rownames(donors_mat[[i]])
        dj <- rownames(donors_mat[[j]])
        common <- intersect(di, dj)
        
        if (length(common) < 3) {
          msg <- sprintf("No/too-few common donors between %s and %s (|common|=%d).",
                         tissue_names[i], tissue_names[j], length(common))
          if (zero_fill_if_no_overlap) {
            message(msg, " Filling CT block with zeros.")
            next
          } else {
            stop(msg, " Consider lowering filters or using consensus WGCNA.")
          }
        }
        
        Mi <- donors_mat[[i]][common, , drop = FALSE]
        Mj <- donors_mat[[j]][common, , drop = FALSE]
        
        if (nrow(Mi) == 0 || nrow(Mj) == 0 || ncol(Mi) == 0 || ncol(Mj) == 0) {
          message(sprintf("Empty block for %s vs %s; filling zeros.", tissue_names[i], tissue_names[j]))
          next
        }
        
        C <- abs(cor(Mi, Mj, use = "pairwise.complete.obs", method = cor_method))^CT_power
        adj_mat[idx_i, idx_j] <- C
        adj_mat[idx_j, idx_i] <- t(C)
      }
    }
  }
  
  saveRDS(adj_mat, "adj_mat.rds")
  adj_mat
}




# -------------------Adjacency matrix from all tissue expression data ----------------------------------------- 
AdjacencyFromExpr_old <- function(tissue_names = NULL, tissue_expr_file_names = NULL, MV_sd_thr = 0.5, TS_power = 6, CT_power = 3, cor_method = "pearson") {
  print("calculating adj matrix from expression")
  total_tissues <- length(tissue_names)
  
  vector_expr<-vector(mode="list", length=total_tissues)
  
  tissue_index_adj<-vector(mode="integer", length=(total_tissues+1))
  
  rc_names<-c()
  
  # Initialize the first index
  for(i in 1:total_tissues) {
    # In each iteration, load the expression data for the tissue and store it in the list
    TS_expr_data<-LoadExprData(tissue_names[i], tissue_expr_file_names[i], MV_sd_thr)
    vector_expr[[i]]<-TS_expr_data
    tissue_index_adj[i+1]<- tissue_index_adj[i]+ncol(vector_expr[[i]])
    rc_names<-c(rc_names,colnames(vector_expr[[i]]))
  }
  
  adj_mat<-matrix(0,nrow=tissue_index_adj[total_tissues+1], ncol=tissue_index_adj[total_tissues+1])
  rownames(adj_mat)<-rc_names
  colnames(adj_mat)<-rc_names
  
  for(i in 1:(total_tissues-1)) {
    adj_mat[(tissue_index_adj[i]+1):tissue_index_adj[i+1], (tissue_index_adj[i]+1):tissue_index_adj[i+1]] <- abs(cor(vector_expr[[i]], method = cor_method))^TS_power
    
    for(j in (i+1):total_tissues) {
      common_Samples <- intersect(rownames(vector_expr[[i]]),rownames(vector_expr[[j]]))
      adj_mat[(tissue_index_adj[i]+1):tissue_index_adj[i+1],(tissue_index_adj[j]+1):tissue_index_adj[j+1]] <- abs(cor(vector_expr[[i]][common_Samples,],vector_expr[[j]][common_Samples,], method = cor_method))^CT_power
      adj_mat[(tissue_index_adj[j]+1):tissue_index_adj[j+1],(tissue_index_adj[i]+1):tissue_index_adj[i+1]] <- t(adj_mat[(tissue_index_adj[i]+1):tissue_index_adj[i+1],(tissue_index_adj[j]+1):tissue_index_adj[j+1]])
    }
  }
  
  adj_mat[(tissue_index_adj[total_tissues]+1):tissue_index_adj[total_tissues+1],(tissue_index_adj[total_tissues]+1):tissue_index_adj[total_tissues+1]] <- abs(cor(vector_expr[[total_tissues]], method = cor_method))^TS_power
  saveRDS(adj_mat, 'adj_mat.rds')
  return(adj_mat)
}


Cross_Tissue_TOM <- function(adj_mat, TOMType = "unsigned") {
  rn <- rownames(adj_mat)
  cn <- colnames(adj_mat)
  TOM <- WGCNA::TOMsimilarity(adj_mat, TOMType = TOMType)
  dimnames(TOM) <- list(rn, cn)
  TOM
  # If memory is still tight, consider lowering max_genes_per_tissue further.
  
}


# ------------------TOM from Adjacency matrix ----------------------- 

Cross_Tissue_TOM_old <- function(adj_mat, block_size = 2000L) {
  message("Computing TOM from adjacency matrix... (old)")
  n <- nrow(adj_mat)
  rn <- rownames(adj_mat)
  cn <- colnames(adj_mat)
  k <- rowSums(adj_mat)
  message("Number of genes: ", n, "; Number of samples: ", ncol(adj_mat))
  TOM_mat <- matrix(0, n, n)
  for (i in seq(1, n, by = block_size)){
    message("Processing block ", i, " to ", min(i + block_size - 1, n))
    i2 <- min(i + block_size - 1, n)
    Ai <- adj_mat[i:i2, , drop = FALSE]
    ki <- k[i:i2]

    for (j in seq(i, n, by = block_size)) {
      j2 <- min(j + block_size - 1, n)
      Aj <- adj_mat[j:j2, , drop = FALSE]
      kj <- k[j:j2]
      
      num <- Ai %*% t(Aj)
      Aij <- adj_mat[i:i2, j:j2, drop = FALSE]
      denom <- outer(ki, kj, pmin) + 1 - Aij

      block <- (Aij + num) / denom
      TOM_mat[i:i2, j:j2] <- block
      if (j > i) TOM_mat[j:j2, i:i2] <- t(block)
      rm(Aj, num, Aij, denom, block)
      gc(FALSE)
    }
    message("Processed block ", i, " to ", min(i + block_size - 1, n))
    rm(Ai)
    gc(FALSE)
  }
  dimnames(TOM_mat) <- list(rn, cn)
  message("TOM computation done.")
  TOM_mat
  # matrix_product <- tcrossprod(adj_mat)
  # # matrix_product <- adj_mat %*% t(adj_mat)
  # print("Matrix Multiplication done...")
  # rsums<-rowSums(adj_mat)
  # D_val<-matrix(rsums,nrow=nrow(adj_mat),ncol=ncol(adj_mat))
  
  # total_col <- ncol(D_val)
  # if(total_col <= 5000) col_limit <- total_col else col_limit <- 5000
  # start_col <- 1
  # end_col <- col_limit
  
  # tmp_D <- matrix(0, nrow = total_col, ncol = total_col)
  
  # repeat {
  #   tmp_D[ ,start_col:end_col] <- pmin(D_val[ ,start_col:end_col], t(D_val)[ ,start_col:end_col])
  #   if(end_col==total_col) break
  #   start_col <- end_col + 1
  #   end_col <- end_col + col_limit
  #   if(end_col > total_col) end_col <- total_col
  # }  
  
  # rm(D_val)
  # rm(rsums)
  # gc()
  # print("Matrix pmin done...")
  
  # TOM_mat <- matrix(NA, nrow=nrow(adj_mat), ncol=ncol(adj_mat))
  # rownames(TOM_mat) <- rownames(adj_mat)
  # colnames(TOM_mat) <- colnames(adj_mat)
  
  # TOM_mat <- (adj_mat + matrix_product) / (tmp_D + 1 - adj_mat)
  
  # print("Tom done...")
  
  # return(TOM_mat) 
  
}

# -------- Optional: blockwise TOM (memory saver) ----------
.compute_TOM_blockwise <- function(adj_mat, TOMType = "unsigned", block_size = 2000L) {
  # Only unsigned implemented (same as WGCNA unsigned)
  n <- nrow(adj_mat)
  rn <- rownames(adj_mat); cn <- colnames(adj_mat)
  k  <- rowSums(adj_mat)
  TOM <- matrix(0, n, n)
  for (i in seq(1, n, by = block_size)) {
    i2 <- min(i + block_size - 1, n)
    Ai <- adj_mat[i:i2, , drop = FALSE]
    for (j in seq(i, n, by = block_size)) {
      j2 <- min(j + block_size - 1, n)
      Aj <- adj_mat[j:j2, , drop = FALSE]
      num <- Ai %*% t(Aj)
      Aij <- adj_mat[i:i2, j:j2, drop = FALSE]
      denom <- outer(k[i:i2], k[j:j2], pmin) + 1 - Aij
      block <- (Aij + num) / denom
      TOM[i:i2, j:j2] <- block
      if (j > i) TOM[j:j2, i:i2] <- t(block)
    }
    rm(Ai); gc()
  }
  dimnames(TOM) <- list(rn, cn)
  TOM
}

Cross_Tissue_TOM <- function(adj_mat, TOMType = "unsigned", blockwise = FALSE, block_size = 2000L) {
  if (!blockwise) {
    rn <- rownames(adj_mat); cn <- colnames(adj_mat)
    TOM <- WGCNA::TOMsimilarity(adj_mat, TOMType = TOMType)
    dimnames(TOM) <- list(rn, cn)
    return(TOM)
  } else {
    message("Using blockwise TOM computation...")
    return(.compute_TOM_blockwise(adj_mat, TOMType = TOMType, block_size = block_size))
  }
}

# -------- Refactored heatmap (no hard-coded paths) --------
network_heatmap <- function(TOM_mat, dynamicColors, restGenes, out_file = "TOMplot.png") {
  if (!requireNamespace("gplots", quietly = TRUE)) {
    warning("gplots not installed; skipping heatmap.")
    return(invisible(NULL))
  }
  library(gplots)
  keep_idx <- which(restGenes)
  if (length(keep_idx) < 4) {
    warning("Too few genes for heatmap; skipped.")
    return(invisible(NULL))
  }
  subTOM <- TOM_mat[keep_idx, keep_idx, drop = FALSE]
  dist_sub <- 1 - subTOM
  h_sub <- hclust(as.dist(dist_sub), method = "average")
  plotTOM <- dist_sub^4
  diag(plotTOM) <- NA
  myheatcol <- colorpanel(250, 'red', "orange", 'lemonchiffon')
  png(out_file, width = 900, height = 900, res = 140)
  WGCNA::TOMplot(plotTOM, h_sub, as.character(dynamicColors[keep_idx]),
                 main = "Network heatmap (module genes)", col = myheatcol)
  dev.off()
  invisible(out_file)
}

#----------------- Plot the resulting clustering tree (dendrogram)-----------------
network_heatmap_old <- function(restGenes, dynamicColors) { 
  
  library(gplots)
  myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
  adj_mat_filter <- AdjacencyFromExpr(c("Adipose", "Muscle", "Brain"), c("/Users/edeneldar/CoExpression_ReProduction/old_scripts/outputs/Adipose - Subcutaneous/young_matrix.csv","/Users/edeneldar/CoExpression_ReProduction/old_scripts/outputs/Muscle - Skeletal/young_matrix.csv", "/Users/edeneldar/CoExpression_ReProduction/old_scripts/outputs/Brain - Cortex/young_matrix.csv"), 0.5)
  
  TOM_mat_filter <- Cross_Tissue_TOM(adj_mat_filter[restGenes, restGenes])
  rm(adj_mat_filter)
  gc()
  
  dist_TOM_mat_filter <- 1-TOM_mat_filter
  rm(TOM_mat_filter)
  gc()
  
  h_TOM_filter <- hclust(as.dist(dist_TOM_mat_filter), method="average")
  plotTOM = dist_TOM_mat_filter^4; 
  # Set diagonal to NA for a nicer plot 
  diag(plotTOM) = NA; 
  
  # Call the plot function 
  sizeGrWindow(9,9)
  png(filename="TOMplot.png")
  TOMplot(plotTOM, h_TOM_filter, as.character(dynamicColors[restGenes]),main = "Network heatmap plot, module genes", col=myheatcol)
  
  #plot()
  TOMplot
  dev.off()
  
}


# ------------------Cross-Tissue Clusters Table from TOM matrix --------------
Clusters_Table <- function(TOM_mat, minClusterSize = 30, plot_heatmap = FALSE) {
  library(dynamicTreeCut)
  
  # Robust names
  gene_names <- colnames(TOM_mat)
  if (is.null(gene_names)) gene_names <- rownames(TOM_mat)
  
  dist_TOM_mat <- 1 - TOM_mat
  # rm(TOM_mat); gc()
  h_TOM <- hclust(as.dist(dist_TOM_mat), method = "average")
  
  sizeGrWindow(12,9)
  plot(h_TOM, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04)
  
  dynamicMods   <- cutreeDynamic(h_TOM, method = "tree", deepSplit = TRUE, minClusterSize = minClusterSize)
  dynamicColors <- labels2colors(dynamicMods)
  restGenes <- (dynamicColors != "grey")
  message("Plot heatmap", sum(restGenes))

  if (plot_heatmap && sum(restGenes) <= 5000) {
    network_heatmap(restGenes, dynamicColors)
  }
  
  module_labels <- sort(setdiff(unique(as.integer(dynamicMods)), 0L))
  if (length(module_labels) == 0L) {
    warning("No non-grey modules detected.")
    empty <- matrix(nrow = 0, ncol = 3)
    colnames(empty) <- c("Cluster ID", "Tissue", "Gene Symbol")
    return(empty)
  }
  
  # Safe split of "Tissue_Gene"
  split_pos  <- regexpr("_", gene_names, fixed = TRUE)
  tissue_vec <- ifelse(split_pos > 0, substr(gene_names, 1, split_pos - 1), "Unknown")
  gene_vec   <- ifelse(split_pos > 0, substr(gene_names, split_pos + 1, nchar(gene_names)), gene_names)
  
  # Build clusters table (preserve names!)
  df_list <- vector("list", length(module_labels))
  for (k in seq_along(module_labels)) {
    lab <- module_labels[k]
    idx <- which(dynamicMods == lab)
    if (length(idx) == 0L) next
    
    df_list[[k]] <- data.frame(
      "Cluster ID"  = rep.int(k, length(idx)),
      "Tissue"      = tissue_vec[idx],
      "Gene Symbol" = gene_vec[idx],
      stringsAsFactors = FALSE,
      check.names = FALSE   # <-- keep "Cluster ID" with a space
    )
  }
  
  clusters_table <- do.call(rbind, Filter(Negate(is.null), df_list))
  clusters_table <- as.matrix(clusters_table)
  return(clusters_table)
}


Clusters_Table_old <- function(TOM_mat, minClusterSize = 30, plot_heatmap = FALSE) {
  
  library('dynamicTreeCut')
  dist_TOM_mat <- 1 - TOM_mat
  rm(TOM_mat); gc()
  
  
  h_TOM <- hclust(as.dist(dist_TOM_mat), method="average")
  
  # Plot the resulting clustering tree (dendrogram)
  sizeGrWindow(12,9)
  plot(h_TOM, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04);
  
  dynamicMods = cutreeDynamic(h_TOM, method = "tree", deepSplit = TRUE, minClusterSize = minClusterSize);
  
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  
  # Plot the dendrogram and colors underneath
  sizeGrWindow(8,6)
  png(filename="plotDendroAndColors.png")
  
  plotDendroAndColors(h_TOM, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  
  #plot()
  dev.off()
  
  
  restGenes= (dynamicColors != "grey")
  if (plot_heatmap && sum(restGenes) <= 5000){
    network_heatmap(restGenes, dynamicColors)
  }
  
  all_gene_names <- h_TOM$labels
  total_clusters <- length(table(dynamicMods)) - 1
  
  vector_clusters <- vector(mode="list", length=(total_clusters))
  for(i in 1:total_clusters) {
    
    index <- which(dynamicMods == i)
    vector_clusters[[i]] <- all_gene_names[index]
    names(vector_clusters[[i]]) <- index
    
  }
  
  clusters_table <- do.call(rbind, lapply(seq_along(vector_clusters), function(i) {data.frame(Cluster_ID=i, Gene_Symbol=vector_clusters[[i]])}))
  
  clusters_table <- as.matrix(clusters_table)
  clusters_table <- cbind(clusters_table, '')
  colnames(clusters_table) <- c('Cluster ID', 'Tissue', 'Gene Symbol')
  
  clusters_table[ ,c(2,3)] <- cbind(unlist(lapply(strsplit(clusters_table[ ,2], split = '_'), function(x) {return(x[1])})), unlist(lapply(strsplit(clusters_table[ ,2], split = '_'), function(x) {return(x[2])})))
  
  return(clusters_table)
}



Clusters_Table_old <- function(TOM_mat, minClusterSize = 30) {
  
  library('dynamicTreeCut')
  dist_TOM_mat <- 1 - TOM_mat
  rm(TOM_mat)
  gc()

  h_TOM <- hclust(as.dist(dist_TOM_mat), method="average")

  # Plot the resulting clustering tree (dendrogram)
  sizeGrWindow(12,9)
  plot(h_TOM, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04);
  
  dynamicMods = cutreeDynamic(h_TOM, method = "tree", deepSplit = TRUE, minClusterSize = minClusterSize);
  
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  
  # Plot the dendrogram and colors underneath
  sizeGrWindow(8,6)
  png(filename="plotDendroAndColors.png")
  
  plotDendroAndColors(h_TOM, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  
  #plot()
  dev.off()
  

  restGenes= (dynamicColors != "grey")
  network_heatmap(restGenes, dynamicColors)
  
  all_gene_names <- h_TOM$labels
  total_clusters <- length(table(dynamicMods)) - 1
 
  vector_clusters <- vector(mode="list", length=(total_clusters))
  for(i in 1:total_clusters) {
    
    index <- which(dynamicMods == i)
    vector_clusters[[i]] <- all_gene_names[index]
    names(vector_clusters[[i]]) <- index
    
  }
  
  clusters_table <- do.call(rbind, lapply(seq_along(vector_clusters), function(i) {data.frame(Cluster_ID=i, Gene_Symbol=vector_clusters[[i]])}))
  
  clusters_table <- as.matrix(clusters_table)
  clusters_table <- cbind(clusters_table, '')
  colnames(clusters_table) <- c('Cluster ID', 'Tissue', 'Gene Symbol')
  
  clusters_table[ ,c(2,3)] <- cbind(unlist(lapply(strsplit(clusters_table[ ,2], split = '_'), function(x) {return(x[1])})), unlist(lapply(strsplit(clusters_table[ ,2], split = '_'), function(x) {return(x[2])})))
  
  return(clusters_table)
}

  

 # -----------------Clusters Deatils Information from clusters table ---------------
Clusters_Details <- function(clusters_table, cluster_type_thr = 0.95) {
  if (is.null(clusters_table) || nrow(clusters_table) == 0) {
    warning("clusters_table is empty; returning empty details.")
    # Build an empty details table with expected columns
    return(matrix(nrow = 0, ncol = 0))
  }
  
  # Normalize column names: turn dots into spaces
  cn <- colnames(clusters_table)
  colnames(clusters_table) <- sub("\\.", " ", cn)
  
  required <- c("Cluster ID", "Tissue", "Gene Symbol")
  missing  <- setdiff(required, colnames(clusters_table))
  if (length(missing) > 0) {
    stop("clusters_table missing columns: ", paste(missing, collapse = ", "))
  }
  
  tissue_names   <- names(table(clusters_table[, "Tissue"]))
  total_clusters <- max(as.numeric(clusters_table[, "Cluster ID"]))
  
  clusters_table_details <- matrix(0, nrow = total_clusters, ncol = (5 + length(tissue_names)))
  colnames(clusters_table_details) <- c('Cluster ID', 'Cluster Size', 'Cluster Type',
                                        'Cluster Tissues', tissue_names, 'Dominant Tissue')
  
  for (i in 1:total_clusters) {
    temp_cluster <- clusters_table[as.numeric(clusters_table[, "Cluster ID"]) == i, , drop = FALSE]
    
    clusters_table_details[i, "Cluster ID"]   <- i
    clusters_table_details[i, "Cluster Size"] <- nrow(temp_cluster)
    
    if (nrow(temp_cluster) == 0) {
      clusters_table_details[i, "Cluster Type"] <- 'NA'
      next
    }
    
    prop_tbl <- round(table(temp_cluster[, "Tissue"]) / nrow(temp_cluster), 2)
    clusters_table_details[i, "Cluster Type"] <- ifelse(max(prop_tbl) >= cluster_type_thr, 'TS', 'CT')
    
    clusters_table_details[i, "Cluster Tissues"] <- paste(names(table(temp_cluster[, "Tissue"])), collapse = ',')
    clusters_table_details[i, names(table(temp_cluster[, "Tissue"]))] <- table(temp_cluster[, "Tissue"])
    clusters_table_details[i, "Dominant Tissue"] <- names(which.max(table(temp_cluster[, "Tissue"])))
  }
  
  clusters_table_details
}



Clusters_Details_old <- function(clusters_table, cluster_type_thr = 0.95) {
  
  tissue_names <- names(table(clusters_table[ ,"Tissue"]))
  total_clusters <- max(as.numeric(clusters_table[ ,"Cluster ID"]))
  
  clusters_table_details <- matrix(0, nrow = total_clusters, ncol = (5+length(tissue_names)))
  colnames(clusters_table_details) <- c('Cluster ID', 'Cluster Size', 'Cluster Type', 'Cluster Tissues', tissue_names, 'Dominant Tissue')
  
  for(i in 1:total_clusters) {
    
    temp_cluster <- clusters_table[as.numeric(clusters_table[ ,"Cluster ID"]) == i, ]
    
    clusters_table_details[i, "Cluster ID"] <- i
    
    clusters_table_details[i, "Cluster Size"] <- nrow(temp_cluster)
    
    if(max(round(table(temp_cluster[ ,"Tissue"])/nrow(temp_cluster), 2)) >= cluster_type_thr) clusters_table_details[i, "Cluster Type"] <- 'TS' else clusters_table_details[i, "Cluster Type"] <- 'CT'
    
    clusters_table_details[i, "Cluster Tissues"] <- paste(names(table(temp_cluster[ ,"Tissue"])), collapse = ',')
    
    clusters_table_details[i, names(table(temp_cluster[ ,"Tissue"]))] <- table(temp_cluster[ ,"Tissue"])
    
    clusters_table_details[i, "Dominant Tissue"] <- names(which.max(round(table(temp_cluster[ ,"Tissue"]))))
    
  }
  
  return(clusters_table_details)
}

XWGCNA_Clusters_autoBeta <- function(
    tissue_names = NULL,
    tissue_expr_file_names = NULL,
    # filtering (forwarded)
    sd_quantile = 0.90,
    max_genes_per_tissue = 4000,
    min_expr = 0,
    min_prop = 0.10,
    TS_power = 6,
    CT_power = 3,
    cor_method = "pearson",
    TOMType = "unsigned",
    minClusterSize = 30,
    cluster_type_thr = 0.95,
    out_prefix = "xwgcna",
    save_intermediates = TRUE,
    plot_heatmap = FALSE,
    auto_beta = TRUE,
    targetR2 = 0.80,
    beta_grid = c(1:10, seq(12, 20, 2)),
    plot_beta_curves = TRUE,
    blockwise_TOM = FALSE,
    TOM_block_size = 2000L
){
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  
  beta_info <- NULL
  if (auto_beta) {
    message("Auto-picking TS/CT betas ...")
    beta_info <- auto_pick_powers(
      tissue_names, tissue_expr_file_names,
      sd_quantile = sd_quantile,
      max_genes_per_tissue = max_genes_per_tissue,
      min_expr = min_expr, min_prop = min_prop,
      cor_method = cor_method,
      beta_grid = beta_grid,
      targetR2 = targetR2
    )
    TS_power <- beta_info$TS_power
    CT_power <- beta_info$CT_power
    message(sprintf("Auto β selected: TS=%d (per tissue: %s); CT=%d",
                    TS_power, paste(beta_info$TS_per_tissue, collapse = ","), CT_power))
    if (plot_beta_curves && !is.null(beta_info$TS_fit_curves)) {
      plot_TS_beta_curves(
        beta_info$TS_fit_curves, TS_power, beta_info$TS_per_tissue,
        targetR2 = targetR2,
        out_prefix = paste0(out_prefix, "_TS_beta_curves"))
    }
    if (plot_beta_curves && !is.null(beta_info$CT_fit_curves)) {
      plot_CT_beta_curves(
        beta_info$CT_fit_curves, CT_power,
        targetR2 = targetR2,
        out_prefix = paste0(out_prefix, "_CT_beta_curves"))
    }
  }

  message(sprintf(
    "Adjacency: %d tissues, up to %d genes/tissue (sd_quantile=%.2f, min_prop=%.2f).",
    length(tissue_names), max_genes_per_tissue, sd_quantile, min_prop
  ))


  # if (auto_beta) {
  #   message("Auto-picking TS/CT betas ...")
  #   pow <- auto_pick_powers(
  #     tissue_names, tissue_expr_file_names,
  #     sd_quantile = sd_quantile,
  #     max_genes_per_tissue = max_genes_per_tissue,
  #     min_expr = min_expr, min_prop = min_prop,
  #     cor_method = cor_method,
  #     beta_grid = beta_grid,
  #     targetR2 = targetR2
  #   )
  #   TS_power <- pow$TS_power
  #   CT_power <- pow$CT_power
  #   message(sprintf("Auto β selected: TS=%d (per tissue: %s); CT=%d",
  #                   TS_power, paste(pow$TS_per_tissue, collapse = ","), CT_power))
  # }
  
  # message(sprintf(
  #   "Adjacency: %d tissues, up to %d genes/tissue (sd_quantile=%.2f, min_prop=%.2f).",
  #   length(tissue_names), max_genes_per_tissue, sd_quantile, min_prop
  # ))
  
  adj_mat <- AdjacencyFromExpr(
    tissue_names = tissue_names,
    tissue_expr_file_names = tissue_expr_file_names,
    sd_quantile = sd_quantile,
    max_genes_per_tissue = max_genes_per_tissue,
    min_expr = min_expr,
    min_prop = min_prop,
    TS_power = TS_power,
    CT_power = CT_power,
    cor_method = cor_method
  )
  
  if (save_intermediates) {
    saveRDS(adj_mat, paste0(out_prefix, "_adjacency.rds"))
  }
  gc()
  
  message("Computing TOM (this may take a while)...")
  # TOM_mat <- Cross_Tissue_TOM_old(adj_mat, TOMType = TOMType, blockwise = blockwise_TOM, block_size = TOM_block_size)
  TOM_mat <- Cross_Tissue_TOM_old(adj_mat)
  if (save_intermediates) saveRDS(TOM_mat, paste0(out_prefix, "_TOM.rds"))

  rm(adj_mat); gc()
  
  message("Clustering on TOM…")
  clusters_table <- Clusters_Table(TOM_mat, minClusterSize = minClusterSize, plot_heatmap = plot_heatmap)
  if (!plot_heatmap) rm(TOM_mat); gc()

    
  clusters_details <- Clusters_Details(clusters_table, cluster_type_thr = cluster_type_thr)
  
  write.table(clusters_table,
              file = paste0(out_prefix, "_Clusters_table.txt"),
              quote = FALSE, sep = '\t', row.names = FALSE)
  write.table(clusters_details,
              file = paste0(out_prefix, "_Clusters_details.txt"),
              quote = FALSE, sep = '\t', row.names = FALSE)
  
  # invisible(list(clusters_table = clusters_table, clusters_details = clusters_details, TS_power = TS_power, CT_power = CT_power, beta_info = beta_info))
  # Return the clusters table, details, and beta info (if auto_beta was used)
  # This allows further analysis or reporting without needing to recompute.
  # If auto_beta was used, beta_info will contain the selected powers and fit curves.
  # If not, beta_info will be NULL.
  clusters_table  <- as.data.frame(clusters_table,  stringsAsFactors = FALSE)
  clusters_details <- as.data.frame(clusters_details, stringsAsFactors = FALSE)

  if (auto_beta) {
    return(list(clusters_table = clusters_table, clusters_details = clusters_details,
                TS_power = TS_power, CT_power = CT_power, beta_info = beta_info))
  } else {
    return(list(clusters_table = clusters_table, clusters_details = clusters_details,
                TS_power = TS_power, CT_power = CT_power))
  }
}


R2_connectivity <- function(adj_mat) {
  k <- rowSums(adj_mat)
  r2 <- checkScaleFree(k)$Rsquared.SFT
  r2_conn <- list(r2 = r2,
                  mean_conn = mean(k),
                  median_conn = median(k),
                  max_conn = max(k),
                  min_conn = min(k))

  tissues <- sub("_.*$", "", colnames(adj_mat))
  counts  <- as.numeric(table(tissues))
  cumulative <- c(0, cumsum(counts))  # boundaries
  TS_conn_mat <- matrix(NA, nrow = length(counts), ncol = 4,
                        dimnames = list(names(table(tissues)), c('mean','median','max','min')))
  CT_conn_mat <- matrix(NA, nrow = (length(counts)*(length(counts)-1))/2, ncol = 4,
                        dimnames = list(NULL, c('mean','median','max','min')))
  ct_idx <- 1
  for (i in seq_along(counts)) {
    rows <- (cumulative[i] + 1):cumulative[i + 1]
    subA <- adj_mat[rows, rows, drop = FALSE]
    kk <- rowSums(subA)
    TS_conn_mat[i, ] <- c(mean(kk), median(kk), max(kk), min(kk))
    if (i < length(counts)) {
      for (j in (i + 1):length(counts)) {
        cols <- (cumulative[j] + 1):cumulative[j + 1]
        block <- adj_mat[rows, cols, drop = FALSE]
        kblock <- rowSums(block)
        CT_conn_mat[ct_idx, ] <- c(mean(kblock), median(kblock), max(kblock), min(kblock))
        ct_idx <- ct_idx + 1
      }
    }
  }
  list(
    r2 = r2_conn$r2,
    mean_conn = r2_conn$mean_conn,
    TS_conn_mat = TS_conn_mat,
    CT_conn_mat = CT_conn_mat
  )
}


XWGCNA_Clusters <- function(
    tissue_names = NULL,
    tissue_expr_file_names = NULL,
    sd_quantile = 0.90,
    max_genes_per_tissue = 4000,
    min_expr = 0,
    min_prop = 0.10,
    # network building knobs
    TS_power = 6,
    CT_power = 3,
    cor_method = "pearson",
    TOMType = "unsigned",
    # clustering / reporting
    minClusterSize = 30,
    cluster_type_thr = 0.95,
    # IO / UX
    out_prefix = "xwgcna",
    save_intermediates = TRUE,
    plot_heatmap = FALSE
){
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  
  message(sprintf(
    "Adjacency: %d tissues, up to %d genes/tissue (sd_quantile=%.2f, min_prop=%.2f).",
    length(tissue_names), max_genes_per_tissue, sd_quantile, min_prop
  ))
  
  adj_mat <- AdjacencyFromExpr(
    tissue_names = tissue_names,
    tissue_expr_file_names = tissue_expr_file_names,
    sd_quantile = sd_quantile,
    max_genes_per_tissue = max_genes_per_tissue,
    min_expr = min_expr,
    min_prop = min_prop,
    TS_power = TS_power,
    CT_power = CT_power,
    cor_method = cor_method
  )
  
  if (save_intermediates) {
    saveRDS(adj_mat, paste0(out_prefix, "_adjacency.rds"))
  }
  gc()
  
  message("Computing TOM (this may take a while)...")
  TOM_mat <- Cross_Tissue_TOM(adj_mat, TOMType = TOMType)
  rm(adj_mat); gc()
  
  message("Clustering on TOM…")
  # pass plot flag down (see tweak to Clusters_Table below)
  clusters_table <- Clusters_Table(TOM_mat, minClusterSize = minClusterSize, plot_heatmap = plot_heatmap)
  rm(TOM_mat); gc()
  
  clusters_details <- Clusters_Details(clusters_table, cluster_type_thr = cluster_type_thr)
  
  write.table(clusters_table,
              file = paste0(out_prefix, "_Clusters_table.txt"),
              quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  write.table(clusters_details,
              file = paste0(out_prefix, "_Clusters_details.txt"),
              quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  
  invisible(list(clusters_table = clusters_table, clusters_details = clusters_details))
}

# ------------------Cross-Tissue Clusters from gene expression data ------------------------

XWGCNA_Clusters_old <- function(tissue_names = NULL, tissue_expr_file_names = NULL, MV_sd_thr = 0.5, cluster_type_thr = 0.95, minClusterSize = 30) {
  
  adj_mat <- AdjacencyFromExpr(tissue_names, tissue_expr_file_names, MV_sd_thr)
  
  TOM_mat <- Cross_Tissue_TOM_old(adj_mat)
  
  clusters_table <- Clusters_Table(TOM_mat, minClusterSize = minClusterSize)
  
  clusters_details <- Clusters_Details(clusters_table, cluster_type_thr = cluster_type_thr)
  
  write.table(clusters_table, file = 'Clusters_table.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  write.table(clusters_details, file = 'Clusters_details.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  
  XWGCNA_clusters = list(clusters_table=clusters_table, clusters_details=clusters_details)
  
  return(XWGCNA_clusters)
}

# ------------------Check Network Scale free and Connectivity ------------------------

checkScaleFree <- function (k, nBreaks = 10, removeFirst = FALSE) 
{
  discretized.k = cut(k, nBreaks)
  dk = tapply(k, discretized.k, mean)
  p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
  breaks1 = seq(from = min(k), to = max(k), length = nBreaks + 1)
  hist1 = hist(k, breaks = breaks1, plot = FALSE, right = TRUE)
  dk2 = hist1$mids
  dk = ifelse(is.na(dk), dk2, dk)
  dk = ifelse(dk == 0, dk2, dk)
  p.dk = ifelse(is.na(p.dk), 0, p.dk)
  log.dk = as.vector(log10(dk))
  if (removeFirst) {
    p.dk = p.dk[-1]
    log.dk = log.dk[-1]
  }
  log.p.dk = as.numeric(log10(p.dk + 1e-09))
  lm1 = lm(log.p.dk ~ log.dk)
  lm2 = lm(log.p.dk ~ log.dk + I(10^log.dk))
  datout = data.frame(Rsquared.SFT = summary(lm1)$r.squared, slope.SFT = summary(lm1)$coefficients[2, 1], truncatedExponentialAdjRsquared = summary(lm2)$adj.r.squared)
  datout
}

R2_connectivity <- function(adj_mat) {
  #--------------------------For scale Free-----------------------------
  k <- rowSums(adj_mat)
  r2 <- checkScaleFree(k)$Rsquared.SFT
  r2_conn <- list(r2=r2, mean_conn=mean(k), median_conn=median(k), max_conn=max(k), min_conn=min(k))
  
  #------------------For connectons-----------------------
  tissue_indexed <- c(0, as.numeric(table(unlist(lapply(strsplit(colnames(adj_mat), split = '_'), function(x) {return(x[1])})))))
  for(i in 2:length(tissue_indexed)) tissue_indexed[i] <- tissue_indexed[i] <- sum(tissue_indexed[i:(i-1)])
  
  TS_conn_mat <- matrix(NA, nrow=(length(tissue_indexed)-1), ncol=4)
  colnames(TS_conn_mat) <- c('mean', 'median', 'max', 'min')
  
  CT_conn_mat <- matrix(NA, nrow=((((length(tissue_indexed)-1)*(length(tissue_indexed)-1))-(length(tissue_indexed)-1))/2), ncol=4)
  colnames(CT_conn_mat) <- c('mean', 'median', 'max', 'min')
  CT_counter <- 1
  
  for(i in 1:(length(tissue_indexed)-1)) {
    temp_adj_mat <- adj_mat[(tissue_indexed[i]+1):tissue_indexed[i+1], (tissue_indexed[i]+1):tissue_indexed[i+1]]
    k <- rowSums(temp_adj_mat)
    TS_conn_mat[i, ] <- c(mean(k), median(k), max(k), min(k))
    
    if(i < length(tissue_indexed)-1) {
      for(j in 1:(length(tissue_indexed)-1-i)) {
        temp_adj_mat <- adj_mat[(tissue_indexed[i]+1):tissue_indexed[i+1], (tissue_indexed[i+j]+1):tissue_indexed[i+j+1]]
        k <- rowSums(temp_adj_mat)
        CT_conn_mat[CT_counter, ] <- c(mean(k), median(k), max(k), min(k))
        CT_counter <- CT_counter + 1
      }
    }
    
  }
  
  
  scale_free_conn_list <- list(r2_conn=r2_conn, TS_conn_mat=TS_conn_mat, CT_conn_mat=CT_conn_mat)
  # save(scale_free_conn_list, file='Scale_free_conn_list.RData')
  
  
  R2_Conn <- list(r2=r2_conn$r2, mean_conn=r2_conn$mean_conn, TS_mean_conn=paste(TS_conn_mat [ ,'mean'], collapse = ', '), mean_TS_mean_conn=mean(TS_conn_mat [ ,'mean']),
                  max_TS_mean_conn=max(TS_conn_mat [ ,'mean']), min_TS_mean_conn=min(TS_conn_mat [ ,'mean']), CT_mean_conn=paste(CT_conn_mat[ ,'mean'], collapse = ', '),
                  mean_CT_mean_conn=mean(CT_conn_mat[ ,'mean']), max_CT_mean_conn=max(CT_conn_mat[ ,'mean']), min_CT_mean_conn=min(CT_conn_mat[ ,'mean']))
  
  return(R2_Conn)
}

# ------------------ Auto-pick beta (TS/CT) ------------------
.choose_power_from_pickSoft <- function(sft, targetR2 = 0.80) {
  df <- sft$fitIndices
  ok <- which(df$SFT.R.sq >= targetR2 & df$slope < 0)
  if (length(ok) > 0) {
    return(min(df$Power[ok]))
  } else {
    return(df$Power[which.max(df$SFT.R.sq)])
  }
}

.pick_TS_power_one <- function(expr, beta_grid, cor_method = "pearson", targetR2 = 0.80) {
  sft <- WGCNA::pickSoftThreshold(
    expr,
    powerVector = beta_grid,
    networkType = "unsigned",
    corFnc = "cor",
    corOptions = list(use = "pairwise.complete.obs", method = cor_method),
    verbose = 0
  )
  .choose_power_from_pickSoft(sft, targetR2 = targetR2)
}

.pick_CT_power_pair <- function(S_abs, beta_grid, targetR2 = 0.80) {
  best_beta <- NA_integer_; best_r2 <- -Inf
  for (b in beta_grid) {
    B <- S_abs^b
    k <- c(rowSums(B), colSums(B))
    r2 <- tryCatch(checkScaleFree(k)$Rsquared.SFT, error = function(e) NA_real_)
    if (!is.na(r2) && r2 >= targetR2) {
      return(b)  
    }
    if (!is.na(r2) && r2 > best_r2) {
      best_r2 <- r2; best_beta <- b
    }
  }
  best_beta
}

auto_pick_powers <- function(
    tissue_names,
    tissue_expr_file_names,
    sd_quantile = 0.90,
    max_genes_per_tissue = 4000,
    min_expr = 0,
    min_prop = 0.10,
    cor_method = "pearson",
    beta_grid = c(1:10, seq(12, 20, 2)),
    targetR2 = 0.80
){
  stopifnot(length(tissue_names) == length(tissue_expr_file_names))
  T <- length(tissue_names)
  
  # 1) טוענים ומסננים כמו ב-pipeline שלך
  expr_list   <- vector("list", T)
  donors_list <- vector("list", T)
  TS_fit_curves <- list()
  for (i in seq_len(T)) {
    X <- LoadExprData(
      tissue_name = tissue_names[i],
      tissue_file_name = tissue_expr_file_names[i],
      sd_quantile = sd_quantile,
      max_genes_per_tissue = max_genes_per_tissue
    )
    expr_list[[i]]   <- X
    donors_list[[i]] <- .aggregate_by_donor(X)
  }
  
  TS_per_tissue <- integer(T)
  for (i in seq_len(T)) {
    sft <- tryCatch(
      WGCNA::pickSoftThreshold(
        expr_list[[i]],
        powerVector = beta_grid,
        networkType = "unsigned",
        corFnc = "cor",
        corOptions = list(use = "pairwise.complete.obs", method = cor_method),
        verbose = 0
      ),
      error = function(e) { message("[TS beta pick] ", tissue_names[i], " error: ", e$message); NULL }
    )
    if (!is.null(sft)) {
      fi <- sft$fitIndices
      fi$tissue <- tissue_names[i]
      TS_fit_curves[[i]] <- fi
      TS_per_tissue[i] <- .choose_power_from_pickSoft(sft, targetR2 = targetR2)
    } else {
      TS_fit_curves[[i]] <- data.frame(Power = beta_grid, SFT.R.sq = NA, slope = NA, tissue = tissue_names[i])
      TS_per_tissue[i] <- 6L
    }
  }
  TS_fit_curves_df <- do.call(rbind, TS_fit_curves)
  TS_power <- as.integer(stats::median(TS_per_tissue))

  CT_betas <- c()
  CT_fit_list <- list()
  if (T >= 2) {
    for (i in 1:(T - 1)) {
      for (j in (i + 1):T) {
        di <- rownames(donors_list[[i]]); dj <- rownames(donors_list[[j]])
        common <- intersect(di, dj)
        if (length(common) < 3) {
          message(sprintf("[CT beta pick] %s-%s: |common donors|=%d → skip",
                          tissue_names[i], tissue_names[j], length(common)))
          next
        }
        Mi <- donors_list[[i]][common, , drop = FALSE]
        Mj <- donors_list[[j]][common, , drop = FALSE]
        S  <- abs(cor(Mi, Mj, use = "pairwise.complete.obs", method = cor_method))
        pair_id <- paste(tissue_names[i], tissue_names[j], sep = "||")
        curve_df <- data.frame(pair = pair_id, Power = beta_grid, Rsquared.SFT = NA_real_)
        for (b_ix in seq_along(beta_grid)) {
          b <- beta_grid[b_ix]
            B <- S^b
            k <- c(rowSums(B), colSums(B))
            r2 <- tryCatch(checkScaleFree(k)$Rsquared.SFT, error = function(e) NA_real_)
            curve_df$Rsquared.SFT[b_ix] <- r2
        }
        CT_fit_list[[length(CT_fit_list) + 1]] <- curve_df
        # choose beta for this pair
        best <- NA_integer_; best_r2 <- -Inf
        for (b_ix in seq_along(beta_grid)) {
          r2 <- curve_df$Rsquared.SFT[b_ix]
          b  <- beta_grid[b_ix]
          if (!is.na(r2) && r2 >= targetR2) { best <- b; break }
          if (!is.na(r2) && r2 > best_r2) { best_r2 <- r2; best <- b }
        }
        CT_betas <- c(CT_betas, best)
        message(sprintf("[CT beta pick] %s: beta=%d", pair_id, best))
      }
    }
  }
  if (length(CT_betas) == 0) {
    CT_power <- 3L
    message("[CT beta pick] no valid pairs → CT_power=3 (fallback)")
  } else {
    CT_power <- as.integer(stats::median(CT_betas))
  }
  CT_fit_curves_df <- if (length(CT_fit_list)) do.call(rbind, CT_fit_list) else
    data.frame(pair = character(), Power = integer(), Rsquared.SFT = numeric())

  list(
    TS_power = TS_power,
    CT_power = CT_power,
    TS_per_tissue = TS_per_tissue,
    CT_per_pair = CT_betas,
    TS_fit_curves = TS_fit_curves_df,
    CT_fit_curves = CT_fit_curves_df
  )
}

plot_TS_beta_curves <- function(ts_df, TS_power, TS_per_tissue, targetR2 = 0.80, out_prefix = "xwgcna_TS_beta") {
  if (!nrow(ts_df)) return(invisible(NULL))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; skipping TS beta plot.")
    return(invisible(NULL))
  }
  library(ggplot2)
  chosen_df <- data.frame(tissue = unique(ts_df$tissue),
                          chosen = TS_per_tissue)
  p <- ggplot(ts_df, aes(Power, SFT.R.sq, color = tissue)) +
    geom_line() + geom_point(size = 1.6) +
    geom_hline(yintercept = targetR2, linetype = "dashed", color = "black") +
    geom_vline(xintercept = TS_power, linetype = "dotted", color = "black") +
    geom_point(data = merge(ts_df, chosen_df, by = "tissue")[
      merge(ts_df, chosen_df, by = "tissue")$Power ==
        merge(ts_df, chosen_df, by = "tissue")$chosen, ],
      aes(Power, SFT.R.sq), color = "black", shape = 21, fill = "yellow", size = 3, stroke = 0.8) +
    labs(title = sprintf("TS power selection (median=%d)", TS_power),
         subtitle = sprintf("Horizontal dashed: target R^2=%.2f; dotted: median chosen", targetR2),
         y = "Scale-free fit R^2") +
    theme_bw()
  ggsave(paste0(out_prefix, ".png"), p, width = 8, height = 5, dpi = 150)
  p
}

plot_CT_beta_curves <- function(ct_df, CT_power, targetR2 = 0.80, out_prefix = "xwgcna_CT_beta") {
  if (!nrow(ct_df)) {
    message("No CT curves to plot.")
    return(invisible(NULL))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; skipping CT beta plot.")
    return(invisible(NULL))
  }
  library(ggplot2)
  p <- ggplot(ct_df, aes(Power, Rsquared.SFT, group = pair, color = pair)) +
    geom_line() + geom_point(size = 1.3) +
    geom_hline(yintercept = targetR2, linetype = "dashed", color = "black") +
    geom_vline(xintercept = CT_power, linetype = "dotted", color = "black") +
    labs(title = sprintf("CT power selection (median=%d)", CT_power),
         subtitle = sprintf("Horizontal dashed: target R^2=%.2f; dotted: median chosen", targetR2),
         y = "Scale-free fit R^2") +
    theme_bw() + theme(legend.position = "bottom", legend.key.size = unit(0.5, "lines"))
  ggsave(paste0(out_prefix, ".png"), p, width = 9, height = 6, dpi = 150)
  p
}

R2_connectivity_old <- function(adj_mat) {
  #--------------------------For scale Free-----------------------------
  k <- rowSums(adj_mat)
  r2 <- checkScaleFree(k)$Rsquared.SFT
  r2_conn <- list(r2=r2, mean_conn=mean(k), median_conn=median(k), max_conn=max(k), min_conn=min(k))
  
  #------------------For connectivty-----------------------
  tissue_indexed <- c(0, as.numeric(table(unlist(lapply(strsplit(colnames(adj_mat), split = '_'), function(x) {return(x[1])})))))
  
  for(i in 2:length(tissue_indexed)) tissue_indexed[i] <- tissue_indexed[i] <- sum(tissue_indexed[i:(i-1)])
  print(tissue_indexed)
  TS_conn_mat <- matrix(NA, nrow=(length(tissue_indexed)-1), ncol=4)
  colnames(TS_conn_mat) <- c('mean', 'median', 'max', 'min')
  
  CT_conn_mat <- matrix(NA, nrow=((((length(tissue_indexed)-1)*(length(tissue_indexed)-1))-(length(tissue_indexed)-1))/2), ncol=4)
  colnames(CT_conn_mat) <- c('mean', 'median', 'max', 'min')
  CT_counter <- 1
  
  for(i in 1:(length(tissue_indexed)-1)) {
    temp_adj_mat <- adj_mat[(tissue_indexed[i]+1):tissue_indexed[i+1], (tissue_indexed[i]+1):tissue_indexed[i+1]]
    k <- rowSums(temp_adj_mat)
    TS_conn_mat[i, ] <- c(mean(k), median(k), max(k), min(k))
    
    if(i < length(tissue_indexed)-1) {
      for(j in 1:(length(tissue_indexed)-1-i)) {
        temp_adj_mat <- adj_mat[(tissue_indexed[i]+1):tissue_indexed[i+1], (tissue_indexed[i+j]+1):tissue_indexed[i+j+1]]
        k <- rowSums(temp_adj_mat)
        CT_conn_mat[CT_counter, ] <- c(mean(k), median(k), max(k), min(k))
        CT_counter <- CT_counter + 1
      }
    }
    
  }
  
  
  scale_free_conn_list <- list(r2_conn=r2_conn, TS_conn_mat=TS_conn_mat, CT_conn_mat=CT_conn_mat)
  # save(scale_free_conn_list, file='Scale_free_conn_list.RData')
  
  print(scale_free_conn_list)
  
  R2_Conn <- list(r2=r2_conn$r2, mean_conn=r2_conn$mean_conn, TS_mean_conn=paste(TS_conn_mat [ ,'mean'], collapse = ', '), mean_TS_mean_conn=mean(TS_conn_mat [ ,'mean']),
                  max_TS_mean_conn=max(TS_conn_mat [ ,'mean']), min_TS_mean_conn=min(TS_conn_mat [ ,'mean']), CT_mean_conn=paste(CT_conn_mat[ ,'mean'], collapse = ', '),
                  mean_CT_mean_conn=mean(CT_conn_mat[ ,'mean']), max_CT_mean_conn=max(CT_conn_mat[ ,'mean']), min_CT_mean_conn=min(CT_conn_mat[ ,'mean']))
  print(R2_Conn)
  return(R2_Conn)
}