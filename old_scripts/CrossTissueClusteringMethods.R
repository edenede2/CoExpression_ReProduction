# This code uses functions from hustal/X-WGCNA Repository
# Repository URL: https://github.com/hustal/X-WGCNA
# Author(s): Husain Ahammad Talukdar
# License: GNU General Public License v3.0 (GPLv3)
# The original code has been modified for this project.
# The modified code follows the same GPLv3 license.


library(WGCNA)


# -------------------Load gene expression data  --------------
LoadExprData<-function(tissue_name, tissue_file_name, MV_sd_thr) {
  
  #rows indicate samples and columns indicate gene symbols
  if(substr(tissue_file_name, start = nchar(tissue_file_name)-3, stop = nchar(tissue_file_name)) == ".RData") {
    
    datExpr <- get(load(tissue_file_name))
  } else if(substr(tissue_file_name, start = nchar(tissue_file_name)-3, stop = nchar(tissue_file_name)) == ".csv") {
    
    datExpr <-read.csv(tissue_file_name)
    datExpr <- as.matrix(datExpr)
    rownames(datExpr) <- datExpr[ ,1]
    datExpr <- datExpr[ ,-1]
    datExpr <- apply(datExpr, c(1,2), as.numeric)
  } else stop("Unsupported input file !!!")
  
  colnames(datExpr) <- paste(tissue_name,"_",colnames(datExpr),sep="")
 
  return(datExpr)
  
}

# -------------------Adjacency matrix from all tissue expression data ----------------------------------------- 
AdjacencyFromExpr <- function(tissue_names = NULL, tissue_expr_file_names = NULL, MV_sd_thr = 0.5, TS_power = 6, CT_power = 3, cor_method = "pearson") {

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

# ------------------TOM from Adjacency matrix ----------------------- 

Cross_Tissue_TOM <- function(adj_mat) {
  
  matrix_product <- adj_mat %*% t(adj_mat)
  print("Matrix Multiplication done...")
  rsums<-rowSums(adj_mat)
  D_val<-matrix(rsums,nrow=nrow(adj_mat),ncol=ncol(adj_mat))
  
  total_col <- ncol(D_val)
  if(total_col <= 5000) col_limit <- total_col else col_limit <- 5000
  start_col <- 1
  end_col <- col_limit
  
  tmp_D <- matrix(0, nrow = total_col, ncol = total_col)
  
  repeat {
    tmp_D[ ,start_col:end_col] <- pmin(D_val[ ,start_col:end_col], t(D_val)[ ,start_col:end_col])
    if(end_col==total_col) break
    start_col <- end_col + 1
    end_col <- end_col + col_limit
    if(end_col > total_col) end_col <- total_col
  }  
  
  rm(D_val)
  rm(rsums)
  gc()
  print("Matrix pmin done...")
  
  TOM_mat <- matrix(NA, nrow=nrow(adj_mat), ncol=ncol(adj_mat))
  rownames(TOM_mat) <- rownames(adj_mat)
  colnames(TOM_mat) <- colnames(adj_mat)
  
  TOM_mat <- (adj_mat + matrix_product) / (tmp_D + 1 - adj_mat)
  
  print("Tom done...")
  
  return(TOM_mat) 
  
}

#----------------- Plot the resulting clustering tree (dendrogram)-----------------
network_heatmap <- function(restGenes, dynamicColors) { 
  
  library(gplots)
  myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
  adj_mat_filter <- AdjacencyFromExpr(c("Adipose", "Muscle", "Brain"), c("Adipose - Subcutaneous_young.csv","Muscle - Skeletal_young.csv", "Brain - Cortex_young.csv"), 0.5)
  
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

Clusters_Table <- function(TOM_mat, minClusterSize = 30) {
  
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

# ------------------Cross-Tissue Clusters from gene expression data ------------------------

XWGCNA_Clusters <- function(tissue_names = NULL, tissue_expr_file_names = NULL, MV_sd_thr = 0.5, cluster_type_thr = 0.95, minClusterSize = 30) {
  
  adj_mat <- AdjacencyFromExpr(tissue_names, tissue_expr_file_names, MV_sd_thr)
  
  TOM_mat <- Cross_Tissue_TOM(adj_mat)
  
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


