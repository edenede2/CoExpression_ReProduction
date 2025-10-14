
# Determines power of cross-tissue modules from reshaped gene expression matrix (tissue_gene x sample)
# Calculates optimal beta values using the scale-free R^2 criterion from WGCNA.
# Finds lowest beta above threhshold

# This code uses functions from skoplev/starnet Repository
# Repository URL: https://github.com/skoplev/starnet
# Author(s): Simon Koplev
# License: GNU General Public License v3.0 (GPLv3)
# The original code has been modified for this project.
# The modified code follows the same GPLv3 license.

opts = list()

# Picking a soft power threshold to get scale-free correlation networks from each tissue alone

opts$powers = seq(0.5, 20, length.out=20)

opts$n_breaks = 50  # discrete bins of connectivity distribution

# opts$abs_cor_min = 0.2  # minimum correlation to be considered
# cor_quant = 0.95


# Parse and check user input
# -------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)  # R-style positional arguments
if (length(args) != 1) {
  stop("Invalid arguments. USAGE: determinePower.R <path to cross-tissue expression object>")
}
opts$emat_file = args[1]

if (!exists("args")) {
   
  data_dir = "adipose-muscle-brain_cortex/var_filter_5000"
  
  # Load imputed recast gene expression matrix
  opts$emat_file =("adipose-muscle-brain_cortex/var_filter_5000/adj_mat.Rds")
}


# Convert table from reshape2
if (exists("expr_recast")) {
  mat = expr_recast[, 3:ncol(expr_recast)]
  mat = data.matrix(mat)
  row_meta = expr_recast[, 1:2]
  row_meta = as.data.frame(row_meta)
}



# Set up R environment
# -------------------------------------------------------------
library(data.table)
library(Matrix)
library(compiler)
enableJIT(3)
library(WGCNA) # Load  WGCNA
options(stringsAsFactors = FALSE)


# Evaluate sequence of soft network cutoffs within tissues
# ------------------------------------------------------------
con_eval = list()
tissue_names = c("Adipose", "Muscle", "Brain")
tissue_expr_file_names = c("Adipose - Subcutaneous_young.csv","Muscle - Skeletal_young.csv", "Brain - Cortex_young.csv")
MV_sd_thr =  0.5
total_tissue <- length(tissue_names)

vector_expr<-vector(mode="list", length=total_tissue)

##----This is to track tissue in whole adjacency matrix--------------------
tissue_index_adj<-vector(mode="integer", length=(total_tissue+1))
rc_names<-c()

for(i in 1:total_tissue) {
  TS_expr_data<-LoadExprData(tissue_names[i], tissue_expr_file_names[i], MV_sd_thr)
  vector_expr[[i]]<-TS_expr_data
 
  tissue_index_adj[i+1]<- tissue_index_adj[i]+ncol(vector_expr[[i]])
  rc_names<-c(rc_names,colnames(vector_expr[[i]]))
  
  message("Estimating power for: ", tissue_names[i])
  power_thresh = pickSoftThreshold(TS_expr_data, 
                                  powerVector=opts$powers, 
                                  verbose = 5)
  
  con_eval[[i]] = power_thresh
  
  gc()
  
}


# Test optimal beta for combinations of tissues. Uses cross-correlations only.
# Between-tissues
# ----------------------------------------------------------
con_eval_pairs = list()

k=1
for(i in 1:(total_tissue-1)) {
  for(j in (i+1):total_tissue) {

    common_Samples <- intersect(rownames(vector_expr[[i]]),rownames(vector_expr[[j]]))

    paired_tissue = cbind(vector_expr[[i]][common_Samples,],vector_expr[[j]][common_Samples,])
    
    message("calculating cross-tissue correlation: ", tissue_names[i], ", ", tissue_names[j])
    
    power_thresh = pickSoftThreshold(paired_tissue, 
                                      powerVector=opts$powers, 
                                     verbose = 5)
    
    con_eval_pairs[[k]] = power_thresh
    k = k +1
    
    gc()
    
  }
}


names(con_eval) = c("Adipose", "Muscle", "Brain")
names(con_eval_pairs) = c("Adipose-Muscle", "Adipose-Brain", "Muscle-Brain")

dir.create("determinePower/output")
save(opts, con_eval, con_eval_pairs, file="con_eval_v2auto.RData")
