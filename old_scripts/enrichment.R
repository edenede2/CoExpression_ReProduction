library(WGCNA)
library(org.Hs.eg.db)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
library(clusterProfiler)
library(purrr)

#####------ GO enrichment of modules-------#######

# Remove GO terms that are prevalent across modules
removePrevalentTerms = function(enrich_tab, max_preval) {
  mod_counts = table(enrich_tab$termID)
  include_idx = mod_counts < max_preval
  message("Excluding GO terms: ", sum(!include_idx))
  include_terms = names(which(include_idx))
  return(enrich_tab[enrich_tab$termID %in% include_terms, ])
}

# Filter out GO terms from WGCNA GO enrichment analysis
filterCommonGOTerms = function(go, max_preval=50) {
  go$bestPTerms$BP$enrichment = removePrevalentTerms(go$bestPTerms$BP$enrichment, max_preval)
  go$bestPTerms$CC$enrichment = removePrevalentTerms(go$bestPTerms$CC$enrichment, max_preval)
  go$bestPTerms$MF$enrichment = removePrevalentTerms(go$bestPTerms$MF$enrichment, max_preval)
  return(go)
}


enrichmentGO = function(modules) {
  require(org.Hs.eg.db)
  require(WGCNA)
  # Map ensembl IDs
  ensembl = modules$Gene.Symbol
  
  # Load map of Ensembl -> ENTREX IDs
  entrez_map = AnnotationDbi::select(org.Hs.eg.db, ensembl, "ENTREZID", "ENSEMBL")
  
  entrez = entrez_map$ENTREZID[match(ensembl, entrez_map$ENSEMBL)]
  
  go_enrich = GOenrichmentAnalysis(between$Cluster.ID,
                                    entrez,
                                    organism="human",
                                    #ontologies = "BP",
                                    removeDuplicates=FALSE,
                                    nBestP=10,
                                    backgroundType = "givenInGO"
                                    #pCut=0.05  
   )
  
  return(go_enrich)
  
}
  
#####------ Kegg enrichment of modules-------#######
  enrichmentKegg = function(modules) {
 
  between1 = merge(modules, entrez_map, by.x=c("Gene.Symbol"), by.y=c("ENSEMBL"))
  kegg_enrich = compareCluster(geneClusters = ENTREZID~Cluster.ID, data=between1, fun = "enrichKEGG" , organism = "hsa"  , keyType = "kegg", pvalueCutoff = 1, qvalueCutoff = 1)
  kegg_enrich = setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  
  return(kegg_enrich)
}


# ---------------------------------------------
Clusters_table = "Clusters_table_old.txt"
Clusters_details <- read.delim("Clusters_details_old.txt",sep="\t", header=T)
Clusters_table <- read.delim(Clusters_table,sep="\t", header=T)

go_enrich = enrichmentGO(Clusters_table)
go_enrich_filter = filterCommonGOTerms(go_enrich, 50)
kegg_enrich = enrichmentKegg(Clusters_table)

# Extract top enrichment terms
top_go_terms_bp = sapply(1:max(Clusters_table$Cluster.ID), function(k) {
  enrichment = go_enrich_filter$bestPTerms$BP$enrichment
  idx = enrichment$module == k
  top_terms = enrichment$termName[idx][1:5]
  top_terms = top_terms[!is.na(top_terms)]
  return(top_terms)
})

top_go_terms_cc = sapply(1:max(Clusters_table$Cluster.ID), function(k) {
  enrichment = go_enrich_filter$bestPTerms$CC$enrichment
  idx = enrichment$module == k
  top_terms = enrichment$termName[idx][1:5]
  top_terms = top_terms[!is.na(top_terms)]
  return(top_terms)
})

top_go_terms_mf = sapply(1:max(between$Cluster.ID), function(k) {
  enrichment = go_enrich_filter$bestPTerms$MF$enrichment
  idx = enrichment$module == k
  top_terms = enrichment$termName[idx][1:5]
  top_terms = top_terms[!is.na(top_terms)]
  return(top_terms)
})


go_tab = data.frame(
  top_go_bp=sapply(top_go_terms_bp, paste, collapse=";"),
  top_go_cc=sapply(top_go_terms_cc, paste, collapse=";"),
  top_go_mf=sapply(top_go_terms_mf, paste, collapse=";")
)



write.table(kegg_enrich@compareClusterResult,
            "kegg_old.csv",
            sep=",",
            quote=FALSE,
            row.names=FALSE)

write.table(go_enrich_filter$bestPTerms$BP$enrichment,
            "go_bp_tab_old.csv",
            sep=",",
            quote=FALSE,
            row.names=FALSE)

write.table(go_enrich_filter$bestPTerms$CC$enrichment,
            "go_cc_tab_old.csv",
            sep=",",
            quote=FALSE,
            row.names=FALSE)

write.table(go_enrich_filter$bestPTerms$MF$enrichment,
            "go_mf_tab_old.csv",
            sep=",",
            quote=FALSE,
            row.names=FALSE)


write.table(go_enrich_filter$enrichmentP,
            "go_pmat_old.tsv",
            sep="\t")

cluster_summary <- data.frame(go_enrich)
fit<- plot(barplot(cluster_summary, showCategory=20))
fit






