library(WGCNA)
library(org.Hs.eg.db)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
library(clusterProfiler)
library(purrr)
library(AnnotationDbi)

#####------ GO enrichment of modules-------#######

normalize_gene_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  # drop Ensembl version suffix
  sub("\\.\\d+$", "", x)
}

guess_keytype_and_keys <- function(ids) {
  ids <- normalize_gene_ids(ids)
  is_ensg <- grepl("^ENSG\\d+", ids)
  frac_ensg <- if (all(is.na(is_ensg))) 0 else mean(is_ensg, na.rm = TRUE)
  if (!is.nan(frac_ensg) && frac_ensg > 0.5) {
    list(keytype = "ENSEMBL", keys = ids)
  } else {
    list(keytype = "SYMBOL", keys = ids)
  }
}

# Helper: map to ENTREZ with validation

map_to_entrez <- function(keys_vec, keytype) {
  k <- unique(keys_vec[!is.na(keys_vec) & nzchar(keys_vec)])
  if (!length(k)) stop("Empty key vector after filtering NAs/empties.")
  # primary attempt
  valid <- tryCatch(AnnotationDbi::keys(org.Hs.eg.db, keytype = keytype),
                    error = function(e) character())
  use <- intersect(k, valid)
  if (length(use) > 0) {
    return(AnnotationDbi::select(org.Hs.eg.db, keys = use, columns = "ENTREZID", keytype = keytype))
  }
  # fallback: try the other common keytype
  alt <- if (identical(keytype, "ENSEMBL")) "SYMBOL" else "ENSEMBL"
  valid_alt <- tryCatch(AnnotationDbi::keys(org.Hs.eg.db, keytype = alt),
                        error = function(e) character())
  use_alt <- intersect(k, valid_alt)
  if (length(use_alt) > 0) {
    message(sprintf("Primary keytype '%s' failed; falling back to '%s'.", keytype, alt))
    return(AnnotationDbi::select(org.Hs.eg.db, keys = use_alt, columns = "ENTREZID", keytype = alt))
  }
  stop(sprintf("No valid keys for keytype='%s' or fallback. First few ids: %s",
               keytype, paste(utils::head(k, 5), collapse = ", ")))
}
`%||%` <- function(a, b) if (!is.null(a)) a else b

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
  stopifnot(all(c("Gene.Symbol", "Cluster.ID") %in% colnames(modules)))

  ids <- normalize_gene_ids(modules$Gene.Symbol)
  id_info <- guess_keytype_and_keys(ids)
  keytype <- id_info$keytype
  keys_vec <- id_info$keys

  
  # # Map ensembl IDs
  # ensembl = modules$Gene.Symbol
  
  # Load map of Ensembl -> ENTREX IDs
  # entrez_map = AnnotationDbi::select(org.Hs.eg.db, ensembl, "ENTREZID", "ENSEMBL")
  entrez_map = map_to_entrez(keys_vec, keytype = keytype)
  map_key_col <- if ("ENSEMBL" %in% colnames(entrez_map)) "ENSEMBL" else "SYMBOL"

  # entrez = entrez_map$ENTREZID[match(ensembl, entrez_map$ENSEMBL)]
  entrez = entrez_map$ENTREZID[match(keys_vec, entrez_map[[map_key_col]])]

  keep <- !is.na(entrez) & nzchar(entrez)
  if (!any(keep)) stop("No genes mapped to ENTREZ IDs.")

  GOenrichmentAnalysis(
    modules$Cluster.ID[keep],
    entrez[keep],
    organism="human",
    removeDuplicates=FALSE,
    nBestP=10,
    backgroundType = "givenInGO"
  )

  # # go_enrich = GOenrichmentAnalysis(between$Cluster.ID,
  # #                                   entrez,
  # #                                   organism="human",
  # #                                   #ontologies = "BP",
  # #                                   removeDuplicates=FALSE,
  # #                                   nBestP=10,
  # #                                   backgroundType = "givenInGO"
  # #                                   #pCut=0.05  
  # #  )
  
  # return(go_enrich)
  
}
  
#####------ Kegg enrichment of modules-------#######
enrichmentKegg = function(modules) {
  require(org.Hs.eg.db)
  stopifnot(all(c("Gene.Symbol", "Cluster.ID") %in% colnames(modules)))
  modules$Gene.ID.norm <- normalize_gene_ids(modules$Gene.Symbol)

  id_info <- guess_keytype_and_keys(modules$Gene.ID.norm)
  keytype <- id_info$keytype
  keys_vec <- id_info$keys

  entrez_map <- map_to_entrez(keys_vec, keytype = keytype)
  map_key_col <- if ("ENSEMBL" %in% colnames(entrez_map)) "ENSEMBL" else "SYMBOL"
  # Drop NAs, deduplicate key->ENTREZ (keep first if 1:many)
  ent_clean <- entrez_map[!is.na(entrez_map$ENTREZID) & nzchar(entrez_map$ENTREZID), c(map_key_col, "ENTREZID")]
  ent_clean <- ent_clean[!duplicated(ent_clean), , drop = FALSE]

  m <- merge(modules, ent_clean, by.x = "Gene.ID.norm", by.y = map_key_col)
  m <- m[!is.na(m$ENTREZID) & nzchar(m$ENTREZID), , drop = FALSE]
  if (!nrow(m)) {
    stop(sprintf(
      "No genes mapped to ENTREZ for KEGG. Example inputs: %s",
      paste(utils::head(unique(modules$Gene.ID.norm), 5), collapse = ", ")
    ))
  }

  kegg_enrich <- compareCluster(
    geneClusters = ENTREZID ~ Cluster.ID,
    data = m,
    fun = "enrichKEGG",
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

  # between1 = merge(modules, entrez_map, by.x=c("Gene.Symbol"), by.y=c("ENSEMBL"))
  # kegg_enrich = compareCluster(geneClusters = ENTREZID~Cluster.ID, data=between1, fun = "enrichKEGG" , organism = "hsa"  , keyType = "kegg", pvalueCutoff = 1, qvalueCutoff = 1)
  # kegg_enrich = setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  
  # return(kegg_enrich)
}


# ---------------------------------------------
Clusters_table = "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/Shared drives/is-SomekhLab/inter_tissue_coexpression/proccessed_data/ROSMAP_full_outputs/xwgcna_rosmap_autobeta_run4_Cluster_table.tsv"
Clusters_details <- read.delim("/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/Shared drives/is-SomekhLab/inter_tissue_coexpression/proccessed_data/ROSMAP_full_outputs/xwgcna_rosmap_autobeta_run4_Cluster_details.tsv",sep="\t", header=T)
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
            "kegg_rosmap_full.csv",
            sep=",",
            quote=FALSE,
            row.names=FALSE)

write.table(go_enrich_filter$bestPTerms$BP$enrichment,
            "go_bp_tab_rosmap_full.csv",
            sep=",",
            quote=FALSE,
            row.names=FALSE)

write.table(go_enrich_filter$bestPTerms$CC$enrichment,
            "go_cc_tab_rosmap_full.csv",
            sep=",",
            quote=FALSE,
            row.names=FALSE)

write.table(go_enrich_filter$bestPTerms$MF$enrichment,
            "go_mf_tab_rosmap_full.csv",
            sep=",",
            quote=FALSE,
            row.names=FALSE)


write.table(go_enrich_filter$enrichmentP,
            "go_pmat_rosmap_full.tsv",
            sep="\t")

cluster_summary <- data.frame(go_enrich)
fit<- plot(barplot(cluster_summary, showCategory=20))
fit






# ---------------------------------------------
clusters_table_path <-
  "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/Shared drives/is-SomekhLab/inter_tissue_coexpression/proccessed_data/ROSMAP_full_outputs/xwgcna_rosmap_autobeta_run5_Cluster_table.txt"
clusters_details <- read.delim(
  "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/Shared drives/is-SomekhLab/inter_tissue_coexpression/proccessed_data/ROSMAP_full_outputs/xwgcna_rosmap_autobeta_run5_Cluster_details.tsv",
  sep = "\t",
  header = TRUE
)
clusters_table <- read.delim(clusters_table_path, sep = "\t", header = TRUE)

go_enrich <- enrichmentGO(clusters_table)
go_enrich_filter <- filterCommonGOTerms(go_enrich, 50)
kegg_enrich <- enrichmentKegg(clusters_table)

# Extract top enrichment terms
n_clu <- max(clusters_table$Cluster.ID, na.rm = TRUE)

top_go_terms_bp <- sapply(seq_len(n_clu), function(k) {
  enrichment <- go_enrich_filter$bestPTerms$BP$enrichment
  idx <- enrichment$module == k
  top_terms <- enrichment$termName[idx][1:5]
  top_terms[!is.na(top_terms)]
})

top_go_terms_cc <- sapply(seq_len(n_clu), function(k) {
  enrichment <- go_enrich_filter$bestPTerms$CC$enrichment
  idx <- enrichment$module == k
  top_terms <- enrichment$termName[idx][1:5]
  top_terms[!is.na(top_terms)]
})

top_go_terms_mf <- sapply(seq_len(n_clu), function(k) {
  enrichment <- go_enrich_filter$bestPTerms$MF$enrichment
  idx <- enrichment$module == k
  top_terms <- enrichment$termName[idx][1:5]
  top_terms[!is.na(top_terms)]
})

go_tab <- data.frame(
  top_go_bp = sapply(top_go_terms_bp, paste, collapse = ";"),
  top_go_cc = sapply(top_go_terms_cc, paste, collapse = ";"),
  top_go_mf = sapply(top_go_terms_mf, paste, collapse = ";")
)

write.table(
  kegg_enrich@compareClusterResult,
  "kegg_rosmap_5000.csv",
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  go_enrich_filter$bestPTerms$BP$enrichment,
  "go_bp_tab_rosmap_5000.csv",
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  go_enrich_filter$bestPTerms$CC$enrichment,
  "go_cc_tab_rosmap_5000.csv",
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  go_enrich_filter$bestPTerms$MF$enrichment,
  "go_mf_tab_rosmap_5000.csv",
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  go_enrich_filter$enrichmentP,
  "go_pmat_rosmap_5000.tsv",
  sep = "\t"
)

# Optional plotting (clusterProfiler barplot expects enrichResult/compareClusterResult)
# barplot(kegg_enrich, showCategory = 20)
