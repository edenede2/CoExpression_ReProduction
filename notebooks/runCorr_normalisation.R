source("/Users/edeneldar/CoExpression_ReProduction/old_scripts/ClusteringBuilding.R")
source("/Users/edeneldar/CoExpression_ReProduction/notebooks/correlation_normalisation_2.R")

res <- compare_correlation_methods(
    tissue_names = c("AdiposeSubcutaneous","AdiposeVisceralOmentum","ArteryTibial","BrainCerebellum","BrainCortex","EsophagusMucosa", "EsophagusMuscularis", "Liver", "MuscleSkeletal", "NerveTibial","SkinNotSunExposedSuprapubic", "SkinSunExposedLowerLeg","Thyroid"),  
    tissue_expr_file_names = c(    
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_AdiposeSubcutaneous.csv",    
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_AdiposeVisceralOmentum.csv",    
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_ArteryTibial.csv",
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_BrainCerebellum.csv",
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_BrainCortex.csv",
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_EsophagusMucosa.csv",
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_EsophagusMuscularis.csv",
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_Liver.csv",
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_MuscleSkeletal.csv",
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_NerveTibial.csv",
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_SkinNotSunExposedSuprapubic.csv",
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_SkinSunExposedLowerLeg.csv",
        "/Users/edeneldar/Library/CloudStorage/GoogleDrive-edenede2@gmail.com/My Drive/Shaked's src/all/Tissues/GTEX_outputs/full/GTEX_full_Thyroid.csv"
    ),  perm_mode   = "analytic",  n_ref_mode  = "global",  n_ref       = "max",  make_plots  = TRUE,  out_prefix  = "gtex_cmp_global_3000",max_genes_per_tissue=2500,  fit_scale_free = TRUE,  sft_target_R2 = 0.8, sft_nBreaks = 12,  sft_power_grid = c(1:10))
