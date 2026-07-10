
# DEA Proteome ------------------------------------------------------------

regions <- c("Cerebellum", "Medulla", "Pons", "Midbrain", "Hypothalamus", "Thalamus", "Hippocampus", "Cortex", "Olfactory")
all_unique_region_pairs <- combn(regions, 2, FUN = paste, collapse =  "_vs_")


diff_prot_vsn_pg_RF_BH_all <- test_diff(prot_vsn_pg_RF, type = "manual", test = all_unique_region_pairs, fdr.type = "BH")
diff_prot_vsn_pg_RF_BH_all_curve <- add_rejections(diff_prot_vsn_pg_RF_BH_all, thresholdmethod = "curve", curvature  = 2, x0_fold = 1)

df_diff_pg_RF_BH_all_curve_StatRes <- get_df_wide(diff_prot_vsn_pg_RF_BH_all_curve)
df_stat_res_p21_proteome_all_curve_demult <- df_diff_pg_RF_BH_all_curve_StatRes %>% separate_longer_delim(c(PG.ProteinGroups, PG.Genes), delim = ";")


saveRDS(df_stat_res_p21_proteome_all_curve_demult, 
        "/data/glycobrain_docs/df_stat_res_p21_proteome_all_curve_demult_230426.rds")
