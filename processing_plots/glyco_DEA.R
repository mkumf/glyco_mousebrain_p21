## DEA ---------------------------------------------------------------------

regions <- c("Cerebellum", "Medulla", "Pons", "Midbrain", "Hypothalamus", "Thalamus", "Hippocampus", "Cortex", "Olfactory")
all_unique_region_pairs <- combn(regions, 2, FUN = paste, collapse =  "_vs_")

#diff_se_glycopep_RF_allpairs <- test_diff(se_glycopep_RF, type = "manual", test  = all_unique_region_pairs, fdr.type = "BH")


diff_se_glycopep_RF_allpairs <- test_diff(se_glycopep_RF, type = "manual", test  = comparisions_old, fdr.type = "BH")

diff_se_glycopep_RF_allpairs_curve <- add_rejections(diff_se_glycopep_RF_allpairs, thresholdmethod = "curve", curvature  = 2, x0_fold = 2)

stat_res_p21_glycopep <- get_df_wide(diff_se_glycopep_RF_allpairs_curve)