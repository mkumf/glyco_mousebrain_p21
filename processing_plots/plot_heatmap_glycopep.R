# heatmap glycopep --------------------------------------------------------


se_glycopep_RF <- readRDS("/data/glycobrain_docs/se_glycopep4008_vsn_RF_230426.rds")

df_glycopep_onlynum <- assay(se_glycopep_RF)

mat_scaled = t(scale(t(as.matrix(df_glycopep_onlynum))))



annotation_col = glycopep_metadata

annotation_col %<>%
  mutate(Region = condition) %>%
  dplyr::select(-c(colname_raw, replicate, timepoint, Batch, condition)) %>%
  column_to_rownames("label")


ann_colors = list(Region = c("Cortex"="#E31A1C",
                             "Olfactory"="#1F78B4",
                             "Cerebellum"="#B2DF8A",
                             "Thalamus"="#1B9E77",
                             "Hippocampus"="#666666",
                             "Midbrain"="#A6CEE3",
                             "Medulla"="#FF61C3",
                             "Hypothalamus"="#FF7F00",
                             "Pons"="#6A3D9A"))

set.seed(123)
ComplexHeatmap::Heatmap(mat_scaled,
                        column_title = "Glycopeptides",
                        column_title_gp = gpar(fontface = "bold"),
                        #name = "log2 abundance\nscaled",  seq(-2, 2, length.out = 11), 
                        col =  colorRamp2(breaks = c(-3.0 ,-1.2 ,-0.8 ,-0.4 ,-0.2 , 0.0 , 0.2 , 0.4 , 0.8 , 1.2 , 3), 
                                          colors = rev(brewer.pal(11, "BrBG"))),
                        cluster_columns = TRUE,
                        clustering_distance_columns = "euclidean",
                        clustering_method_columns = "complete",
                        # column_split = 4,
                        cluster_rows = FALSE,
                        #        row_split = row_split, 
                        row_km = 4, 
                        #border = TRUE, 
                        #cluster_row_slices = TRUE,
                        row_km_repeats = 100,
                        show_row_dend = FALSE,
                        # row_dend_width = unit(10, "mm"),
                        #        row_km = 8, # equivalent to cutree_rows = 8
                        show_row_names = FALSE,
                        row_title_rot = 0,
                        #        row_names_gp = gpar(fontsize = 1),
                        column_names_gp = gpar(fontsize = 8),
                        
                        column_dend_height = unit(5, "mm"),
                        #        top_annotation = col_colors_compartments,
                        #        left_annotation = row_ha,
                        use_raster = FALSE,
                        top_annotation = HeatmapAnnotation(df = annotation_col, 
                                                           col = ann_colors,
                                                           annotation_name_side = "left",
                                                           border = TRUE),
                        heatmap_legend_param = list(
                          title = "abundance\nlog2, scaled",
                          legend_direction = "vertical",
                          border = "black"))