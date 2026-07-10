# missval heatmap ---------------------------------------------------------

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

df <-  assay(se_glycopep_filt) %>% data.frame(.)
missval <- df[apply(df, 1, function(x) any(is.na(x))), ]
missval <- ifelse(is.na(missval), 0, 1)
ht2 = ComplexHeatmap::Heatmap(missval, 
                              column_title = "Glycopeptides missing value patterns",
                              col = c("white", "black"),
                              column_names_side = "bottom", 
                              show_row_names = FALSE,
                              show_column_names = TRUE, 
                              cluster_rows = TRUE,
                              cluster_columns = TRUE,
                              clustering_distance_columns = 'binary',
                              clustering_distance_rows = 'binary',
                              # clustering_distance_columns = 'euclidean',
                              #clustering_distance_rows = 'euclidean',
                              clustering_method_columns = "ward.D2",
                              clustering_method_rows = "ward.D2",
                              name = "Glycopeptide\ndetected",
                              column_names_gp = gpar(fontsize = 8),
                              top_annotation = HeatmapAnnotation(df = annotation_col, 
                                                                 col = ann_colors,
                                                                 annotation_name_side = "left",
                                                                 border = TRUE),
                              heatmap_legend_param = list(at = c(0,1), 
                                                          labels = c("absent", "valid"),
                                                          legend_direction = "vertical",
                                                          border = "black"),
                              rect_gp = gpar(type = "none"), # Deaktiviert das Standard-Zeichnen der Kacheln
                              layer_fun = function(j, i, x, y, w, h, fill) {
                                grid.rect(x = x, 
                                          y = y, 
                                          width = w,
                                          height = h, 
                                          gp = gpar(fill = fill,
                                                    col = fill,
                                                    lex = 0,
                                                    linejoin = "mitre")
                                )},
                              use_raster = FALSE)

draw(ht2, merge_legends = TRUE)