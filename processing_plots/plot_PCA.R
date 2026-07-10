# PCA 2 -------------------------------------------------------------------

se_glycopep_RF <- readRDS("/data/glycobrain_docs/se_glycopep4008_vsn_RF_230426.rds")

df_glycopep_RF_meta <- get_df_wide(se_glycopep_RF) %>%
  dplyr::select(c(PG.Genes, PEP.PeptidePosition, PEP.StrippedSequence, epitope))

assay_p21_glycopep <- assay(se_glycopep_RF)

assay_p21_glycopep <- cbind(assay_p21_glycopep, df_glycopep_RF_meta)

assay_p21_glycopep %<>%
  mutate(end.pos = nchar(PEP.StrippedSequence) + as.numeric(PEP.PeptidePosition)-1) %>%
  dplyr::filter(epitope == "1Tn") %>%
  mutate(feature.name = paste0(PG.Genes, "_", PEP.PeptidePosition, "-", end.pos, "_", epitope)) %>%
  dplyr::select(-c(PG.Genes, PEP.PeptidePosition, PEP.StrippedSequence, epitope, end.pos)) %>%
  as.data.frame() %>%
  { mat <- as.matrix(.)
  rownames(mat) <- .$feature.name
  mat } %>%
  t() %>%
  as.data.frame() %>%  
  filter(row_number() < n()) %>%
  mutate(row.names = rownames(.)) %>%
  separate(row.names, into = "Region", sep = "_", extra = "drop")


assay_p21_glycopep <- type.convert(assay_p21_glycopep, as.is = TRUE)

res.pca <-  prcomp(assay_p21_glycopep[1:(length(assay_p21_glycopep)-1)])

eig_matrix <- get_eigenvalue(res.pca)

# 2. Extract variance for Dim 1 and Dim 2, then round the numbers
pc1_var <- round(eig_matrix[1, 2], 1) # Row 1, Column 2 (Variance %)
pc2_var <- round(eig_matrix[2, 2], 1) # Row 2, Column 2 (Variance %)

region_colors <- c("Cortex"="#E31A1C",
                   "Olfactory"="#1F78B4",
                   "Cerebellum"="#B2DF8A",
                   "Thalamus"="#1B9E77",
                   "Hippocampus"="#666666",
                   "Midbrain"="#A6CEE3",
                   "Medulla"="#FF61C3",
                   "Hypothalamus"="#FF7F00",
                   "Pons"="#6A3D9A")

#pca_all / (pca_1T | pca_1Tn)

pca_1Tn <- fviz_pca_biplot(res.pca,
                           habillage= assay_p21_glycopep$Region,
                           title = "PCA: Glycopeptides, 1Tn only",
                           pointshape = 19,
                           pointsize = 2,
                           label = "var",
                           labelsize = 3,      
                           mean.point = FALSE,
                           repel = TRUE,
                           addEllipses=TRUE,
                           ellipse.type="confidence",
                           ellipse.level=0.95,
                           ellipse.alpha = 0.33,
                           #ellipse.border.remove = TRUE,
                           legend.title = "Region",
                           axes.linetype = NA,
                           select.var = list(contrib = 15),
                           geom.var = c("arrow", "text"),
                           arrowsize = 0.25,
                           col.var = "black",
                           ggtheme = theme_bw()) +
  scale_color_manual(values = region_colors) +
  scale_fill_manual(values = region_colors) +
  # guides(color = "none",
  guides(color = guide_legend(override.aes = list(shape = 19, size = 3, fill = NA ,linetype = NA)),
         fill = "none") + 
  labs(x = paste0("PC1 (", pc1_var, "%)"), 
       y = paste0("PC2 (", pc2_var, "%)"))

#pca_all pca_1T pca_1Tn


pca_all / (pca_1T | pca_1Tn) + 
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a') &               
  theme(plot.tag = element_text(face = "bold", size = 14)) 
