# library GO treeplot -----------------------------------------------------


mBrain_DIA_Library_entire_raw <- read.xlsx("/data/mBrain_DIA_Library_Cumulative_Statistics_reduced.xlsx",
                                           sheet = "Entire Brain non red")

mBrain_DIA_Library_compartments_raw <- read.xlsx("/data/mBrain_DIA_Library_Cumulative_Statistics_reduced.xlsx",
                                                 sheet = "Compartments_details")

entire_uniprot <- mBrain_DIA_Library_entire_raw %>%
  dplyr::select(UniProtIds) %>%
  tidyr::separate_rows(UniProtIds, sep = ";") %>%
  mutate(UniProtIds = str_trim(UniProtIds, side = "left")) %>%
  distinct()

compartments_uniprot <- mBrain_DIA_Library_compartments_raw %>%
  dplyr::select(Master.Protein.Accessions) %>%
  tidyr::separate_rows(Master.Protein.Accessions, sep = ";") %>%
  mutate(Master.Protein.Accessions = str_trim(Master.Protein.Accessions, side = "left")) %>%
  distinct()

colnames(entire_uniprot) <- names(compartments_uniprot)
glyco_all_uniprot <- data.frame(Uniprot = unique(rbind(entire_uniprot, compartments_uniprot)))
colnames(glyco_all_uniprot) <- "PG.ProteinAccessions"


protein_uniprot <- read.xlsx("/data/aging_p21/FINAL/2025_gBrain_Aging_p21_peptidome_NotNorm_HybridDIA_Report.xlsx",
                             cols = 1) %>% distinct()

protein_uniprot_unique <- protein_uniprot %>%
  dplyr::select(PG.ProteinAccessions) %>%
  tidyr::separate_rows(PG.ProteinAccessions, sep = ";") %>%
  mutate(PG.ProteinAccessions = str_trim(PG.ProteinAccessions, side = "left")) %>%
  distinct()

protein_background_comb <- rbind(protein_uniprot_unique, glyco_all_uniprot)

ego_mf <- clusterProfiler::enrichGO(gene = glyco_all_uniprot$PG.ProteinAccessions,
                                    OrgDb         = org.Mm.eg.db,
                                    universe = protein_background_comb$PG.ProteinAccessions,
                                    keyType       = "UNIPROT",
                                    ont           = "MF",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.01,
                                    qvalueCutoff  = 0.05,
                                    readable = TRUE)
ego_mf_simp <- clusterProfiler::simplify(ego_mf, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")
ego_mf_simp_termsim <- enrichplot::pairwise_termsim(ego_mf_simp)



ego_bp <- clusterProfiler::enrichGO(gene = glyco_all_uniprot$PG.ProteinAccessions,
                                    OrgDb         = org.Mm.eg.db,
                                    universe = protein_background_comb$PG.ProteinAccessions,
                                    keyType       = "UNIPROT",
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.01,
                                    qvalueCutoff  = 0.05,
                                    readable = TRUE)
ego_bp_simp <- clusterProfiler::simplify(ego_bp, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")
ego_bp_simp_termsim <- enrichplot::pairwise_termsim(ego_bp_simp)



ego_cc <- clusterProfiler::enrichGO(gene = glyco_all_uniprot$PG.ProteinAccessions,
                                    OrgDb         = org.Mm.eg.db,
                                    universe = protein_background_comb$PG.ProteinAccessions,
                                    keyType       = "UNIPROT",
                                    ont           = "CC",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.01,
                                    qvalueCutoff  = 0.05,
                                    readable = TRUE)
ego_cc_simp <- clusterProfiler::simplify(ego_cc, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")
ego_cc_simp_termsim <- enrichplot::pairwise_termsim(ego_cc_simp)



p_mf <- enrichplot::treeplot(ego_mf_simp_termsim,
                             showCategory = 20,
                             cluster_method = "average",
                             cladelab_offset = 8, 
                             tiplab_offset = .3,
                             fontsize_tiplab = 2.5,
                             fontsize_cladelab = 3,
                             label_format = 20) + 
  ggtree::hexpand(.2) +
  labs(title = "          Molecular function")
p_mf$data$label <- stringr::str_wrap(p_mf$data$label, width = 23)

for (i in seq_along(p_mf$layers)) {
  if ("geom" %in% names(p_mf$layers[[i]]) && inherits(p_mf$layers[[i]]$geom, "GeomText")) {
    p_mf$layers[[i]]$aes_params$lineheight <- 0.8
  }
}

p_mf$layers[[1]]$show.legend <- FALSE
p_mf$layers[[2]]$show.legend <- FALSE

p_mf <- p_mf + theme(legend.position = "bottom",
                     legend.box = "horizontal",
                     legend.direction = "horizontal",
                     legend.title = element_text(size = 7), 
                     legend.text = element_text(size = 7)) +
  scale_size_continuous(breaks = scales::breaks_extended(n = 3)) + 
  guides(color = guide_colorbar(title.position = "top",
                                barwidth = 8,
                                barheight = 0.5,
                                reverse = TRUE,
                                theme = theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))))




p_bp <- enrichplot::treeplot(ego_bp_simp_termsim,
                             showCategory = 20,
                             cluster_method = "average",
                             cladelab_offset = 8, 
                             tiplab_offset = .3,
                             fontsize_tiplab = 2.5,
                             fontsize_cladelab = 3,
                             label_format = 20) + 
  ggtree::hexpand(.2) +
  labs(title = "          Biological process")

p_bp$data$label <- stringr::str_wrap(p_bp$data$label, width = 23)

for (i in seq_along(p_bp$layers)) {
  if ("geom" %in% names(p_bp$layers[[i]]) && inherits(p_bp$layers[[i]]$geom, "GeomText")) {
    p_bp$layers[[i]]$aes_params$lineheight <- 0.8
  }
}

p_bp$layers[[1]]$show.legend <- FALSE
p_bp$layers[[2]]$show.legend <- FALSE

p_bp <- p_bp + theme(legend.position = "bottom",
                     legend.box = "horizontal",
                     legend.direction = "horizontal",
                     legend.title = element_text(size = 7), 
                     legend.text = element_text(size = 7)) +
  scale_size_continuous(breaks = scales::breaks_extended(n = 3)) + 
  guides(color = guide_colorbar(title.position = "top",
                                barwidth = 8,
                                barheight = 0.5,
                                reverse = TRUE,
                                theme = theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))))





p_cc <- enrichplot::treeplot(ego_cc_simp_termsim,
                             showCategory = 20,
                             cluster_method = "average",
                             cladelab_offset = 8, 
                             tiplab_offset = .3,
                             fontsize_tiplab = 2.5,
                             fontsize_cladelab = 3,
                             label_format = 20) + 
  ggtree::hexpand(.2) +
  labs(title = "          Cellular compartment") 

p_cc$data$label <- stringr::str_wrap(p_cc$data$label, width = 23)

for (i in seq_along(p_cc$layers)) {
  if ("geom" %in% names(p_cc$layers[[i]]) && inherits(p_cc$layers[[i]]$geom, "GeomText")) {
    p_cc$layers[[i]]$aes_params$lineheight <- 0.8
  }
}

p_cc$layers[[1]]$show.legend <- FALSE
p_cc$layers[[2]]$show.legend <- FALSE

p_cc <- p_cc + theme(legend.position = "bottom",
                     legend.box = "horizontal",
                     legend.direction = "horizontal",
                     legend.title = element_text(size = 7), 
                     legend.text = element_text(size = 7)) +
  scale_size_continuous(breaks = scales::breaks_extended(n = 3)) + 
  guides(color = guide_colorbar(title.position = "top",
                                barwidth = 8,
                                barheight = 0.5,
                                reverse = TRUE,
                                theme = theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))))



p_mf + p_bp + p_cc +  plot_spacer() +
  # plot_layout(guides = 'collect') +
  #  plot_annotation(theme = theme(legend.position = "bottom")) +
  plot_annotation(tag_levels = 'a') &               
  theme(plot.tag = element_text(face = "bold", size = 14)) 