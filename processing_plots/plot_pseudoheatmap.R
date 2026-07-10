# pseudoheatmap -----------------------------------------------------------

df_stat_res_p21_proteome_all_curve_demult <- readRDS("/data/glycobrain_docs/df_stat_res_p21_proteome_all_curve_demult_230426.rds")
stat_res_p21_glycopep_demult_class <- readRDS("/data/glycobrain_docs/stat_res_p21_glycopep4008_demult_class_230426.rds")

#se_prot_pep_vsn_RF <- readRDS("/data/se_prot_pep_vsn_RF_230526.rds")
#df_prot_pep_vsn_RF <- get_df_wide(se_prot_pep_vsn_RF)


overlap_glyco_prot_accessions <- intersect(stat_res_p21_glycopep_demult_class$PG.ProteinAccessions, df_prot_pep_vsn_RF$PG.ProteinAccessions)
glyco_no_prot_accessions <- setdiff(stat_res_p21_glycopep_demult_class$PG.ProteinAccessions, df_prot_pep_vsn_RF$PG.ProteinAccessions)


jet.colors <- colorRampPalette(c("black", "#00007F", "blue", "#007FFF", "cyan", "green", "yellow", "#FF7F00", "red", "#7F0000", "purple"))

#length(unique(stat_res_p21_glycopep_demult_class$PG.ProteinAccessions))

#bw.colors <- colorRampPalette(c("grey95", "black"))

for(current_accession in glyco_no_prot_accessions){
  
  #  current_accession  <- "Q3UDR8"
  Accession <- current_accession
  #  print(current_accession)
  uniprot_info_df <- rbioapi::rba_uniprot_proteins(accession = Accession)
  gene_name <- as.character(unlist(uniprot_info_df$gene$name[1]))
  long_name <- as.character(unlist(uniprot_info_df$protein$recommendedName[1]))
  prot_data <- drawProteins::feature_to_dataframe(drawProteins::get_features(Accession))
  protein_sequence_full_split <- unlist(strsplit(uniprot_info_df$sequence$sequence[1], split = ""))
  
  df_Glycopep_POI <- stat_res_p21_glycopep_demult_class[grep(paste0(Accession), stat_res_p21_glycopep_demult_class$PG.ProteinAccessions), ]
  
  Glycopep_region_cols <- unique(sapply(strsplit(grep("_\\d+$", colnames(df_Glycopep_POI), value = TRUE),"_", fixed = T), `[`, 1))
  #  df_region_means <- as.data.frame(sapply(region_cols, function (x) rowMeans(df_Glycopep_imp_pg_GSimp_POI[grep(x, names(df_Glycopep_imp_pg_GSimp_POI))])))
  df_Glycopep_region_means <- Glycopep_region_cols %>%
    # Use imap to apply a function to each element in region_cols and keep track of the names
    imap(~ {
      # Select columns that match the current region_col
      selected_cols <- df_Glycopep_POI %>%
        dplyr::select(matches(paste0("^", .x, "_\\d+$")))
      # Calculate row means
      if (ncol(selected_cols) > 0) {
        rowMeans(2^selected_cols, na.rm = TRUE)
      } else {
        NA
      }
    }) %>%
    # Convert the list to a data frame and set the column names
    bind_cols() %>%
    as.data.frame() %>%
    setNames(Glycopep_region_cols)
  
  plotting_df_Glycopep <- cbind(df_Glycopep_region_means,
                                PG.ProteinAccessions = df_Glycopep_POI$PG.ProteinAccessions,
                                PG.Genes =  df_Glycopep_POI$PG.Genes,
                                PEP.PeptidePosition =  df_Glycopep_POI$PEP.PeptidePosition,
                                PEP.StrippedSequence =  df_Glycopep_POI$PEP.StrippedSequence,
                                Glycan = df_Glycopep_POI$Glycan,
                                epitope = df_Glycopep_POI$epitope)
  
  plotting_df_Glycopep_long <- plotting_df_Glycopep  %>% pivot_longer(cols = contains(Glycopep_region_cols), 
                                                                      names_to = "region", 
                                                                      values_to = "abundance")
  # plotting_df_Glycopep_long <- plotting_df_long[order(plotting_df_long$abundance, decreasing = FALSE),]
  
  test_df_glycopep <- plotting_df_Glycopep_long %>%
    mutate(peptide_length = nchar(PEP.StrippedSequence)) %>%
    uncount(peptide_length, .id = "offset") %>%
    mutate(position = as.numeric(PEP.PeptidePosition) + offset - 1) %>%
    group_by(position, region, Glycan) %>%
    summarize(abundance = sum(abundance), .groups = "drop") %>%
    complete(position = 1:length(protein_sequence_full_split),
             region = unique(plotting_df_Glycopep_long$region),
             Glycan = unique(plotting_df_Glycopep_long$Glycan),
             fill = list(abundance = 0)) %>%
    arrange(position, region, Glycan)
  
  test_df_glycopep_log <- test_df_glycopep %>%
    #  filter(Glycan == "Hex(1)HexNAc(1)") %>%
    mutate(abundance = case_when(
      abundance != 0 ~ log2(abundance),
      TRUE ~ 0))
  
  test_df_glycopep_log$region <- factor(test_df_glycopep_log$region, levels = c("Cerebellum",
                                                                                "Medulla",
                                                                                "Pons",
                                                                                "Midbrain",
                                                                                "Hypothalamus",
                                                                                "Thalamus",
                                                                                "Hippocampus",
                                                                                "Cortex",
                                                                                "Olfactory")) 
  
  
  p21_protpep_abundance_delog_POI <- df_prot_pep_vsn_RF[grep(paste0(Accession), df_prot_pep_vsn_RF$PG.ProteinAccessions),]
  
  region_cols_pep_prot <- unique(sapply(strsplit(grep("_\\d+$", colnames(p21_protpep_abundance_delog_POI), value = TRUE),"_", fixed = T), `[`, 1))
  
  df_region_means_pep_prot <- region_cols_pep_prot %>%
    # Use imap to apply a function to each element in region_cols and keep track of the names
    imap(~ {
      # Select columns that match the current region_col
      selected_cols <- p21_protpep_abundance_delog_POI %>%
        dplyr::select(matches(paste0("^", .x, "_\\d+$"))) 
      # Calculate row means
      if (ncol(selected_cols) > 0) {
        rowMeans(2^selected_cols, na.rm = TRUE)
      } else {
        NA
      }
    }) %>%
    # Convert the list to a data frame and set the column names
    bind_cols() %>%
    as.data.frame() %>%
    setNames(region_cols_pep_prot)
  
  plotting_df_pep_prot <- cbind(df_region_means_pep_prot,
                                PG.ProteinAccessions = p21_protpep_abundance_delog_POI$PG.ProteinAccessions,
                                PG.Genes =  p21_protpep_abundance_delog_POI$PG.Genes,
                                PEP.PeptidePosition =  p21_protpep_abundance_delog_POI$PEP.PeptidePosition,
                                PEP.StrippedSequence = p21_protpep_abundance_delog_POI$PEP.StrippedSequence)
  
  plotting_df_long_pep_prot <- plotting_df_pep_prot  %>% pivot_longer(cols = contains(region_cols_pep_prot), names_to = "region",values_to = "abundance")
  
  test_df_prot <- plotting_df_long_pep_prot %>%
    mutate(peptide_length = nchar(PEP.StrippedSequence)) %>%
    uncount(peptide_length, .id = "offset") %>%
    mutate(position = as.numeric(PEP.PeptidePosition) + offset - 1) %>%
    group_by(position, region) %>%
    summarize(abundance = sum(abundance), .groups = "drop") %>%
    complete(position = 1:length(protein_sequence_full_split),
             region = unique(plotting_df_long_pep_prot$region),
             fill = list(abundance = 0)) %>%
    arrange(position, region)  
  
  test_df_fixed_prot <- test_df_prot %>%
    #  filter(Glycan == "Hex(1)HexNAc(1)") %>%
    mutate(abundance = case_when(
      abundance != 0 ~ log2(abundance),
      TRUE ~ 0))
  
  test_df_fixed_prot$Glycan <- "Peptidome"
  
  test_df_fixed_prot$region <- factor(test_df_fixed_prot$region, levels = c("Cerebellum",
                                                                            "Medulla",
                                                                            "Pons",
                                                                            "Midbrain",
                                                                            "Hypothalamus",
                                                                            "Thalamus",
                                                                            "Hippocampus",
                                                                            "Cortex",
                                                                            "Olfactory")) 
  
  test_df_fixed_prot <- test_df_fixed_prot %>%
    filter(abundance != 0)
  test_df_glycopep_log_nonzero <- test_df_glycopep_log %>%
    filter(abundance != 0)
  
  hline_data <- data.frame(region = factor(rep(c("Cerebellum", "Medulla", "Pons", "Midbrain", "Hypothalamus", "Thalamus", "Hippocampus", "Cortex", "Olfactory"), 13), 
                                           levels = c("Cerebellum", "Medulla", "Pons", "Midbrain", "Hypothalamus", "Thalamus", "Hippocampus", "Cortex", "Olfactory")),
                           Glycan = rep(c("Peptidome", 
                                          "Hex(1)HexNAc(1)", "Hex(2)HexNAc(2)", "Hex(3)HexNAc(3)", "Hex(4)HexNAc(4)", "HexNAc", 
                                          "Hex(5)HexNAc(5)", "HexNAc(2)", "HexNAc(3)", "HexNAc(4)", "Hex(1)HexNAc(2)", "Hex(1)HexNAc(3)", "Hex(2)HexNAc(3)"), each=9))
  hline_data <- hline_data %>%
    filter(Glycan %in% append(unique(test_df_glycopep_log$Glycan), "Peptidome")) %>%
    dplyr::filter(Glycan != "Peptidome")
  
  
  
  ### integrate TMT -----------------------------------------------------------
  
  positional_data_comb_filt <- positional_data_comb %>%
    dplyr::filter(Accession == paste0(current_accession)) %>%
    dplyr::select(-c(sugar, AA))  %>%
    dplyr::rename(accession = Accession,
                  begin = position) %>%
    mutate(end = begin,
           type = "CARBOHYD",
           description = "O-GalNAc",
           length = 0,
           entryName = "_MOUSE",
           taxid = 10090,
           order = 1)
  
  prot_data <- rbind(prot_data, positional_data_comb_filt)
  # jet.colors <- colorRampPalette(c("black", "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000", "purple"))
  
  ## draw ggplot -------------------------------------------------------------
  
  g <- ggplot() +
    geom_tile(data=test_df_glycopep_log, aes(x=factor(position), y=region, fill=abundance, color=abundance), show.legend = FALSE) +
    new_scale("fill") +
    new_scale("color") +
    geom_hline(data=hline_data, aes(yintercept = region), color="grey", linewidth = 0.2, show.legend = FALSE) +
    geom_tile(data=test_df_glycopep_log_nonzero, aes(x=factor(position), y=region, fill=abundance, color=abundance)) +
    #scale_fill_gradientn(colors = jet.colors(11), name = expression(log[2]~abundance),  limits = c(min(df_p21_gylcopep_abundance), max(df_p21_gylcopep_abundance)), breaks = scales::breaks_extended(n = 6)) +
    #scale_color_gradientn(colors = jet.colors(11),  guide = "none",  limits = c(min(df_p21_gylcopep_abundance), max(df_p21_gylcopep_abundance))) +                     
    scale_fill_gradientn(colors = jet.colors(11), name = expression(log[2]~abundance),  limits = c(11.38456, 24.71704), breaks = scales::extended_breaks(n = 6), oob = squish) +
    scale_color_gradientn(colors = jet.colors(11),  guide = "none",  limits = c(11.38456, 24.71704), oob = squish) +                     
    #scale_fill_viridis(option="turbo", name = expression(log[2]~abundance),  limits = c(bottom_z, top_z), breaks = scales::breaks_extended(n = 6)) +
    #scale_color_viridis(option="turbo",  guide = "none",  limits = c(bottom_z, top_z)) +
    new_scale("fill") +
    new_scale("color") +
    #  geom_tile(data = test_df_fixed_prot, aes(x=position, y=region, fill=abundance, color =abundance)) +
    scale_fill_gradientn(colors = bw.colors(2), name = expression(log[2]~abundance),  limits = c(10.3729, 22.61196), breaks = scales::extended_breaks(n = 6), oob = squish) +
    scale_color_gradientn(colors = bw.colors(2),  guide = "none",  limits = c(10.3729, 22.61196), oob = squish) +                     
    facet_wrap(~ factor(Glycan, levels=c("Hex(1)HexNAc(1)", "Hex(2)HexNAc(2)", "Hex(3)HexNAc(3)", "Hex(4)HexNAc(4)", "HexNAc",
                                         "Hex(5)HexNAc(5)", "HexNAc(2)", "HexNAc(3)", "HexNAc(4)", "Hex(1)HexNAc(2)", "Hex(1)HexNAc(3)", "Hex(2)HexNAc(3)")), ncol = 1,  strip.position="right") +
    ggtitle(paste(Accession, gene_name, long_name, sep = " - ")) +
    scale_y_discrete(limits=c("Cerebellum", "Medulla", "Pons", "Midbrain", "Hypothalamus", "Thalamus", "Hippocampus", "Cortex", "Olfactory")) +
    scale_x_discrete(breaks = factor(seq(1,length(protein_sequence_full_split))),
                     limits = factor(seq(1,length(protein_sequence_full_split))),
                     labels = protein_sequence_full_split,
                     expand = c(0, 0)) +
    theme_bw() +
    theme(axis.text.x = element_text(size=2),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          strip.text.y = element_text(color = "black"),
          axis.ticks.x = element_line(linewidth = 0.2),
          panel.grid = element_blank())
  
  prot_data[prot_data$type == "TRANSMEM","description"] <- "Transmembrane"
  prot_data[prot_data$description == "N-linked (GlcNAc...) asparagine","description"] <- "N-GlcNAc"
  prot_data[prot_data$description == "Phosphoserine","description"] <- "Phosphorylation"
  prot_data[prot_data$description == "Phosphothreonine","description"] <- "Phosphorylation"
  prot_data[prot_data$description == "GPI-anchor amidated serine","description"] <- "GPI-anchor"
  
  p <- ggplot(prot_data) +
    ylim(0.8, 1.15) +
    labs(x = "Amino acid number", y = "") +
    geom_rounded_rect(data = prot_data[prot_data$type == "SIGNAL" | prot_data$type == "PROPEP",],
                      mapping=ggplot2::aes(xmin=begin, xmax=end, ymin=order-0.05, ymax=order+0.05),
                      colour = "grey", fill = "grey", linewidth = 0.5,  alpha = 1.0, radius=.1) +
    geom_rounded_rect(data = prot_data[prot_data$type == "CHAIN",],
                      mapping=ggplot2::aes(xmin=begin, xmax=end, ymin=order-0.05, ymax=order+0.05),
                      colour = "black", fill = "grey", linewidth = 0.5,alpha = 1.0, radius=.1) +
    geom_rect(data = prot_data[prot_data$type == "TOPO_DOM" | prot_data$type == "TRANSMEM",],
              mapping=ggplot2::aes(xmin=begin, xmax=end, ymin=order-0.20, ymax=order-0.15, fill=description),
              alpha = 0.75, linetype=0) + 
    scale_fill_manual(name = "Localisation",
                      breaks=c("Cytoplasmic", "Transmembrane", "Extracellular", "Intragranular", "Vesicular", "Lumenal"),
                      values = c("Cytoplasmic" = "#F8766D","Extracellular" = "#00B8E7", "Transmembrane" = "#7CAE00", "Intragranular" = "#CD0BBC", "Vesicular" = "#CD0BBC", "Lumenal"="#CD0BBC")) +
    new_scale("fill") +
    geom_rounded_rect(data = prot_data[prot_data$type == "DOMAIN",],
                      mapping=ggplot2::aes(xmin=begin, xmax=end, ymin=order-0.1, ymax=order+0.10, fill=description, radius=.1),
                      alpha = 1.0, color="black") + 
    ggthemes::scale_fill_tableau("Tableau 20", name = "Domains") +
    new_scale("fill") +
    geom_point(data =  prot_data[prot_data$type == "MOD_RES" | prot_data$type == "CARBOHYD" | prot_data$type == "LIPID",],
               aes(x = begin,  y = order+0.07, shape=description, fill=description),
               #             shape = 21,
               colour = "black",
               size = 4) +
    scale_shape_manual(name = "PTM", values = c('Phosphorylation'=21, 'N-GlcNAc'=22, 'O-GalNAc'=22, "GPI-anchor"=25)) +
    scale_fill_manual(name = "PTM", values = c('Phosphorylation'="orange", 'N-GlcNAc'="#00A9FF", 'O-GalNAc'="gold", "GPI-anchor"="pink")) +
    #   theme(legend.key=element_rect(fill="white")) +
    new_scale("fill") +
    ggpattern::geom_rect_pattern(data = prot_data[prot_data$description == "Disordered",],
                                 aes(xmin=begin, xmax=end, ymin=order-0.05, ymax=order+0.05, fill = description),
                                 color = NA,  pattern_fill = "ivory4", pattern_color = "ivory4", 
                                 alpha = 0, pattern_density = 0.35, pattern = 'stripe', pattern_key_scale_factor=.1) + # pattern_spacing = 0.012
    scale_fill_manual(values = c("black"), name = "Structure") +
    theme(legend.key.size = unit(0.5, 'cm')) +
    scale_x_continuous(expand = c(0, 0), limits = c(1,length(protein_sequence_full_split))) + 
    theme_minimal() + 
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(colour = "ivory4"),
          panel.grid.minor.x = element_line(colour = "ivory4"),
          panel.background = element_blank(),
          legend.background = element_rect(color = NA),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0))
  
  
  p <- p + scale_x_continuous(expand = c(0, 0), limits = c(1,length(protein_sequence_full_split))) + 
    theme_minimal() + 
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(colour = "ivory4"),
          panel.grid.minor.x = element_line(colour = "ivory4"),
          panel.background = element_blank(),
          legend.background = element_rect(color = NA),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0))
  
  #plot_grid(g, p, ncol = 1, align = "v", axis = "lr", rel_heights = c(length(unique(plotting_df_long_Glycopep_AA_centric_select$Glycan))+2, 1.5))
  g / p + plot_layout(ncol = 1, guides = "collect", 
                      heights = unit(c(length(unique(test_df_glycopep_log$Glycan))+1, 1), c('null'))) 
  
  
  protein_width <- (length(protein_sequence_full_split)*0.75) + 50
  plot_height <- (length(unique(test_df_glycopep_log$Glycan))*40) + 80
  
  filename_plot <- paste0(Accession, "_", gene_name,"_glycopep_pseudoheatmap", ".pdf")
  #  ggsave(filename_plot,
  #        path = "/data/pseudoheatmaps_070626/",
  #         width = protein_width, height = plot_height, units = "mm",limitsize = FALSE,
  #         create.dir = TRUE)
  
  #  test_df_glycopep_log$AA <- rep(protein_sequence_full_split, each = 9*length(unique(test_df_glycopep_log$Glycan)))
  
  #  test_df_glycopep_log <- test_df_glycopep_log %>% filter(abundance != 0)
  
  #  filename_data <- paste0("/Users/mkummerfeld/Downloads/web_data/non_zero/", Accession, "_", gene_name,"_data_pseudoheatmap", ".csv")
  #  write.csv(test_df_glycopep_log, filename_data, row.names = FALSE)
  
  #  prot_data_export <- prot_data %>%
  #   dplyr::select(-c(accession, entryName, taxid, order))
  
  #  path_prot_data_export <- "/data/glycoprot_anno_080626/"
  #  filename_prot_data_export <- paste0(path_prot_data_export, Accession, "_", gene_name,"_annotation_data", ".csv")
  
  #   if(!dir.exists(path_prot_data_export)){ 
  #   dir.create(path_prot_data_export)
  #  }
  
  
  #  write.csv(prot_data_export, filename_prot_data_export, row.names = FALSE, quote = FALSE)
  
  pseudoheat_data_export <- test_df_glycopep_log %>% 
    dplyr::filter(abundance != 0) %>%
    mutate(abundance = round(abundance, 2))
  
  path_pseudoheat_data_export <- "/data/pseudoheat_data_080626/"
  
  filename_pseudoheat_data_export <- paste0(path_pseudoheat_data_export, Accession, "_", gene_name,"_data_pseudoheatmap", ".csv")
  
  
  if(!dir.exists(path_pseudoheat_data_export)){ 
    dir.create(path_pseudoheat_data_export)
  }
  
  write.csv(pseudoheat_data_export, filename_pseudoheat_data_export, row.names = FALSE,  quote = FALSE)
}
