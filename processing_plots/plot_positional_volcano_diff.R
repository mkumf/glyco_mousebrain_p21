## diff plot 2 -------------------------------------------------------------

#stat_res_p21_glycopep_demult_class <- readRDS("/data/df_stat_res_p21_glycopep_processed_demult_170326.rds")

stat_res_p21_glycopep_demult_class <- readRDS("/data/glycobrain_docs/stat_res_p21_glycopep4008_demult_class_230426.rds")

df_diff_se_p21_proteome_all_curve <- readRDS("/data/glycobrain_docs/df_stat_res_p21_proteome_all_curve_demult_230426.rds")



#Cx3cl1
#Vgf
#Igfbp5
#Ptprz1

#stat_res_p21_glycopep_demult_class <- stat_res_p21_glycopep  %>%
#  separate_longer_delim(c("PG.ProteinAccessions", "PG.Genes", "PEP.PeptidePosition"), delim=";") %>%
#  separate_longer_delim(c("PEP.PeptidePosition"), delim=",") %>%
#  mutate(PEP.PeptidePosition = as.numeric(PEP.PeptidePosition))

stat_res_p21_glycopep_demult_class_onlycanonical <- stat_res_p21_glycopep_demult_class %>%
  dplyr::filter(epitope %in% c("1T", "2T", "3T", "4T", "1Tn"))

all_overlap <- intersect(stat_res_p21_glycopep_demult_class$PG.Genes, df_diff_se_p21_proteome_all_curve$PG.Genes) 

#c("Cx3cl1", "Vgf", "Igfbp5", "Syt5", "Cntn1", "Ptprz1")


for (gene in all_overlap){
    gene_of_interest <- gene  
 # gene_of_interest <- "Aplp1"
  Accession <- stat_res_p21_glycopep_demult_class[stat_res_p21_glycopep_demult_class$PG.Genes == gene_of_interest, "PG.ProteinAccessions"][1]
  current_accession <- Accession
  
  uniprot_info_df <- rbioapi::rba_uniprot_proteins(accession = Accession)
  gene_name <- as.character(unlist(uniprot_info_df$gene$name[1]))
  long_name <- as.character(unlist(uniprot_info_df$protein$recommendedName[1]))
  prot_data <- drawProteins::feature_to_dataframe(drawProteins::get_features(Accession))
  protein_sequence_full_split <- unlist(strsplit(uniprot_info_df$sequence$sequence[1], split = ""))  
  
  for (comp_loop in all_unique_region_pairs){
    comparison <- "Medulla_vs_Cortex" # Cerebellum vs Olfactory  Medulla vs Cortex
    
    # comparison <- comp_loop
    
    #Accession <- stat_res_p21_glycopep_demult_class[stat_res_p21_glycopep_demult_class$PG.Genes == gene_of_interest, "PG.ProteinAccessions"][1]
    
    #comparison <- loop_df_list[[i, "comparison"]]
    #gene_of_interest <- loop_df_list[[i, "PG.Genes"]]
    #Accession <- loop_df_list[[i, "accessions"]] 
    
    print(paste(comparison, gene_of_interest))
    
    comparison_parts <- strsplit(comparison, "_vs_")[[1]]
    
    df_log2.delta <- stat_res_p21_glycopep_demult_class %>%
      dplyr::select(which(!duplicated(names(.)))) %>%
      dplyr::select(-contains("_vs_")) %>%
      mutate(mean.1 = rowMeans(2^across(contains(comparison_parts[1])), na.rm = TRUE),
             mean.2 = rowMeans(2^across(contains(comparison_parts[2])), na.rm = TRUE)) %>%
      mutate(delta = mean.1 - mean.2) %>%
      mutate(log2.delta = case_when(delta > 0  ~ log2(abs(delta)),
                                    delta < 0  ~ -log2(abs(delta)),
                                    delta == 0 ~ 0,
                                    TRUE ~ NA_real_ )) %>%
      select(log2.delta)
    
    
    
    df_volcano_glyco <- data.frame(log2FC = stat_res_p21_glycopep_demult_class[[paste0(comparison, "_diff")]],
                                   padj = stat_res_p21_glycopep_demult_class[[paste0(comparison, "_p.adj")]],
                                   signifcant = stat_res_p21_glycopep_demult_class[[paste0(comparison, "_significant")]],
                                   Glycan = stat_res_p21_glycopep_demult_class$Glycan,
                                   epitope = stat_res_p21_glycopep_demult_class$epitope,
                                   name = stat_res_p21_glycopep_demult_class$name,
                                   PG.Genes =  stat_res_p21_glycopep_demult_class$PG.Genes,
                                   PEP.PeptidePosition = stat_res_p21_glycopep_demult_class$PEP.PeptidePosition,
                                   PEP.StrippedSequence = stat_res_p21_glycopep_demult_class$PEP.StrippedSequence,
                                   negLog10Padj = -log10(stat_res_p21_glycopep_demult_class[[paste0(comparison, "_p.adj")]]),
                                   PG.ProteinAccessions = stat_res_p21_glycopep_demult_class$PG.ProteinAccessions,
                                   log2.delta = df_log2.delta$log2.delta,
                                   cleavage.classification = stat_res_p21_glycopep_demult_class$cleavage.classification)
    
    sigma <- sd(df_volcano_glyco$log2FC)
    
    x0.fold <- 1     # user-chosen multiplier
    c <- 2           # curvature parameter
    x0 <- x0.fold * sigma
    
    cutoff_fun_glyco <- function(x, c, x0) {c / (abs(x) - x0)}
    
    df_volcano_glyco$signifcant <- with(df_volcano_glyco, 
                                        ifelse((abs(log2FC) > x0) & (negLog10Padj > cutoff_fun_glyco(log2FC, c, x0)), "TRUE", "FALSE"))
    
    #   df_volcano_glyco <- df_volcano_glyco %>%
    #     filter(Glycan %in% c("Hex(1)HexNAc(1)", "Hex(2)HexNAc(2)", "Hex(3)HexNAc(3)", "Hex(4)HexNAc(4)", "HexNAc"))  #%>%
    #dplyr::filter(cleavage.classification == "Full specific") #Semi-specific
    
    #   df_diff_se_p21_proteome_all_curve_demult
    
    df_volcano_protein <- data.frame(log2FC_prot = df_diff_se_p21_proteome_all_curve[[paste0(comparison, "_diff")]],
                                     padj_prot = df_diff_se_p21_proteome_all_curve[[paste0(comparison, "_p.adj")]],
                                     signifcant_prot = df_diff_se_p21_proteome_all_curve[[paste0(comparison, "_significant")]],
                                     name_prot = df_diff_se_p21_proteome_all_curve$name,
                                     negLog10Padj_prot = -log10(df_diff_se_p21_proteome_all_curve[[paste0(comparison, "_p.adj")]]),
                                     PG.Genes_prot = df_diff_se_p21_proteome_all_curve$PG.Genes,
                                     PG.ProteinGroups = df_diff_se_p21_proteome_all_curve$PG.ProteinGroups)
    
    sigma_prot <- sd(df_volcano_protein$log2FC_prot)
    
    x0.fold_prot <- 1     # user-chosen multiplier
    c_prot <- 2           # curvature parameter
    x0_prot <- x0.fold_prot * sigma_prot
    
    cutoff_fun_prot <- function(x, c_prot, x0_prot) {c_prot / (abs(x) - x0_prot)}
    
    df_volcano_protein$signifcant_prot <- with(df_volcano_protein, 
                                               ifelse((abs(log2FC_prot) > x0_prot) & 
                                                        (negLog10Padj_prot > cutoff_fun_prot(log2FC_prot, c_prot, x0_prot)), "TRUE", "FALSE"))
    
    df_volcano_protein_label <- df_volcano_protein %>%
      mutate(significant_label=case_when(
        signifcant_prot == "TRUE" & log2FC_prot > 0 ~ "Prot. up",
        signifcant_prot == "TRUE" & log2FC_prot < 0 ~ "Prot. down",
        TRUE ~ "not sign."))
    
    df_volcano_glyco_prot_comb <- merge(df_volcano_glyco, df_volcano_protein,
                                        by.x = "PG.ProteinAccessions",
                                        by.y = "PG.ProteinGroups",
                                        all.x = TRUE)
    
    #    df_volcano_glyco_prot_comb <- df_volcano_glyco_prot_comb %>%
    #      mutate(Glycan_facet=case_when(
    #        Glycan=="Hex(1)HexNAc(1)" ~ "Hex(1)HexNAc(1)",
    #       Glycan=="Hex(2)HexNAc(2)" ~ "Hex(2)HexNAc(2)",
    #        Glycan=="Hex(3)HexNAc(3)" ~ "Hex(3)HexNAc(3)",
    #        Glycan=="Hex(4)HexNAc(4)" ~ "Hex(4)HexNAc(4)",
    #        Glycan=="HexNAc" ~ "HexNAc",
    # TRUE ~ "HexNAc"
    #        TRUE ~ "Misc. glycans"
    #      ))
    
    df_volcano_glyco_prot_comb[is.na(df_volcano_glyco_prot_comb$signifcant_prot ),]$signifcant_prot <- "missing"
    
    df_volcano_glyco_prot_comb_reorder <- df_volcano_glyco_prot_comb
    
    #  df_volcano_glyco_prot_comb_reorder <- rbind(df_volcano_glyco_prot_comb[df_volcano_glyco_prot_comb$signifcant_prot == "missing" | df_volcano_glyco_prot_comb$signifcant_prot == FALSE, ],
    #                                             df_volcano_glyco_prot_comb[df_volcano_glyco_prot_comb$signifcant_prot == TRUE, ])
    
    df_volcano_glyco_prot_comb_reorder <- df_volcano_glyco_prot_comb_reorder %>%
      mutate(DE_precursor=case_when(
        signifcant_prot == "TRUE" & log2FC_prot > 0 ~ "Prot. up",
        signifcant_prot == "TRUE" & log2FC_prot < 0 ~ "Prot. down",
        signifcant_prot == "FALSE" & log2FC_prot < 0 ~ "not sign.",
        TRUE ~ "missing"))
    
    df_volcano_glyco_prot_comb_reorder <- df_volcano_glyco_prot_comb_reorder %>%
      mutate(pep_significance_flag=case_when(
        signifcant == "TRUE" & log2FC > 0 ~ "Glycopep. up",
        signifcant == "TRUE" & log2FC < 0 ~ "Glycopep. down",
        TRUE ~ "not sign."))
    
    df_volcano_glyco_prot_comb_reorder$gene_pos_label <- paste(df_volcano_glyco_prot_comb_reorder$PG.Genes, df_volcano_glyco_prot_comb_reorder$PEP.PeptidePosition, sep="_")
    df_volcano_glyco_prot_comb_reorder$DE_precursor <- factor(df_volcano_glyco_prot_comb_reorder$DE_precursor , 
                                                              levels=c("Prot. up", "not sign.", "Prot. down", "missing"))
    
    df_volcano_glyco_prot_comb_reorder <- df_volcano_glyco_prot_comb_reorder %>%
      mutate(facet_group = DE_precursor) 
    #%>%
    #    mutate(Glycan_facet = case_when(Glycan %in% c("Hex(1)HexNAc(1)", "Hex(2)HexNAc(2)","Hex(3)HexNAc(3)","Hex(4)HexNAc(4)") ~ "1T-4T",
    #                                TRUE ~ "1Tn"))
    #   
    
    df_volcano_glyco_prot_comb_reorder$sig_group <- interaction(df_volcano_glyco_prot_comb_reorder$DE_precursor,
                                                                df_volcano_glyco_prot_comb_reorder$signifcant, 
                                                                sep = "_")
    
    df_volcano_glyco_prot_comb_reorder <- df_volcano_glyco_prot_comb_reorder %>%
      mutate(GOI=case_when(PG.Genes == gene_of_interest ~ gene_of_interest,
                           TRUE ~ NA))  %>%
      arrange(., Glycan)
    
    
    curve_limit <- (c / max(df_volcano_glyco_prot_comb_reorder$negLog10Padj)) + x0
    curve_limit_prot <- (c / max(df_volcano_protein_label$negLog10Padj)) + x0_prot
    
    gene_of_interest_color <- c("gold")
    names(gene_of_interest_color) <- gene_of_interest
    
    #df_volcano_glyco_prot_comb_reorder %<>% dplyr::filter(cleavage.classification == "Semi-specific") #Full specific Semi-specific
    
    glyco_volcano_plot <- ggplot(data=df_volcano_glyco_prot_comb_reorder[df_volcano_glyco_prot_comb_reorder$PG.Genes %in% gene_of_interest, ], 
                                 aes(x = log2FC, y = negLog10Padj, fill = PEP.PeptidePosition)) +
      geom_point(shape = 21, stroke = .25, size = 4) +
      scale_fill_viridis_c() +
      guides(fill = "none") +
      scale_y_continuous(limits = c(NA, 20), oob = scales::squish) +
      scale_x_continuous(limits = c(-6, 6), oob = scales::squish) +
      # ylim(NA, ceiling(max(df_volcano_glyco_prot_comb_reorder$negLog10Padj)/5)*5) +
      # scale_x_continuous(breaks = pretty(df_volcano_glyco_prot_comb_reorder$log2FC,  n = 7, bounds = FALSE)) +
      stat_function(fun = function(x) cutoff_fun_glyco(x, c, x0), 
                    color="black", linewidth = .25, linetype = "dashed", xlim=c(curve_limit, 6)) +
      stat_function(fun = function(x) cutoff_fun_glyco(x, c, x0), 
                    color="black", linewidth = .25, linetype = "dashed", xlim=c(-6, -(curve_limit))) +
      labs(subtitle = paste0("Glycopeptides: ", gene_of_interest, " - ", Accession, ", " ,comparison_parts[1], " vs ", comparison_parts[2]), 
           x=expression(paste("Fold change, ", log[2])), 
           y= expression(-log[10] ~ p[adj])) +
      #    guides(fill = guide_legend(override.aes = list(size = 4))) +
      theme_bw() +
      #    facet_grid(Glycan_facet ~ facet_group) +
      #   facet_grid( ~ Glycan_facet) +
      theme(plot.subtitle = element_text(size = 8, face="bold"),
            strip.text.x = element_text(size = 8),
            strip.text.y = element_text(size = 8),
            legend.title = element_text(size = 8),
            axis.title = element_text(size = 8))
    
    volcano_prot_plot <- ggplot(data=df_volcano_protein_label[df_volcano_protein_label$PG.Genes_prot %in% gene_of_interest, ],
                                aes(x = log2FC_prot, y = negLog10Padj_prot, fill = significant_label)) +
      geom_point(shape = 21, stroke = .25, size = 4, color = "black", ) +
      scale_fill_manual(values = c("Prot. up"="#b2182b", "not sign." = "lightgrey", "Prot. down"="#2166ac")) +
      #scale_fill_manual(values = gene_of_interest_color,  name = "GOI") +
      guides(fill = "none") +
      ylim(0, 20) +
      xlim(-6, 6) +
      # ylim(NA, ceiling(max(df_volcano_glyco_prot_comb_reorder$negLog10Padj)/5)*5) +
      # scale_x_continuous(breaks = pretty(df_volcano_glyco_prot_comb_reorder$log2FC,  n = 7, bounds = FALSE)) +
      stat_function(fun = function(x) cutoff_fun_prot(x, c_prot, x0_prot), 
                    color="black", linewidth = .25, linetype = "dashed", xlim=c(curve_limit_prot, 6)) +
      stat_function(fun = function(x) cutoff_fun_prot(x, c_prot, x0_prot), 
                    color="black", linewidth = .25, linetype = "dashed", xlim=c(-6, -(curve_limit_prot))) +
      #  labs(subtitle = paste0("Proteome: ", gsub("_vs_", " vs. ", paste0(comparison))), x=expression(paste("Fold change, ", log[2])), y= expression(-log[10] ~ p[adj])) +
      labs(subtitle = paste("Proteome:", gene_of_interest), x=expression(paste("Fold change, ", log[2])), y= expression(-log[10] ~ p[adj])) +
      # guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1.0))) +
      theme_classic() +
      #    facet_grid(Glycan_facet ~ facet_group) +
      #    facet_grid( ~ Glycan_facet) +
      theme(plot.subtitle = element_text(size = 8),
            strip.text.x = element_text(size = 8),
            strip.text.y = element_text(size = 8),
            legend.title = element_text(size = 8),
            axis.title = element_text(size = 8))
    
    multi_volcano_plot <- glyco_volcano_plot + (volcano_prot_plot / plot_spacer()  ) +
      plot_layout(widths = c(3, 1))
    
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
    
    #uniprot_info_df <- rbioapi::rba_uniprot_proteins(accession = Accession)
    #gene_name <- as.character(unlist(uniprot_info_df$gene$name[1]))
    #long_name <- as.character(unlist(uniprot_info_df$protein$recommendedName[1]))
    #prot_data <- drawProteins::feature_to_dataframe(drawProteins::get_features(Accession))
    #protein_sequence_full_split <- unlist(strsplit(uniprot_info_df$sequence$sequence[1], split = ""))
    
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
      
      #   theme(legend.key=element_rect(fill="white")) +
      new_scale("fill") +
      ggpattern::geom_rect_pattern(data = prot_data[prot_data$description == "Disordered",],
                                   aes(xmin=begin, xmax=end, ymin=order-0.05, ymax=order+0.05, fill = description),
                                   color = NA,  pattern_fill = "ivory4", pattern_color = "ivory4", 
                                   alpha = 0, pattern_density = 0.35, pattern = 'stripe', pattern_key_scale_factor=.1) + # pattern_spacing = 0.012
      scale_fill_manual(values = c("black"), name = "Structure") +
      
      new_scale("fill") +
      geom_point(data =  prot_data[prot_data$type == "MOD_RES" | prot_data$type == "CARBOHYD" | prot_data$type == "LIPID",],
                 aes(x = begin,  y = order+0.07, shape=description, fill=description),
                 #             shape = 21,
                 colour = "black",
                 size = 4) +
      scale_shape_manual(name = "PTM", values = c('Phosphorylation'=21, 'N-GlcNAc'=22, 'O-GalNAc'=22, "GPI-anchor"=25)) +
      scale_fill_manual(name = "PTM", values = c('Phosphorylation'="orange", 'N-GlcNAc'="#00A9FF", 'O-GalNAc'="gold", "GPI-anchor"="pink")) +
      
      
      
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
    
    
    p_round <- p + scale_x_continuous(expand = c(0, 0), limits = c(1,length(protein_sequence_full_split))) + 
      theme_minimal() + 
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(colour = "ivory4"),
            panel.grid.minor.x = element_line(colour = "ivory4"),
            panel.background = element_blank(),
            legend.background = element_rect(color = NA),
            axis.title.x = element_text(size = 8),
            plot.margin = margin(t = 0,  # Top margin
                                 r = 0,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = 0))
    
    
    pos_first_glycopep <- min(as.numeric(df_volcano_glyco_prot_comb_reorder[df_volcano_glyco_prot_comb_reorder$PG.Genes == gene_of_interest, "PEP.PeptidePosition"]))
    pos_last_glycopep <- max(as.numeric(df_volcano_glyco_prot_comb_reorder[df_volcano_glyco_prot_comb_reorder$PG.Genes == gene_of_interest, "PEP.PeptidePosition"]))
    
    mydf <- data.frame(id = rep("Position", pos_last_glycopep-pos_first_glycopep+1), pos = seq(pos_first_glycopep, pos_last_glycopep))
    
    position_indication_plot <- ggplot(mydf, aes(x=pos, y=id, color=pos,  fill=pos)) +
      geom_tile() +
      viridis::scale_color_viridis(name = "Position\nGlycopeptide") +
      viridis::scale_fill_viridis() +
      guides(fill = "none") +
      scale_x_discrete(breaks = factor(seq(1,length(protein_sequence_full_split))),
                       limits = factor(seq(1,length(protein_sequence_full_split))),
                       labels = protein_sequence_full_split,
                       expand = c(0, 0)) +
      theme_minimal() + 
      theme(axis.title.y=element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x =  element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.background = element_blank(),
            legend.background = element_rect(color = NA),
            plot.margin = margin(t = 0,  # Top margin
                                 r = 0,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = 0))
    
    plotting_df <- df_volcano_glyco_prot_comb_reorder %>%
      dplyr::filter(PG.Genes == paste0(gene_of_interest)) %>%
      # dplyr::filter(epitope %in% c("1T", "2T", "1Tn")) %>%
      #   mutate(epitope = factor(epitope, levels = c("1T", "2T", "1Tn"))) %>%
      mutate(PEP.PeptidePosition = as.numeric(PEP.PeptidePosition)) %>%
      mutate(num.AA = str_count(PEP.StrippedSequence)) %>%
      mutate(end.position = (num.AA + PEP.PeptidePosition)) %>%
      tidyr::complete(pep_significance_flag, epitope) %>%
      mutate(epitope_facet = "all")
    
    plotting_df_1T <-  plotting_df %>%
      dplyr::filter(epitope %in% c("1T")) %>%
      mutate(epitope_facet = "1T")
    
    plotting_df_1Tn <-  plotting_df %>%
      dplyr::filter(epitope %in% c("1Tn")) %>%
      mutate(epitope_facet = "1Tn")
    
    plotting_df <- rbind(plotting_df, plotting_df_1T, plotting_df_1Tn)
    
    plotting_df %<>%  mutate(epitope_facet = factor(epitope_facet, levels = c("all", "1T", "1Tn"))) %>%
      mutate(pep_significance_flag = factor(pep_significance_flag, levels = c("Glycopep. up", "not sign.", "Glycopep. down")))
    
    
    plotting_df <- plotting_df[order(plotting_df$signifcant), ]
    
    #plotting_df %<>% dplyr::filter(cleavage.classification == "Semi-specific")
    
    # values = c("Glycopep. up"="#b2182b", "not sign." = "lightgrey", "Glycopep. down"="#2166ac"),
    diffplot <- ggplot(plotting_df, aes(xmin=PEP.PeptidePosition, xmax=end.position, ymin=0, ymax=log2.delta, alpha=0.5, fill = pep_significance_flag)) +
      geom_rect(linewidth = 0, show.legend = TRUE) +
      scale_alpha(guide = 'none') +
      scale_fill_manual(name = "Regulation",
                        values = c("Glycopep. up"="#b2182b", "not sign." = "lightgrey", "Glycopep. down"="#2166ac"),
                        breaks = c("Glycopep. up", "not sign.", "Glycopep. down"),
                        labels = c("up", "n.s.", "down"),
                        drop = FALSE) +
      scale_x_discrete(breaks = factor(seq(1,length(protein_sequence_full_split))),
                       limits = factor(seq(1,length(protein_sequence_full_split))),
                       labels = protein_sequence_full_split,
                       expand = c(0, 0)) +
      #scale_y_continuous(limits = c(-max(abs(plotting_df$log2FC), na.rm = TRUE), max(abs(plotting_df$log2FC), na.rm = TRUE))) +
      scale_y_continuous(limits = c(-ceiling(max(abs(plotting_df$log2.delta), na.rm = TRUE)), ceiling(max(abs(plotting_df$log2.delta), na.rm = TRUE)))) +
      geom_hline(yintercept=0) +
      theme_bw() + 
      labs(y= expression(log[2] ~ Delta ~ abundance)) +
      theme(#axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        #     panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(),
        legend.background = element_rect(color = NA),
        axis.title.y = element_text(size = 8),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) +
      facet_wrap(epitope_facet ~ ., strip.position = "right", ncol = 1)
    
    
    prot_schematic_element <- position_indication_plot / p_round  +
      plot_layout(heights = c(1, 3))
    
    (free(multi_volcano_plot, type = "label") / 
        prot_schematic_element / 
        free(diffplot, type = "label")) +
      plot_layout(guides = 'collect',
                  heights = c(2, 1, 2)) +
      plot_annotation(theme = theme(legend.position = "right"),
                      tag_levels = 'a') &               
      theme(plot.tag = element_text(face = "bold", size = 10),
            legend.margin = margin(5, 0, 5, 0),    # Alle Ränder auf 0 setzen
            legend.spacing = unit(0, "cm"),        # Abstand zwischen Legenden minimieren
            # legend.spacing.x = unit(0.2, "cm"),
            #legend.spacing.y = unit(0.2, "cm"),
            legend.key.size = unit(0.25, "cm"))     
    
    
    
    filename_plot <- paste0(gene_of_interest, "_", Accession, "_", comparison, "_glycopep_pos_volc_diffplotv2", ".pdf")
    path <- paste0("/data/", gene_of_interest, "/")  
    
    ggsave(filename_plot,
           plot = final_plot,
           path = path,
           width = 297, height = 210, units = "mm",
           create.dir = TRUE)
  }  } 