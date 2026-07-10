# regulation state heatmap ------------------------------------------------


stat_res_p21_glycopep_demult_class <- readRDS("/data/glycobrain_docs/stat_res_p21_glycopep4008_demult_class_230426.rds")
df_stat_res_p21_proteome_all_curve_demult <- readRDS("/data/glycobrain_docs/df_stat_res_p21_proteome_all_curve_demult_230426.rds")


glyco_prot_overlap_accessions <- intersect(stat_res_p21_glycopep_demult_class$PG.ProteinAccessions, 
                                           df_stat_res_p21_proteome_all_curve_demult$PG.ProteinGroups)


df_classification_state <- data.frame(PG.ProteinAccessions = glyco_prot_overlap_accessions)
df_classification_state_independ_reg <- df_classification_state
df_classification_state_prot_selective_reg <- df_classification_state
df_classification_state_glyco_selective_reg <- df_classification_state
df_classification_state_co_reg <- df_classification_state
df_classification_state_opposite_reg <- df_classification_state
df_classification_state_not_reg <- df_classification_state
df_classification_state_prot_reg <- df_classification_state

df_result <- data.frame()

regions <- c("Cerebellum", "Medulla", "Pons", "Midbrain", "Hypothalamus", "Thalamus", "Hippocampus", "Cortex", "Olfactory")
all_unique_region_pairs <- combn(regions, 2, FUN = paste, collapse =  "_vs_")

for(comparison in all_unique_region_pairs) {
  # comparison <- opposite_reg_plot_list[current_row,2]
  # gene_of_interest <- opposite_reg_plot_list[current_row,1]
  #  comparison <- "Pons_vs_Hippocampus"
  #  gene_of_interest <- "Igfbp5"
  #print(paste(current_row, comparison, gene_of_interest))
  # print(paste(comparison, gene_of_interest))
  
  df_volcano_glyco <- data.frame(log2FC = stat_res_p21_glycopep_demult_class[[paste0(comparison, "_diff")]],
                                 padj = stat_res_p21_glycopep_demult_class[[paste0(comparison, "_p.adj")]],
                                 epitope = stat_res_p21_glycopep_demult_class$epitope,
                                 name = stat_res_p21_glycopep_demult_class$name,
                                 PG.Genes =  stat_res_p21_glycopep_demult_class$PG.Genes,
                                 PEP.PeptidePosition = stat_res_p21_glycopep_demult_class$PEP.PeptidePosition,
                                 PEP.StrippedSequence = stat_res_p21_glycopep_demult_class$PEP.StrippedSequence,
                                 negLog10Padj = -log10(stat_res_p21_glycopep_demult_class[[paste0(comparison, "_p.adj")]]),
                                 PG.ProteinAccessions = stat_res_p21_glycopep_demult_class$PG.ProteinAccessions)
  
  sigma <- sd(df_volcano_glyco$log2FC)
  x0.fold <- 1     # user-chosen multiplier
  c <- 2           # curvature parameter
  x0 <- x0.fold * sigma
  # Define curve function
  cutoff_fun_glyco <- function(x, c, x0) {c / (abs(x) - x0)}
  
  df_volcano_glyco$significant <- with(df_volcano_glyco, 
                                       ifelse((abs(log2FC) > x0) & (negLog10Padj > cutoff_fun_glyco(log2FC, c, x0)), "TRUE", "FALSE"))
  
  glycopeptides_min10pc <- stat_res_p21_glycopep_demult_class %>%
    dplyr::select(1:38, 42, 43) %>%
    dplyr::select(contains(str_split_i(comparison, "_vs_", 1)), contains(str_split_i(comparison, "_vs_", 2)),
                  name, PG.ProteinAccessions,PG.Genes,Glycan,epitope) %>%
    # filter(epitope %in% c("1T","2T","3T","4T","1Tn")) %>%
    mutate(mean = apply(dplyr::select(., where(is.numeric)), 1, function(x) mean(2^x))) %>%
    ungroup() %>%
    group_by(PG.ProteinAccessions) %>%
    mutate(mean_10pc = .10*max(mean)) %>%
    filter(mean >= mean_10pc) %>%
    ungroup() %>%
    dplyr::select(name) %>%
    as.list()
  
  df_volcano_prot <- data.frame(log2FC_prot = df_stat_res_p21_proteome_all_curve_demult[[paste0(comparison, "_diff")]],
                                padj_prot = df_stat_res_p21_proteome_all_curve_demult[[paste0(comparison, "_p.adj")]],
                                PG.Genes =  df_stat_res_p21_proteome_all_curve_demult$PG.Genes,
                                negLog10Padj_prot = -log10(df_stat_res_p21_proteome_all_curve_demult[[paste0(comparison, "_p.adj")]]),
                                PG.ProteinAccessions = df_stat_res_p21_proteome_all_curve_demult$PG.ProteinGroups)
  
  sigma_prot <- sd(df_volcano_prot$log2FC)
  x0.fold_prot <- 1     # user-chosen multiplier
  c_prot <- 2           # curvature parameter
  x0_prot <- x0.fold_prot * sigma_prot
  # Define curve function
  cutoff_fun_prot <- function(x, c_prot, x0) {c_prot / (abs(x) - x0_prot)}
  
  df_volcano_prot$significant_prot <- with(df_volcano_prot, 
                                           ifelse((abs(log2FC_prot) > x0_prot) & (negLog10Padj_prot > cutoff_fun_prot(log2FC_prot, c_prot, x0_prot)), "TRUE", "FALSE"))
  
  df_volcano_glyco <- df_volcano_glyco %>%
    #  filter(epitope %in% c("1T","2T","3T","4T","1Tn")) %>%
    # filter(epitope %in% c("1T","2T","3T","4T")) %>%
    #filter(epitope %in% c("1Tn")) %>%
    filter(PG.ProteinAccessions %in% glyco_prot_overlap_accessions) %>%
    filter(name %in% glycopeptides_min10pc$name)
  
  df_volcano_prot <- df_volcano_prot %>%
    filter(PG.ProteinAccessions %in% glyco_prot_overlap_accessions)
  
  df_volcano_glyco_prot_comb <- merge(df_volcano_glyco, df_volcano_prot, by.x ="PG.ProteinAccessions", by.y="PG.ProteinAccessions")
  
  df_volcano_glyco_prot_comb <- df_volcano_glyco_prot_comb %>%
    mutate(significant_label_prot=case_when(significant_prot == "TRUE" & log2FC_prot > 0 ~ "up",
                                            significant_prot == "TRUE" & log2FC_prot < 0 ~ "down",
                                            TRUE ~ "not")) %>%
    mutate(significant_label_glyco=case_when(significant == "TRUE" & log2FC > 0 ~ "up",
                                             significant == "TRUE" & log2FC < 0 ~ "down",
                                             TRUE ~ "not"))
  
  list_classified_genes_independ_reg  <- df_volcano_glyco_prot_comb %>%
    filter(significant != significant_prot) %>%
    #filter(significant == "FALSE" & significant_prot == "TRUE") %>%
    # dplyr::rename(PG.Genes = PG.Genes.x) %>%
    dplyr::select(PG.ProteinAccessions) %>%
    distinct() %>%
    mutate(state = 1) %>%
    dplyr::rename("{comparison}" := state)
  
  df_classification_state_independ_reg <- merge(df_classification_state_independ_reg, list_classified_genes_independ_reg, 
                                                by.x ="PG.ProteinAccessions", by.y="PG.ProteinAccessions", all.x = TRUE)
  
  list_classified_genes_glyco_selective_reg  <- df_volcano_glyco_prot_comb %>%
    filter(significant == "TRUE" & significant_prot == "FALSE") %>%
    #filter(significant == "FALSE" & significant_prot == "TRUE") %>%
    # dplyr::rename(PG.Genes = PG.Genes.x) %>%
    dplyr::select(PG.ProteinAccessions) %>%
    distinct() %>%
    mutate(state = 3) %>%
    dplyr::rename("{comparison}" := state)
  
  df_classification_state_glyco_selective_reg <- merge(df_classification_state_glyco_selective_reg, list_classified_genes_glyco_selective_reg, 
                                                       by.x ="PG.ProteinAccessions", by.y="PG.ProteinAccessions", all.x = TRUE)
  
  
  list_classified_genes_prot_selective_reg  <- df_volcano_glyco_prot_comb %>%
    # filter(significant == "TRUE" & significant_prot == "FALSE") %>%
    filter(significant == "FALSE" & significant_prot == "TRUE") %>%
    # dplyr::rename(PG.Genes = PG.Genes.x) %>%
    dplyr::select(PG.ProteinAccessions) %>%
    distinct() %>%
    mutate(state = 5) %>%
    dplyr::rename("{comparison}" := state)
  
  df_classification_state_prot_selective_reg <- merge(df_classification_state_prot_selective_reg, list_classified_genes_prot_selective_reg, 
                                                      by.x ="PG.ProteinAccessions", by.y="PG.ProteinAccessions", all.x = TRUE)
  
  
  
  
  
  list_classified_genes_co_reg  <- df_volcano_glyco_prot_comb %>%
    filter(significant_prot == "TRUE" & significant_label_prot == significant_label_glyco) %>%
    #  dplyr::rename(PG.Genes = PG.Genes.x) %>%
    dplyr::select(PG.ProteinAccessions) %>%
    distinct() %>%
    mutate(state = 1) %>%
    dplyr::rename("{comparison}" := state)
  
  df_classification_state_co_reg <- merge(df_classification_state_co_reg, list_classified_genes_co_reg, 
                                          by.x ="PG.ProteinAccessions", by.y="PG.ProteinAccessions", all.x = TRUE)
  
  list_classified_genes_opposite_reg  <- df_volcano_glyco_prot_comb %>%
    filter(significant_prot == "TRUE" & significant == "TRUE" & significant_label_prot != significant_label_glyco) %>%
    # dplyr::rename(PG.Genes = PG.Genes.x) %>%
    dplyr::select(PG.ProteinAccessions) %>%
    distinct() %>%
    mutate(state = 1) %>%
    dplyr::rename("{comparison}" := state)
  
  df_classification_state_opposite_reg <- merge(df_classification_state_opposite_reg, list_classified_genes_opposite_reg, 
                                                by.x ="PG.ProteinAccessions", by.y="PG.ProteinAccessions", all.x = TRUE)
  
  list_classified_genes_not_reg  <- df_volcano_glyco_prot_comb %>%
    filter(significant_prot == "FALSE" & significant == "FALSE") %>%
    #   dplyr::rename(PG.Genes = PG.Genes.x) %>%
    dplyr::select(PG.ProteinAccessions) %>%
    distinct() %>%
    mutate(state = 1) %>%
    dplyr::rename("{comparison}" := state)
  
  df_classification_state_not_reg <- merge(df_classification_state_not_reg, list_classified_genes_not_reg, 
                                           by.x ="PG.ProteinAccessions", by.y="PG.ProteinAccessions", all.x = TRUE)
  
  list_classified_genes_prot_reg  <- df_volcano_glyco_prot_comb %>%
    filter(significant_prot == "TRUE") %>%
    #  dplyr::rename(PG.Genes = PG.Genes.x) %>%
    dplyr::select(PG.ProteinAccessions) %>%
    distinct() %>%
    mutate(state = 1) %>%
    dplyr::rename("{comparison}" := state)
  
  df_classification_state_prot_reg <- merge(df_classification_state_prot_reg, list_classified_genes_prot_reg, 
                                            by.x ="PG.ProteinAccessions", by.y="PG.ProteinAccessions", all.x = TRUE)
  
  n_coreg <- df_volcano_glyco_prot_comb %>%
    filter(significant_prot == "TRUE" & significant_label_prot == significant_label_glyco) %>%
    dplyr::select(PG.ProteinAccessions) %>%
    distinct() %>% nrow()
  
  n_independ_reg <- df_volcano_glyco_prot_comb %>%
    filter(significant_label_prot != significant_label_glyco) %>%
    dplyr::select(PG.ProteinAccessions) %>%
    distinct() %>% nrow()
  
  n_prot_reg <- df_volcano_glyco_prot_comb %>%
    filter(significant_prot == TRUE) %>%
    dplyr::select(PG.ProteinAccessions) %>%
    distinct() %>% nrow()
  
  current_row <- data.frame(comparison = comparison,
                            n_coreg = n_coreg,
                            n_independ_reg = n_independ_reg,
                            n_prot_reg = n_prot_reg)
  df_result <- rbind(df_result, current_row)    
  
  
  print(paste(comparison))
}




# annotated htmp ----------------------------------------------------------



df_classification_state_independ_reg <- df_classification_state_independ_reg %>% 
  column_to_rownames("PG.ProteinAccessions") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))


#df_classification_independ_reg_rowanno <- df_classification_state_independ_reg %>%
#  column_to_rownames("PG.ProteinAccessions") %>%
#  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
#  mutate(Prot_independ_reg=case_when(rowSums(.) > 0 ~ TRUE,
#                                     TRUE ~ FALSE)) %>%
#  select(c(Prot_independ_reg))


df_classification_state_rowanno <- df_classification_state_independ_reg %>%
  mutate(n_reg = rowSums(.)) %>%
  mutate(n_non_reg = 36 - n_reg) %>%
  select(c(n_reg, n_non_reg))

df_classification_oppositereg_rowanno <- df_classification_state_opposite_reg %>%
  column_to_rownames("PG.ProteinAccessions") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(Prot_opposite_reg=case_when(rowSums(.) > 0 ~ TRUE,
                                     TRUE ~ FALSE)) %>%
  select(c(Prot_opposite_reg))

df_classification_pro_coreg_rowanno <- df_classification_state_co_reg %>%
  column_to_rownames("PG.ProteinAccessions") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(Prot_co_reg=case_when(rowSums(.) > 0 ~ TRUE,
                               TRUE ~ FALSE)) %>%
  select(c(Prot_co_reg))

df_classification_state_nonreg_rowanno <- df_classification_state_not_reg  %>%
  column_to_rownames("PG.ProteinAccessions") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(n_nonreg = rowSums(.)) %>%
  mutate(n_non_nonreg = 36 - n_nonreg) %>%
  select(c(n_nonreg, n_non_nonreg))

df_classification_state_protreg_rowanno <- df_classification_state_prot_reg %>%
  column_to_rownames("PG.ProteinAccessions") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .))) %>%
  mutate(n_prot_reg = rowSums(.)) %>%
  mutate(n_prot_non_reg = 36 - n_prot_reg) %>%
  select(c(n_prot_reg, n_prot_non_reg))


sig_genes_list <- data.frame(first = str_split_i(all_unique_region_pairs, "_vs_", 1),
                             second  = str_split_i(all_unique_region_pairs, "_vs_", 2))

sig_genes_list <- cbind(sig_genes_list, df_result)

sig_genes_list_anno <- sig_genes_list %>%
  column_to_rownames("comparison")

#heatmap_rowanno_genes <- df_stat_res_p21_proteome_all_curve_demult %>%
#  dplyr::filter(PG.ProteinGroups %in% glyco_prot_overlap_accessions) %>%
#  dplyr::select(PG.ProteinGroups, PG.Genes) %>%
# column_to_rownames("PG.ProteinGroups")


#heatmap_rowanno_genes <- df_volcano_glyco_prot_comb %>%
#  dplyr::rename(PG.Genes = PG.Genes.x) %>%
#  dplyr::select(c(PG.ProteinAccessions, PG.Genes)) %>%
#  distinct() %>%
# column_to_rownames("PG.ProteinAccessions") #%>% nrow()

#heatmap_rowanno_genes["Q8C985", "PG.Genes"] <- "Nrxn3b"


metadata_regulation_ht <- data.frame(comparison = all_unique_region_pairs,
                                     Region_1 = str_split_i(all_unique_region_pairs, "_vs_", 1),
                                     Region_2 = str_split_i(all_unique_region_pairs, "_vs_", 2))

metadata_regulation_ht <- metadata_regulation_ht %>%
  column_to_rownames("comparison")



df_classification_state_prot_selective_reg <- df_classification_state_prot_selective_reg %>% 
  column_to_rownames("PG.ProteinAccessions") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))

df_classification_state_glyco_selective_reg <- df_classification_state_glyco_selective_reg %>% 
  column_to_rownames("PG.ProteinAccessions") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))

df_classification_state_selective_reg <- df_classification_state_glyco_selective_reg + df_classification_state_prot_selective_reg

rowanno_prot_selective_reg <- df_classification_state_prot_selective_reg %>%
  mutate(n_prot_reg = rowSums(.)/5) %>%
  dplyr::select(n_prot_reg)

rowanno_glyco_selective_reg <- df_classification_state_glyco_selective_reg %>%
  mutate(n_glyco_reg = rowSums(.)/3) %>%
  dplyr::select(n_glyco_reg)

df_classification_state_rowanno <- cbind(rowanno_glyco_selective_reg, rowanno_prot_selective_reg) %>%
  mutate(n_non_reg = 36 - (n_glyco_reg + n_prot_reg))

heatmap_rowanno_genes <- df_stat_res_p21_proteome_all_curve_demult %>%
  dplyr::filter(PG.ProteinGroups %in% glyco_prot_overlap_accessions) %>%
  dplyr::select(PG.ProteinGroups, PG.Genes) %>%
  column_to_rownames("PG.ProteinGroups")

heatmap_rowanno_genes = heatmap_rowanno_genes[rownames(df_classification_state_selective_reg), ]

my_right_annotation = rowAnnotation(names = anno_text(heatmap_rowanno_genes, gp = gpar(fontsize = 2)),
                                    opposite_reg = df_classification_oppositereg_rowanno$Prot_opposite_reg,
                                    co_reg = df_classification_pro_coreg_rowanno$Prot_co_reg,
                                    gp = gpar(col = "white", lwd = 1),
                                    simple_anno_size = unit(4, "mm"), 
                                    annotation_label = c(NA,"Discordant","Concordant"),
                                    annotation_name_gp = gpar(fontsize = 6),
                                    annotation_name_rot = 45,
                                    annotation_name_side = "bottom",
                                    col = list(opposite_reg =  c("TRUE"="red", "FALSE" = "white"),
                                               co_reg = c("TRUE"="darkgreen", "FALSE" = "white")),
                                    show_legend = FALSE, 
                                    annotation_legend_param = list(opposite_reg = list(title = "Glycopeptides\noppposite\nregulated"),
                                                                   co_reg = list(title = "Glycopeptides\nco-regulated")),
                                    border = FALSE)


lgd1 = Legend(
  title = "Occupancy\nchange", 
  at = c(3, 5, 0), 
  labels = c("Glyco selective", "Protein selective" ,"None"),
  legend_gp = gpar(fill = c("#ec7014", "#662506", "lightgrey")),
  ncol = 1
)

lgd2 = Legend(
  title = "Region", 
  at = c("Olfactory", "Cortex", "Hippocampus", "Thalamus", "Hypothalamus", "Midbrain", "Pons", "Medulla", "Cerebellum"),
  legend_gp = gpar(fill = c("#1F78B4", "#E31A1C", "#666666", "#1B9E77", "#FF7F00", "#A6CEE3", "#6A3D9A", "#FF61C3", "#B2DF8A")),
  ncol = 1
)

lgd3 = Legend(
  title = "Independently\nregulated\nGlycoproteins", 
  col_fun = colorRamp2(breaks = seq(from=0, to=max(sig_genes_list_anno$n_independ_reg, na.rm = TRUE), length.out=10), 
                       colors = c("white", brewer.pal(9, "YlOrBr"))),
  border = "black",
  ncol = 1
)

lgd4 = Legend(
  title = "Regulated\nProteins", 
  col_fun = colorRamp2(breaks = seq(from=0, to=max(sig_genes_list_anno$n_prot_reg, na.rm = TRUE), length.out=10), 
                       colors = c("white", brewer.pal(9, "RdPu"))),
  border = "black",
  ncol = 1
)


ht <- Heatmap(as.matrix(df_classification_state_selective_reg),               
              rect_gp = gpar(col = "white", lwd = 1),               
              col = c("0" = "lightgrey", "3" = "#ec7014", "5" = "#662506"),
              name = "classification",               
              na_col = "white",               
              cluster_rows = TRUE,               
              cluster_columns = TRUE,               
              clustering_distance_columns = 'euclidean',               
              clustering_distance_rows = 'euclidean',               
              clustering_method_columns = "ward.D2",               
              clustering_method_rows = "ward.D2",              
              split = 4,            
              show_row_names = FALSE,               
              row_title_rot = 0,               
              row_dend_side = "left",               
              column_dend_side = "top",               
              column_names_gp = gpar(fontsize = 3),               
              column_names_side = "bottom",               
              column_names_rot = 45,               
              column_title = "Overview of glycoprotein paired regulation states",               
              
              top_annotation = HeatmapAnnotation(                 
                Comparison = cbind(metadata_regulation_ht$Region_2, metadata_regulation_ht$Region_1),                 
                "Independently\nregulated Glycoproteins" = sig_genes_list_anno$n_independ_reg,                 
                "Regulated\nProteins" = sig_genes_list_anno$n_prot_reg,                  
                col = list(Comparison = c("Cortex"="#E31A1C", "Olfactory"="#1F78B4", "Cerebellum"="#B2DF8A",                                            
                                          "Thalamus"="#1B9E77", "Hippocampus"="#666666", "Midbrain"="#A6CEE3",                                            
                                          "Medulla"="#FF61C3", "Hypothalamus"="#FF7F00", "Pons"="#6A3D9A"),                            
                           "Independently\nregulated Glycoproteins" = colorRamp2(breaks = seq(from=0, to=max(sig_genes_list_anno$n_independ_reg, na.rm = TRUE),length.out=10),                                                         
                                                                                 colors = c("white", brewer.pal(9, "YlOrBr"))),                            
                           "Regulated\nProteins" = colorRamp2(breaks = seq(from=0, to=max(sig_genes_list_anno$n_prot_reg, na.rm = TRUE),length.out=10),                                                     
                                                              colors = c("white", brewer.pal(9, "RdPu")))),                 
                show_legend = FALSE, 
                annotation_name_side = "left",                  
                annotation_name_gp = gpar(fontsize = 6),                 
                border = TRUE),                
              right_annotation = my_right_annotation,                
              show_heatmap_legend = FALSE) + 
  
  rowAnnotation("Independent\nregulation" = anno_barplot(df_classification_state_rowanno, gp = gpar(fill = c("n_glyco_reg" = "#ec7014", "n_prot_reg" = "#662506", n_non_reg = "lightgrey"),lwd = NA),  width = unit(1.25, "cm")), annotation_name_rot = 0, annotation_name_gp = gpar(fontsize = 5), annotation_name_side = "top") +   
  rowAnnotation("No\nregulation" = anno_barplot(df_classification_state_nonreg_rowanno, gp = gpar(fill = c(n_nonreg = "ivory4", n_non_nonreg = "lightgrey"),lwd = NA),  width = unit(1.25, "cm")), annotation_name_rot = 0, annotation_name_gp = gpar(fontsize = 5), annotation_name_side = "top") +   
  rowAnnotation("Protein\nregulation" = anno_barplot(df_classification_state_protreg_rowanno, gp = gpar(fill = c(n_prot_reg = "#49006a", n_prot_non_reg = "lightgrey"),lwd = NA),  width = unit(1.25, "cm")), annotation_name_rot = 0, annotation_name_gp = gpar(fontsize = 5), annotation_name_side = "top")


combined_single_column_legend = packLegend(lgd2, lgd3, lgd4, lgd1, 
                                           direction = "vertical", 
                                           max_height = unit(150, "cm"))

draw(ht, 
     ht_gap = unit(1, "mm"), 
     heatmap_legend_side = "right", 
     heatmap_legend_list = combined_single_column_legend)




