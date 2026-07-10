# p21 proteome 2026 -----------------------------------------------------------

proteome_complete_raw <- read.xlsx("/data/glycobrain_docs/250922_gBrain_Aging_p21_90min_prot_HybridDIA_notNorm_Report.xlsx")


proteome_complete_raw_clean <- proteome_complete_raw %>%
  #dplyr::select(-c(PG.ProteinDescriptions, PG.ProteinNames, PG.CellularComponent, PG.BiologicalProcess, PG.MolecularFunction)) %>%
  dplyr::select(-c(PG.ProteinDescriptions, PG.ProteinNames)) %>%
  dplyr::select(!contains(c("Entire","RestBrain"))) 

colnames(proteome_complete_raw_clean) = str_remove_all(colnames(proteome_complete_raw_clean), "\\[|\\]")
#colnames(proteome_complete_raw_clean) <- gsub("MidBrain", "Midbrain", colnames(proteome_complete_raw_clean))
#colnames(proteome_complete_raw_clean) <- gsub("MidB", "Midb", colnames(proteome_complete_raw_clean))

prot_metadata <- data.frame(colname_raw = names(proteome_complete_raw_clean)[grep("htrms.PG.Quantity", colnames(proteome_complete_raw_clean))])
#prot_metadata$condition <- str_extract(prot_metadata$colname_raw, "(?<=Comp_v5b_\\(300621\\)_).*(?=_50%_try_prot)")
prot_metadata$condition <- str_extract(prot_metadata$colname_raw, "(?<=50k_DIA_).*(?=_P21_WT)")
prot_metadata$replicate <- str_extract(prot_metadata$colname_raw, "(?<=_WT)\\d+")
#prot_metadata$replicate <-  as.numeric(str_extract(prot_metadata$colname_raw, "(?<=#)\\d+"))


region_replacement_pairs <- c("Cor" = "Cortex", 
                              "OB" = "Olfactory", 
                              "Cer" = "Cerebellum",
                              "Tha" = "Thalamus",
                              "Hip" = "Hippocampus",
                              "MidB" = "Midbrain",
                              "Med" = "Medulla",
                              "Hyp" = "Hypothalamus")

prot_metadata <- prot_metadata %>% 
  mutate(across("condition", ~ str_replace_all(., region_replacement_pairs)))


prot_metadata$label <- paste(prot_metadata$condition, prot_metadata$replicate, sep = "_")

prot_sample_columns <- grep("htrms.PG.Quantity", colnames(proteome_complete_raw_clean))

proteome_complete_raw_clean[proteome_complete_raw_clean == "NaN"] <- NA
proteome_complete_raw_clean[prot_sample_columns] <- sapply(proteome_complete_raw_clean[prot_sample_columns],as.numeric)

colnames(proteome_complete_raw_clean)[prot_sample_columns] <- prot_metadata$label

proteome_complete_raw_clean_unique <- DEP2::make_unique(proteome_complete_raw_clean, "PG.Genes", "PG.ProteinGroups", delim = ";")

data_se_prot <- DEP2::make_se(proteome_complete_raw_clean_unique, 
                              columns = prot_sample_columns, 
                              expdesign = prot_metadata)

#df_data_se_prot <- get_df_wide(data_se_prot)
#write.xlsx(df_data_se_prot,"/data/p21_proteome_raw_050726.xlsx")

prot_filter_pg <- DEP2::filter_se(data_se_prot, thr = 0)

#df_prot_filter_pg <- get_df_wide(prot_filter_pg)
#write.xlsx(df_prot_filter_pg,"/data/p21_proteome_filtered_050726.xlsx")

prot_vsn_pg <- normalize_vsn(prot_filter_pg)
set.seed(4)
prot_vsn_pg_RF <- DEP2::impute(prot_vsn_pg, ntree = 200, mtry = 10, fun = "RF", verbose = TRUE)

saveRDS(prot_vsn_pg_RF, 
        "/data/glycobrain_docs/prot_p21_vsn_pg_RF_230426.rds")

