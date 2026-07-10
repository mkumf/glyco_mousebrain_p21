## load data ---------------------------------------------------------------

glycopep_raw <- read.xlsx("/data/aging_p21/FINAL/2024_Lumos_mBrain_p21_Aging_LB2_NotNorm_Report.xlsx")

glycopep_raw[glycopep_raw$PG.ProteinAccession == "Q9D8N1", ]$PG.Genes <- "C11orf24"
glycopep_raw[glycopep_raw$PG.ProteinAccession == "P0DI97", ]$PG.ProteinAccessions <- "Q9CS84"
glycopep_raw[glycopep_raw$PG.ProteinAccession == "Q9Z0F1", ]$PG.ProteinAccessions <- "Q6R0H7"

glycopep_raw_clean <- glycopep_raw %>%
  dplyr::select(-c("PG.ProteinDescriptions", "[2].240229_Lumos_nLC1200_gBrain_Aging_glyco_(250124)_1h_C18_3u_40w_30k_DIA_Cor_P21_WT3_R1.raw.PEP.MS2Quantity"))

colnames(glycopep_raw_clean) = str_remove_all(colnames(glycopep_raw_clean), "\\[|\\]")

glycopep_raw_clean_nonox <- glycopep_raw_clean[!grepl("Oxidation", glycopep_raw_clean$EG.ModifiedSequence), ]

#glycopep_raw_clean_nonox <- glycopep_raw_clean

glycopep_raw_epitope <- glycopep_raw_clean_nonox %>%
  mutate(Glycan = sub(".*?\\[([^ ]+).*", "\\1", .$EG.ModifiedSequence)) %>%
  mutate(epitope=case_when(Glycan=="Hex(1)HexNAc(1)" ~ "1T",
                           Glycan=="Hex(2)HexNAc(2)" ~ "2T",
                           Glycan=="Hex(3)HexNAc(3)" ~ "3T",
                           Glycan=="Hex(4)HexNAc(4)" ~ "4T",
                           Glycan=="Hex(5)HexNAc(5)" ~ "5T",
                           Glycan=="HexNAc" ~ "1Tn",
                           Glycan=="HexNAc(2)" ~ "2Tn",
                           Glycan=="HexNAc(3)" ~ "3Tn",
                           Glycan=="HexNAc(4)" ~ "4Tn",
                           Glycan=="HexNAc(5)" ~ "5Tn",
                           Glycan=="Hex(1)HexNAc(2)" ~ "1T_1Tn",
                           Glycan=="Hex(1)HexNAc(3)" ~ "1T_2Tn",
                           Glycan=="Hex(2)HexNAc(3)" ~ "2T_1Tn",
                           TRUE ~ "misc.")) %>%
  relocate(Glycan, .after = EG.ModifiedSequence)  %>%
  relocate(epitope, .after = Glycan)

### QScore filt -------------------------------------------------------------

glycopep_qscore_filt <- read.xlsx("/data/glycobrain_docs/df_p21_glycopep_processed_demult_190126_QualityScoreFilt.xlsx")

glycopep_qscore_pass_list <- glycopep_qscore_filt %>%
  dplyr::filter(Quality.Score.Filter == "+")  %>%
  mutate(glycopeptide.ID = paste0(PEP.StrippedSequence, "_", epitope))  %>%
  dplyr::select(glycopeptide.ID) %>%
  unlist()

glycopep_qscore_pass <- glycopep_raw_epitope %>%
  mutate(glycopeptide.ID = paste0(PEP.StrippedSequence, "_", epitope)) %>%
  dplyr::filter(glycopeptide.ID %in% glycopep_qscore_pass_list) %>%
  dplyr::select(-glycopeptide.ID)

glycopep_filt <- glycopep_qscore_pass %>%
  # mutate(n.glycans = str_extract(epitope, "\\d+$")) %>%
  mutate(sum.glycans = str_extract_all(epitope, "\\d") %>%
           lapply(as.numeric) %>%
           sapply(sum)) %>%
  mutate(sum.acceptors = str_count(PEP.StrippedSequence, "[STY]")) %>%
  dplyr::filter(sum.glycans <= sum.acceptors)


#glycopep_filt <- glycopep_raw_epitope

# DEP2 --------------------------------------------------------------------

glycopep_filt_unique <- DEP2::make_unique(glycopep_filt, "PG.Genes", "PG.ProteinAccessions", delim = ";")

glycopep_metadata <- data.frame(colname_raw = names(glycopep_filt)[grep("raw.PEP.MS2Quantity", colnames(glycopep_filt))])
glycopep_metadata$condition <- str_remove(glycopep_metadata$colname_raw, ".*?_DIA_")
glycopep_metadata$condition <- gsub("\\_.*", "", glycopep_metadata$condition)
glycopep_metadata$replicate <- str_remove(glycopep_metadata$colname_raw, ".*?P21_WT")
glycopep_metadata$replicate <- gsub("\\_.*", "", glycopep_metadata$replicate)
#glycopep_metadata$colname_raw <- str_remove(glycopep_metadata$colname_raw, ".*?\\ ")
glycopep_metadata$timepoint <- c(rep("March", 19), rep("Jan", 16))
glycopep_metadata$Batch <- c(rep(1, 19), rep(2, 16))

region_replacement_pairs <- c("Cor" = "Cortex", 
                              "OB" = "Olfactory", 
                              "Cer" = "Cerebellum",
                              "Tha" = "Thalamus",
                              "Hip" = "Hippocampus",
                              "MidB" = "Midbrain",
                              "Med" = "Medulla",
                              "Hyp" = "Hypothalamus")

glycopep_metadata <- glycopep_metadata %>% 
  mutate(across("condition", ~ str_replace_all(., region_replacement_pairs)))

glycopep_metadata$label <- paste(glycopep_metadata$condition, glycopep_metadata$replicate, sep = "_")

glycopep_sample_columns <- grep("raw.PEP.MS2Quantity", colnames(glycopep_filt_unique))
glycopep_filt_unique[glycopep_filt_unique == "NaN"] <- NA
glycopep_filt_unique[glycopep_sample_columns] <- sapply(glycopep_filt_unique[glycopep_sample_columns],as.numeric)

colnames(glycopep_filt_unique)[glycopep_sample_columns] <- glycopep_metadata$label
rownames(glycopep_filt_unique) = glycopep_filt_unique$name

se_glycopep <- DEP2::make_se(glycopep_filt_unique, columns = glycopep_sample_columns, expdesign = glycopep_metadata)

#df_se_glycopep <- get_df_wide(se_glycopep)

#write.xlsx(df_se_glycopep,"/data/p21_glycopep4008_raw_050726.xlsx")

#saveRDS(se_glycopep, "/data/glycobrain_docs/se_glycopep4008_220626.rds")

se_glycopep_filt <- DEP2::filter_se(se_glycopep, thr = 0)

#df_se_glycopep_filt <- get_df_wide(se_glycopep_filt)

#write.xlsx(df_se_glycopep_filt,"/data/p21_glycopep4008_filtered_050726.xlsx")

se_glycopep_vsn <- normalize_vsn(se_glycopep_filt)

#DEP2::plot_pca(se_glycopep_vsn,  indicate = c("condition", "timepoint"), n = 1000, label = TRUE)
#plot_normalization(se_glycopep_filt, se_glycopep_vsn)

### Harmonizer --------------------------------------------------------------

harmonizR(se_glycopep_vsn, 
          output_file = "/data/glycobrain_docs/glycopep_QScorefilt_vsn_Harmn_ouput1",
          verbosity = 4, algorithm = "ComBat", ComBat_mode = 3)

cured_tsv <- read_tsv("/data/glycobrain_docs/glycopep_QScorefilt_vsn_Harmn_ouput1.tsv")
colnames(cured_tsv)[1] <- "name"

df_se_glycopep_vsn <- get_df_wide(se_glycopep_vsn)
df_glycopep_meta_columns <- df_se_glycopep_vsn[c("name", 
                                                 "PG.ProteinAccessions", 
                                                 "PG.Genes", 
                                                 "PEP.StrippedSequence", 
                                                 "PEP.PeptidePosition" , 
                                                 "EG.ModifiedSequence",  
                                                 "ID",
                                                 "Glycan",
                                                 "epitope")]

df_post_harmonizer <- merge(df_glycopep_meta_columns,
                            cured_tsv,
                            by.x = "name",
                            by.y = "name",
                            all.x = FALSE,
                            all.y = TRUE)

df_post_harmonizer_sample_columns <- grep("_", colnames(df_post_harmonizer))

#experimental_design <- glycopep_metadata

se_glycopep_post_harmonizer_vsnfirst <- DEP2::make_se(df_post_harmonizer, 
                                                      columns = df_post_harmonizer_sample_columns, 
                                                      expdesign = glycopep_metadata, 
                                                      log2transform = FALSE)

## RF impute ---------------------------------------------------------------

set.seed(4)
se_glycopep_RF <- DEP2::impute(se_glycopep_post_harmonizer_vsnfirst, 
                               fun = "RF", ntree = 300, mtry = 12, verbose = TRUE, replace = FALSE, decreasing = TRUE)

saveRDS(se_glycopep_RF, 
        "/data/glycobrain_docs/se_glycopep4008_vsn_RF_230426.rds")