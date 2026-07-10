# ETD data ----------------------------------------------------------------

NoL_data_raw <- read.xlsx("/data/mBrain_GlycoSites_Summ_paper1_210526.xlsx",
                          sheet = "NoL")

TMT_data_raw <- read.xlsx("/data/mBrain_GlycoSites_Summ_paper1_210526-1.xlsx",
                          sheet = "TMT-canonical")



TMT_data <- TMT_data_raw %>%
  dplyr::select(-c(id, `N-`, Sequence, `C-`, Entry, Gene, UniProt.Name, local.position)) %>% #glimpse()
  dplyr::rename(Accession = Accessions,
                sugar = Glyco,
                AA = AminoAcid,
                position = Global.Position) %>%
  mutate(across(where(is.character), str_trim))

NoL_data <- NoL_data_raw %>%
  dplyr::select(-Description) %>%
  separate_longer_delim(c("Modifications"), delim = "];") %>%
  separate_wider_delim(c("Modifications"), delim = " [", names = c("sugar", "site")) %>%
  separate_longer_delim(c("site"), delim = "; ") %>%
  mutate(site = str_remove_all(site, "\\]")) %>%
  dplyr::filter(str_detect(site, "\\(")) %>%
  separate_wider_delim(c("site"), delim = "(", names = c("acceptor", "score")) %>% 
  mutate(score =  as.numeric(str_remove_all(score, "\\)"))) %>%
  dplyr::filter(score > 90) %>%
  mutate(across(where(is.character), str_trim)) %>% 
  dplyr::select(-score) %>% 
  distinct() %>%
  mutate(AA = str_remove_all(acceptor, "[0-9]"),
         position = as.numeric(str_remove_all(acceptor, "[A-Z]"))) %>%
  dplyr::select(-acceptor)

positional_data_comb <- rbind(NoL_data, TMT_data) %>%
  distinct()

positional_data_comb %>%
  dplyr::select(-sugar) %>% 
  distinct() %>% 
  nrow()