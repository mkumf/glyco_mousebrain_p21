# plot normalization ------------------------------------------------------
#se_glycopep <- readRDS("/data/glycobrain_docs/se_glycopep4008_220626.rds")
se_glycopep_RF <- readRDS("/data/glycobrain_docs/se_glycopep4008_vsn_RF_230426.rds")

region_colors <- c("Cortex"="#E31A1C",
                   "Olfactory"="#1F78B4",
                   "Cerebellum"="#B2DF8A",
                   "Thalamus"="#1B9E77",
                   "Hippocampus"="#666666",
                   "Midbrain"="#A6CEE3",
                   "Medulla"="#FF61C3",
                   "Hypothalamus"="#FF7F00",
                   "Pons"="#6A3D9A")

df <- assay(se_glycopep) %>%
  data.frame() %>%
  tidyr::gather(ID, val) %>%
  left_join(., data.frame(colData(se_glycopep)), by = "ID")

condition_order <- c("Olfactory", "Cortex", "Hippocampus", "Thalamus", "Hypothalamus", "Midbrain", "Pons", "Medulla", "Cerebellum")

df$condition <- factor(df$condition, levels = condition_order)
df <- df %>%
  arrange(condition, ID) %>%
  mutate(ID = factor(ID, levels = unique(ID)))


p_raw <- ggplot(df, aes(x = ID, y = val, fill = condition)) +
  geom_boxplot(notch = TRUE, na.rm = TRUE, outlier.size = 1, outlier.alpha = 0.5, outlier.shape = 16) +
  scale_y_continuous(limits = c(0, 30)) +
  scale_fill_manual(values = region_colors,  name = "Region") +
  # coord_flip() +
  labs(title = "Glycopeptide abundance, raw, filtered", x = "Sample", y = expression(log[2]~"abundance")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#p_processed
(p_raw + p_processed) + plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a') &           
  theme(plot.tag = element_text(face = "bold", size = 14)) 
