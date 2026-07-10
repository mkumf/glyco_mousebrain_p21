# CV curves ---------------------------------------------------------------


se_glycopep_RF <- readRDS("/data/glycobrain_docs/se_glycopep4008_vsn_RF_230426.rds")

df_assay_se_glycopep <- assay(se_glycopep_RF)

data_long <- data.frame(df_assay_se_glycopep) %>%
  tibble::rownames_to_column(var = "Gene") %>%
  #  pivot_longer( everything(),names_to = "region", values_to = "abundance") %>%
  pivot_longer(cols = -Gene,
               names_to = c("Region", "sample"),
               names_sep = "_",
               values_to = "intensity")

#data_long <- data_long %>% mutate(Region = str_replace_all(Region, region_replacement_pairs))

data_long$intensity <- 2^as.numeric(data_long$intensity)

# Calculate CV for each feature within each condition
cv_data <- data_long %>%
  group_by(Region, Gene) %>%
  summarise(
    mean_intensity = mean(as.numeric(intensity), na.rm = TRUE),
    sd_intensity = sd(as.numeric(intensity),  na.rm = TRUE),
    cv = sd_intensity / mean_intensity
  ) %>%
  ungroup() %>% na.omit()

# Function to calculate the exact percentage of data below each unique CV value
calculate_exact_percentages <- function(cv_data) {
  unique_cvs <- unique(cv_data$cv)
  percentages <- sapply(unique_cvs, function(threshold) {
    mean(cv_data$cv <= threshold, na.rm = TRUE) * 100
  })
  data.frame(cv_threshold = unique_cvs, percentage = percentages)
}

# Calculate exact percentages for each condition
percentage_data <- cv_data %>%
  group_by(Region) %>%
  reframe(calculate_exact_percentages(cur_data())) %>%
  unnest(cols = c(cv_threshold, percentage))

threshold_data <- percentage_data[percentage_data$percentage <= 80, ]
threshold_data_min <- percentage_data[percentage_data$percentage >= 80, ]
# Find the curve with the highest X value at the threshold
last_curve_to_reach_80 <- threshold_data[threshold_data$cv_threshold == max(threshold_data$cv_threshold), ]
first_curve_to_reach_80 <- threshold_data_min[threshold_data_min$cv_threshold == min(threshold_data_min$cv_threshold), ]


percentage_data$Region <- factor(percentage_data$Region, levels = c("Olfactory",
                                                                    "Cortex",
                                                                    "Hippocampus",
                                                                    "Thalamus",
                                                                    "Hypothalamus",
                                                                    "Midbrain",
                                                                    "Pons",
                                                                    "Medulla",
                                                                    "Cerebellum"))
color_regions <- c("Cortex"="#E31A1C",
                   "Olfactory"="#1F78B4",
                   "Cerebellum"="#B2DF8A",
                   "Thalamus"="#1B9E77",
                   "Hippocampus"="#666666",
                   "Midbrain"="#A6CEE3",
                   "Medulla"="#FF61C3",
                   "Hypothalamus"="#FF7F00",
                   "Pons"="#6A3D9A")

# Create the plot with continuous values for percentage
ggplot(percentage_data, aes(x = cv_threshold, y = percentage, color = Region)) +
  geom_hline(linewidth = .25, yintercept = 80, linetype = "dashed") +
  # geom_vline(linewidth = .25, xintercept = first_curve_to_reach_80$cv_threshold, linetype = "dashed") +
  # geom_vline(linewidth = .25, xintercept = last_curve_to_reach_80$cv_threshold, linetype = "dashed") +
  geom_line(linewidth = 0.5) +
  geom_point(data = last_curve_to_reach_80, aes(x = cv_threshold, y = percentage, fill = Region), shape = 21, color = "black") +
  geom_label(data = last_curve_to_reach_80, aes(label = paste0(format(round(cv_threshold, digits = 3), nsmall = 3))), color = "black", nudge_y = -5, nudge_x = +0.05, alpha = .75) +
  geom_point(data = first_curve_to_reach_80, aes(x = cv_threshold, y = percentage, fill = Region), shape = 21, color = "black") +
  geom_label(data = first_curve_to_reach_80, aes(label = paste0(format(round(cv_threshold, digits = 3), nsmall = 3))), color = "black", nudge_y = +5, nudge_x = -0.05,  alpha = .75) +
  labs(title = "Glycopeptide cumulative coefficient of variation (CV) distribution",
       x = "CV",
       y = " Cumulative CV, %",
       color = "Region") +
  scale_x_continuous(limits = c(0, 1.2), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100), labels = scales::label_number(suffix = "%")) +
  scale_color_manual(values = color_regions) +
  scale_fill_manual(values = color_regions, guide = "none") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5)))