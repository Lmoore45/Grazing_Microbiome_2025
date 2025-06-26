##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R") # Custom script to compute VIPs

##### Step 2: Load Data
# Load metadata and ITS ASV table
metadata <- read_csv("1.1_Metadata.csv")
asv_ITS <- read_csv("1.3_ITS_ASV_Table.csv")

##### Step 3: Preprocess ASV Table
# Convert ITS ASV table to long and then wide format
asv_long <- asv_ITS %>% dplyr::select(-Taxon) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata and ASV table for matching samples
metadata_filtered <- metadata %>% dplyr::select(Sample, maom_c) %>% drop_na() %>% filter(Sample %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# Prepare matrix, replace zeros, apply CLR transformation
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
asv_matrix[asv_matrix == 0] <- 0.0001
asv_clr <- clr(asv_matrix)

# Combine with sample names
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure ASV and metadata rows match in order
metadata_filtered <- metadata_filtered %>% arrange(match(Sample, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable
# Log-transform MAOM carbon, check for normality
response_variable <- log(metadata_filtered$maom_c)
hist(response_variable, main = "Histogram of Log-Transformed MAOM", xlab = "Log(MAOM)", breaks = 20)
print(shapiro.test(response_variable))

##### Step 7: Run sPLS Regression
# Define model structure and sparsity
# Adjust ncomp and keepX based on cross-validation results
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
ncomp <- 3
keepX <- rep(25, ncomp)
spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold repeated CV to determine optimal number of components
perf_result <- perf(spls_result, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)
plot(perf_result)

##### Step 9: Predicted vs Observed
# Evaluate model performance visually and quantitatively
predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

##### Step 10: Extract and Plot VIP Scores
# Filter for ASVs with VIP > 1 for interpretation and downstream analysis
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

ggplot(df_vip[1:22,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 25 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

##### Step 11: Calculate R2 for High VIP ASVs
# Quantify strength of association between individual ASVs and response
vip_high <- names(vip_scores)
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Clean and merge ITS taxonomy table with VIP and R2 results
taxonomy <- asv_ITS %>%  dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  mutate(across(everything(), ~replace_na(., "Unknown")))

plot_data <- df_vip %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Generate main interpretation plot with curved arrows
maoc_plot <- ggplot() +
  geom_segment(aes(x = 0, xend = max(plot_data$VIP, na.rm = TRUE) * 1.1, y = 0, yend = 0), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max(plot_data$R2, na.rm = TRUE) * 1.1), arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8, color = "black") +
  geom_curve(data = plot_data, aes(x = VIP, y = 0, xend = 0, yend = R2, color = Phylum), curvature = 0.4, alpha = 1, linewidth = 0.8) +
  geom_point(data = plot_data, aes(x = VIP, y = 0, fill = Phylum), shape = 21, size = 4, color = "black") +
  geom_point(data = plot_data, aes(x = 0, y = R2, fill = Phylum), shape = 21, size = 4, color = "black") +
  scale_fill_paletteer_d("ltc::crbhits", name = "Phylum") +
  scale_color_paletteer_d("ltc::crbhits", guide = "none") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  labs(title = "VIP Scores vs R² Values by Phylum", x = "VIP Score", y = "Correlation to Mineral-Associated Organic Carbon (R2)") +
  scale_x_continuous(limits = c(0, max(plot_data$VIP, na.rm = TRUE) * 1.1)) +
  scale_y_continuous(limits = c(0, max(plot_data$R2, na.rm = TRUE) * 1.1))

print(maoc_plot)

##### Step 14: Save Top VIP Relative Abundance (Sum and Mean)
# Normalize ITS table to relative abundance
asv_ITS_relab <- asv_ITS %>% mutate(across(-ASV, ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip$ASV
filtered_asv <- asv_ITS_relab %>% filter(ASV %in% top_asvs)

# Reshape and align with metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count) %>%
  dplyr::slice(-89)  # Drop problematic sample if necessary

transposed_asv_maoc <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, maom_c), by = "Sample")

# Identify slope direction per ASV
asv_relationships <- transposed_asv_maoc %>%
  pivot_longer(cols = -c(Sample, maom_c), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(maom_c ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope_Direction = ifelse(Coefficient >= 0, "Positive", "Negative"))

plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope_Direction == "Negative"]))

# Calculate summed and mean abundances per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Optional: write to CSV
# write_csv(asv_sums, "ITs_summed_relab_vip.csv")
# write_csv(asv_means, "ITs_mean_relab_vip.csv")
