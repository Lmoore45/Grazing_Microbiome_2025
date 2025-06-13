##### Step 1: Load Libraries and Source VIP Function
library(ggplot2)
library(tidyverse)
library(compositions)  # For CLR transformation
library(mixOmics)
library(paletteer)
source("6.0_sPLS_VIP.R")  # Custom script to compute VIPs

##### Step 2: Load Data
metadata <- read_csv("1.1_Metadata.csv")
asv_16S <- read_csv("1.2_16S_ASV_Table.csv")

##### Step 3: Preprocess ASV Table
# Pivot ASV table to long format, then wide with samples as rows
asv_long <- asv_16S %>% dplyr::select(-Taxon) %>% pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Abundance")
asv_wide <- asv_long %>% pivot_wider(names_from = ASV, values_from = Abundance)

# Filter metadata to match ASV samples and vice versa
metadata_filtered <- metadata %>% dplyr::select(Sample, maom_c) %>% drop_na() %>% filter(Sample %in% asv_wide$Sample)
asv_wide_filtered <- asv_wide %>% filter(Sample %in% metadata_filtered$Sample)

##### Step 4: Normalize ASV Data with CLR
# CLR transformation requires no zeros, so replace with pseudocount
asv_numeric <- asv_wide_filtered %>% dplyr::select(-Sample)
asv_matrix <- as.matrix(asv_numeric)
asv_matrix[asv_matrix == 0] <- 0.0001
asv_clr <- clr(asv_matrix)

# Re-add sample column
asv_clr_df <- cbind(Sample = asv_wide_filtered$Sample, as.data.frame(asv_clr))

##### Step 5: Align Samples
# Ensure metadata and ASV data are aligned by sample
metadata_filtered <- metadata_filtered %>% arrange(match(Sample, asv_clr_df$Sample))
stopifnot(all(metadata_filtered$Sample == asv_clr_df$Sample))

##### Step 6: Prepare Response Variable
# Log-transform MAOC and check normality
response_variable <- log(metadata_filtered$maom_c)
hist(response_variable, main = "Histogram of Log-Transformed MAOM", xlab = "Log(MAOM)", breaks = 20)
print(shapiro.test(response_variable))

##### Step 7: Run sPLS Regression
# This section defines the key parameters for the sPLS model.
# 'ncomp' specifies the number of components used to summarize variation in the predictor matrix.
# 'keepX' determines how many ASVs are retained for each component, allowing sparsity in variable selection.
# These parameters should be refined based on cross-validation to balance model interpretability and predictive performance.
ASV_table <- as.matrix(asv_clr_df %>% dplyr::select(-Sample))
ncomp <- 3
keepX <- rep(25, ncomp)
spls_result <- spls(X = ASV_table, Y = response_variable, ncomp = ncomp, keepX = keepX, mode = "regression", scale = TRUE)
print(spls_result)

##### Step 8: Cross-Validate Components
# Use 10-fold cross-validation (repeated 10 times) to evaluate model performance across components.
# The perf function returns Mean Squared Error of Prediction (MSEP).
# The component with the lowest MSEP typically offers the best trade-off between predictive accuracy and model simplicity.
# These results help guide the selection of an optimal number of components (ncomp) for interpretation.
perf_result <- perf(spls_result, validation = "Mfold", folds = 10, progressBar = TRUE, nrepeat = 10)
plot(perf_result)

##### Step 9: Predicted vs Observed
# Extract and plot predicted values vs observed to assess model fit
predicted <- predict(spls_result, newdata = ASV_table)$predict[, 1, 1]
df_pred_obs <- data.frame(Observed = response_variable, Predicted = predicted)

# Plot predicted vs observed with linear fit
ggplot(df_pred_obs, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Observed (Component 1)", x = "Observed", y = "Predicted") +
  theme_minimal()

print(paste("R-squared:", summary(lm(Predicted ~ Observed, data = df_pred_obs))$r.squared))

##### Step 10: Extract and Plot VIP Scores
# Extract VIP scores for the first component and keep only ASVs with VIP > 1
vip_scores_full <- mixOmics::vip(spls_result)[, "comp1"]
vip_scores <- vip_scores_full[vip_scores_full > 1]
vip_sorted <- sort(vip_scores, decreasing = TRUE)
df_vip <- data.frame(ASV = names(vip_sorted), VIP = vip_sorted)

# Plot top 25 VIP ASVs
ggplot(df_vip[1:25,], aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 25 VIP Scores (Component 1)", x = "ASV", y = "VIP") +
  theme_minimal()

##### Step 11: Calculate R2 for High VIP ASVs
# Evaluate how strongly individual ASVs are correlated with the response variable
vip_high <- names(vip_scores[vip_scores > 1])
r2_results <- purrr::map_dfr(vip_high, ~{
  asv_abund <- ASV_table[, .x]
  r2 <- summary(lm(response_variable ~ asv_abund))$r.squared
  tibble(ASV = .x, R2 = r2)
})

##### Step 12: Integrate Taxonomy
# Split taxonomy strings and join with VIP and R2 data
taxonomy <- asv_16S %>% dplyr::select(1:2)
taxonomy_split <- taxonomy %>%
  separate(Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~str_remove(., ".*__"))) %>%
  replace_na(list(Domain="Unknown", Phylum="Unknown", Class="Unknown", Order="Unknown", Family="Unknown", Genus="Unknown", Species="Unknown"))

plot_data <- df_vip %>% dplyr::rename(ASV = ASV) %>% left_join(taxonomy_split, by = "ASV") %>% left_join(r2_results, by = "ASV")

##### Step 13: Plot VIP vs R2 with Taxonomy
# Use original curved arrow plot with phylum color palette for interpretability
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
# Normalize to relative abundance and extract top ASVs by VIP
asv_16S_relab <- asv_16S %>% dplyr::select(-Taxon) %>% mutate(across(-ASV, ~ . / sum(. , na.rm = TRUE)))
top_asvs <- df_vip %>% filter(ASV %in% vip_high) %>% pull(ASV)
filtered_asv <- asv_16S_relab %>% filter(ASV %in% top_asvs)

# Reshape and join metadata
transposed_asv <- filtered_asv %>%
  pivot_longer(cols = -ASV, names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = ASV, values_from = Count)

transposed_asv_maoc <- transposed_asv %>% left_join(metadata %>% dplyr::select(Sample, maom_c), by = "Sample")

# Calculate slope direction for each ASV
asv_relationships <- transposed_asv_maoc %>%
  pivot_longer(cols = -c(Sample, maom_c), names_to = "ASV", values_to = "RelAbundance") %>%
  group_by(ASV) %>%
  summarize(Coefficient = coef(lm(maom_c ~ RelAbundance))[2], .groups = "drop") %>%
  mutate(Slope = ifelse(Coefficient >= 0, "Positive", "Negative"))

plot_data <- plot_data %>% left_join(asv_relationships, by = "ASV")
transposed_asv <- transposed_asv %>% dplyr::select(-one_of(asv_relationships$ASV[asv_relationships$Slope == "Negative"]))

# Calculate sum and mean abundance of positively associated ASVs per sample
asv_sums <- transposed_asv %>% rowwise() %>% mutate(Sum_VIP = sum(c_across(-Sample), na.rm = TRUE)) %>% ungroup()
asv_means <- transposed_asv %>% rowwise() %>% mutate(Mean_VIP = mean(c_across(-Sample), na.rm = TRUE)) %>% ungroup()

# Export as needed
# write_csv(asv_sums, "16s_summed_relab_vip.csv")
# write_csv(asv_means, "16s_mean_relab_vip.csv")
