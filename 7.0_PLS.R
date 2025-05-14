# DISCLOSURE: Due to rancher privacy agreements, precise geographic coordinates used in the manuscript analyses are not included in this repository. The version of the code shared here uses county-level centroids as placeholders. 
# As a result, the results and figure from this analysis differ slightly from those in the published manuscript. However, all primary findings and interpreations remain consistent.


# Load necessary libraries
library(tidyverse)           
library(mixOmics)          
library(pls)                  
library(corrplot)            
library(ggplot2)              
library(ggthemes)             

# Source custom VIP scoring script
source("/Users/lmoore/Desktop/Documents/CSU/PhD/Dissertation/Grazing Microbiome/R_Code/Manuscript_Code/VIP.R")


##### Step 1: Load Metadata
# Read in the metadata file containing sample-level environmental, microbial, and geographic data.
metadata <- read_csv("/Users/lmoore/Desktop/Documents/CSU/PhD/Dissertation/Grazing Microbiome/R_Code/Manuscript_Code/metadata.csv")

##### Step 2: Build Predictor Matrix
# Subset and prepare predictor variables for the PLS model.
# Add total_n and maom_n to observe collinearity in the correlation matrix.
combined_predictors <- metadata %>% 
  dplyr::select(Sample, Richness_16s, Evenness_16s, Shannons_16s, 
                Evenness_ITS, Shannons_ITS, Richness_ITS, Biomass, 
                FB_Ratio_Normal, vip_sum_16S_maoc_relab, vip_sum_ITS_maoc_relab, 
                pom_c, pom_n, c_n_ratio, pH, percent_clay, percent_sand, percent_silt, 
                Adaptive_Score, mean_annual_temp, mean_annual_precip, 
                PET_mm_day, radiation, Latitude, Longitude) %>% 
  drop_na() %>% 
  mutate(across(-Sample, as.numeric))

##### Step 3: Screen for Collinearity
# Generate a Spearman correlation matrix to identify highly collinear predictors (r > 0.9).
# Remove redundant variables prior to running the PLS model.
cor_matrix <- cor(combined_predictors %>% dplyr::select(-Sample), 
                  use = "pairwise.complete.obs", method = "spearman")


##### Step 4: Align Datasets
# Ensure predictors and metadata align on shared samples by intersecting sample IDs.
common_samples <- intersect(metadata$Sample, combined_predictors$Sample)
metadata_aligned <- metadata %>% filter(Sample %in% common_samples)
combined_predictors <- combined_predictors %>% 
  filter(Sample %in% common_samples) %>% 
  dplyr::select(-Sample) %>% 
  as.data.frame()

##### Step 5: Define and Normalize Response Variable
# Select the response variable (e.g., MAOC) and apply log transformation to normalize its distribution.
# Confirm normality using a Shapiro-Wilk test.
response_variable <- metadata_aligned$maom_c
response_variable <- log(response_variable)

hist(response_variable, main = "Histogram of Scaled Response", 
     xlab = "Log-transformed MAOM", breaks = 20)
shapiro_test <- shapiro.test(response_variable)
print(shapiro_test)

##### Step 7: Run PLS Regression
# Begin by testing 20 components to evaluate how well R2 and RMSEP explain variation in the response variable
# The optimal number of components is selected based on where R2 platues, and RMSEP stops decreasing, minimizing overfitting
# Once identified, update 'opimal_ncomp' and refit the final model
# Extract VIP scores from teh final model to evaluate variable importance.

optimal_ncomp <- 10
pls_model <- plsr(response_variable ~ ., 
                  data = combined_predictors, 
                  validation = "CV", 
                  scale = TRUE, 
                  method = "oscorespls", 
                  ncomp = optimal_ncomp)

summary(pls_model)
R2(pls_model, estimate = "CV")
RMSEP(pls_model)

vip_scores <- VIP(pls_model)

##### Step 8: Extract VIP Scores
# Create a table of variable importance scores (VIP) from the final PLS model.
# VIP scores summarize how strongly each predictor contributes to explaining the response.
plot_data <- data.frame(
  Variable = colnames(vip_scores),
  VIP_Score = vip_scores[optimal_ncomp, ]
)

##### Step 9: Compute Variable Significance
# For each predictor, run a univariate linear regression against the response variable.
# Extract R² and p-values, and apply FDR correction for multiple comparisons.
predictors <- colnames(combined_predictors)
stats_df <- data.frame(Variable = character(), p_value = numeric(), R2 = numeric())

for (var in predictors) {
  lm_model <- lm(response_variable ~ combined_predictors[[var]])
  r2_val <- summary(lm_model)$r.squared
  p_val <- summary(lm_model)$coefficients[2, 4]
  stats_df <- rbind(stats_df, data.frame(Variable = var, p_value = p_val, R2 = r2_val))
}

stats_df$p_value <- p.adjust(stats_df$p_value, method = "fdr")
plot_data <- left_join(plot_data, stats_df, by = "Variable")

##### Step 10: Annotate Significance for Plotting
# Categorize variables by significance thresholds (FDR-adjusted p-values) for visualization.
plot_data <- plot_data %>%
  mutate(
    p_value_category = case_when(
      p_value <= 0.001 ~ "p ≤ 0.001",
      p_value <= 0.01 ~ "p ≤ 0.01",
      p_value <= 0.05 ~ "p ≤ 0.05",
      TRUE ~ "p > 0.05"
    ),
    Color = case_when(
      p_value <= 0.001 ~ "#CBC106",
      p_value <= 0.01 ~ "#389CA7",
      p_value <= 0.05 ~ "#FFA500",
      TRUE ~ "gray10"
    )
  )

##### Step 11: Plot VIP Scores
# Generate a dot plot of VIP scores colored by p-value and sized by R².
plot <- ggplot(plot_data, aes(x = VIP_Score, y = reorder(Variable, VIP_Score), size = R2)) +
  geom_point(aes(color = p_value_category, fill = p_value_category),
             shape = 21, stroke = 0.8, color = "black") +
  geom_vline(xintercept = 1.0, linetype = "dotted", color = "black") +
  scale_color_manual(name = "P-value",
                     values = c("p ≤ 0.001" = "#CBC106", "p ≤ 0.01" = "#389CA7", 
                                "p ≤ 0.05" = "#FFA500", "p > 0.05" = "gray10")) +
  scale_fill_manual(name = "P-value",
                    values = c("p ≤ 0.001" = "#CBC106", "p ≤ 0.01" = "#389CA7", 
                               "p ≤ 0.05" = "#FFA500", "p > 0.05" = "gray10")) +
  scale_size_continuous(name = "R²", range = c(1, 7)) +
  labs(x = "VIP Scores", y = NULL, title = "PLS VIP Scores") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80")
  )

print(plot)
