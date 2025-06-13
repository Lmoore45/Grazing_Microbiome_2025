# DISCLOSURE: Due to rancher privacy agreements, precise geographic coordinates used in the manuscript analyses are not included in this repository. The version of the code shared here uses county-level centroids as placeholders. 
# As a result, the results and figure from this analysis differ slightly from those in the published manuscript. However, all primary findings and interpreations remain consistent.

##### Step 1: Load Required Libraries
library(tidyverse)
library(vegan)
library(ggplot2)
library(car)
library(VennDiagram)

##### Step 2: Load Metadata
metadata <- read_csv("1.1_Metadata.csv")

##### Step 3: Define Predictor Groups
# Microbiome predictors
microbiome_predictors <- metadata %>%
  dplyr::select(Sample, Richness_16s, Evenness_16s, Shannons_16s, Biomass, FB_Ratio_Normal,
                Richness_ITS, Evenness_ITS, Shannons_ITS, vip_sum_16S_maoc_relab, vip_sum_ITS_maoc_relab) %>%
  mutate(across(-Sample, as.numeric))

# Physicochemical predictors
physicochemical_predictors <- metadata %>%
  dplyr::select(Sample, pom_c, total_n, maom_n, pom_n, pH, percent_clay, percent_sand, percent_silt, c_n_ratio)

# Management predictors
management_predictors <- metadata %>% dplyr::select(Sample, Adaptive_Score)

# Climate predictors
climate_predictors <- metadata %>%
  dplyr::select(Sample, mean_annual_temp, mean_annual_precip, PET_mm_day, radiation, Latitude, Longitude)

##### Step 4: Remove NA Values and Align Samples
# Drop NA values
microbiome_predictors <- na.omit(microbiome_predictors)
physicochemical_predictors <- na.omit(physicochemical_predictors)
management_predictors <- na.omit(management_predictors)
climate_predictors <- na.omit(climate_predictors)

# Align by common samples
common_samples <- Reduce(intersect, list(
  microbiome_predictors$Sample,
  physicochemical_predictors$Sample,
  management_predictors$Sample,
  climate_predictors$Sample
))

microbiome_predictors <- microbiome_predictors %>% filter(Sample %in% common_samples)
physicochemical_predictors <- physicochemical_predictors %>% filter(Sample %in% common_samples)
management_predictors <- management_predictors %>% filter(Sample %in% common_samples)
climate_predictors <- climate_predictors %>% filter(Sample %in% common_samples)

##### Step 5: Set Row Names and Remove Sample Columns
set_rownames <- function(df) {
  df <- as.data.frame(df)               # Convert tibble to data.frame
  rownames(df) <- df$Sample             # Assign row names
  df <- df[, !names(df) %in% "Sample"]  # Drop Sample column
  return(df)
}

microbiome_predictors <- set_rownames(microbiome_predictors)
physicochemical_predictors <- set_rownames(physicochemical_predictors)
management_predictors <- set_rownames(management_predictors)
climate_predictors <- set_rownames(climate_predictors)

##### Step 6: Scale Predictor Variables
scale_df <- function(df) {
  scaled <- as.data.frame(scale(df))
  rownames(scaled) <- rownames(df)
  return(scaled)
}

microbiome_scaled <- scale_df(microbiome_predictors)
physicochemical_scaled <- scale_df(physicochemical_predictors)
management_scaled <- scale_df(management_predictors)
climate_scaled <- scale_df(climate_predictors)

##### Step 7: Prepare Response Variable
response_variable_df <- metadata %>%
  dplyr::select(Sample, maom_c) %>%
  filter(Sample %in% rownames(microbiome_scaled), maom_c > 0) %>%  # Only keep positive values
  mutate(maom_c = log(maom_c))

response_variable <- response_variable_df$maom_c
names(response_variable) <- response_variable_df$Sample

hist(response_variable, main = "Histogram of Log-Transformed MAOM", xlab = "Log(MAOM)", breaks = 20)
print(shapiro.test(response_variable))

##### Step 8: Check Pairwise Correlations 
# This step identifies highly collinear variables (|r| > 0.9) within each predictor category.
# Collinearity can distort regression estimates and inflate the apparent influence of certain predictors.
check_correlation <- function(df, label) {
  cor_matrix <- cor(df, use = "pairwise.complete.obs")
  high_corr <- which(abs(cor_matrix) > 0.9, arr.ind = TRUE)
  high_corr <- high_corr[high_corr[,1] != high_corr[,2], ]
  cat("\nHighly correlated pairs in", label, ":\n")
  print(high_corr)
}

check_correlation(microbiome_scaled, "Microbiome")
check_correlation(physicochemical_scaled, "Physicochemical")
check_correlation(management_scaled, "Management")
check_correlation(climate_scaled, "Climate")

##### Step 9: Remove Collinear Predictors Based on Prior Inspection
# Based on the correlation matrices (Step 8), we manually remove variables that are highly collinear (|r| > 0.9).
# Removing these variables helps prevent redundancy and instability in downstream regression models.
# This step improves model interpretability and avoids inflated variance in coefficient estimates.
# You may need to revisit this step after calculating VIF scores in Step 10 to further refine the predictor set.
microbiome_clean <- microbiome_scaled %>% select(-Richness_16s)
physicochemical_clean <- physicochemical_scaled %>% select(-pom_n, -percent_silt, -maom_n, -total_n)
climate_clean <- climate_scaled %>% select(-Latitude, -radiation, -PET_mm_day)

##### Step 10: Calculate Variance Inflation Factor (VIF)
# VIF quantifies how much the variance of a regression coefficient is inflated due to multicollinearity.
# As a rule of thumb, variables with VIF > 5 may introduce instability in the model due to collinearity with other predictors.
# If any variables exceed this threshold, return to Step 9 and remove those variables from the predictor set.
# It's important to re-check VIF after any removal to ensure remaining variables are stable.

check_vif <- function(df, y, label) {
  model <- lm(y ~ ., data = df)
  vif_scores <- car::vif(model)
  cat("\n", label, "VIF:\n")
  print(vif_scores)
}

check_vif(microbiome_clean, response_variable, "Microbiome")
check_vif(physicochemical_clean, response_variable, "Physicochemical")
check_vif(climate_clean, response_variable, "Climate")

##### Step 11: Confirm Sample Alignment
common_samples <- Reduce(intersect, list(
  rownames(microbiome_clean),
  rownames(physicochemical_clean),
  rownames(climate_clean),
  rownames(management_scaled),
  names(response_variable)
))

##### Step 12: Variance Partitioning
varpart_results <- varpart(
  response_variable,
  microbiome_clean,
  physicochemical_clean,
  climate_clean,
  management_scaled
)

print(varpart_results)
plot(varpart_results, main = "Variance Partitioning of Predictors")

##### Step 13: Extract and Summarize Variance Fractions
# Results from 'final_df' were used to generate the upset plot included in the manuscript.
# This table combines both individual and combined variance contributions from each predictor category.

fract <- varpart_results$part$fract
indfract <- varpart_results$part$indfract

individual_names <- c("Microbiome (X1 | X2+X3+X4)",
                      "Physicochemical (X2 | X1+X3+X4)",
                      "Climate (X3 | X1+X2+X4)",
                      "Management (X4 | X1+X2+X3)")

indfract_df <- data.frame(
  Fraction = rownames(indfract)[1:4],
  Predictor_Category = individual_names,
  Df = indfract$Df[1:4],
  R.square = indfract$Adj.R.square[1:4],
  Adj.R.square = NA,
  Testable = indfract$Testable[1:4]
)

main_names <- c("Microbiome", "Physicochemical", "Climate", "Management",
                "Microbiome + Physicochemical", "Microbiome + Climate",
                "Microbiome + Management", "Physicochemical + Climate",
                "Physicochemical + Management", "Climate + Management",
                "Microbiome + Physicochemical + Climate",
                "Microbiome + Physicochemical + Management",
                "Microbiome + Climate + Management",
                "Physicochemical + Climate + Management", "All Combined")

fract_df <- data.frame(
  Fraction = rownames(fract),
  Predictor_Category = main_names,
  Df = fract$Df,
  R.square = fract$R.square,
  Adj.R.square = fract$Adj.R.square,
  Testable = fract$Testable
)

final_df <- rbind(fract_df, indfract_df)
print(final_df)

##### Step 14: Redundancy Analysis (RDA) and Permutation Test
# Test whether adding Microbiome predictors significantly improves model over Physicochemical-only
rda_phys <- rda(response_variable ~ ., data = physicochemical_clean)
rda_phys_micro <- rda(response_variable ~ ., data = cbind(physicochemical_clean, microbiome_clean))

anova_result <- anova(rda_phys, rda_phys_micro, permutations = 9999)
print(anova_result)

cat("\nRDA Adjusted R² Values:\n")
cat("Physicochemical-only: R² =", round(RsquareAdj(rda_phys)$r.squared, 3),
    ", Adj-R² =", round(RsquareAdj(rda_phys)$adj.r.squared, 3), "\n")
cat("Phys + Microbiome:    R² =", round(RsquareAdj(rda_phys_micro)$r.squared, 3),
    ", Adj-R² =", round(RsquareAdj(rda_phys_micro)$adj.r.squared, 3), "\n")

