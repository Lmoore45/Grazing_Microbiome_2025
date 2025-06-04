# PCA Analysis for Soil Physicochemical and Microbial Data
# Author: [Laura Moore]
# Date: [2025-02-10]
# Description: This script performs PCA using FactoMineR and visualizes results with factoextra and ggplot2.



# Load necessary packages
library(FactoMineR)  # For performing PCA
library(factoextra)  # For visualizing PCA results
library(ggplot2)     # For general plotting
library(dplyr)       # For data manipulation
library(tidyverse)   # For general data wrangling
library(ggpubr)      # For arranging plots
library(grid)        # For graphical functions


# Load metadata
data <- read.csv("metadata.csv")

# Filter out the outlier and select relevant variables
metadata <- data %>%
  filter(Sample != "WYCO_0923_05_PR_D3_B_C03_XT") %>%
  select(Sample, pH, total_c, maom_c, pom_c, total_n, maom_n, pom_n,
         percent_clay, percent_sand, percent_silt,
         mean_annual_temp, mean_annual_precip, PET, radiation, Location)

# Remove missing values
metadata[metadata == "#N/A"] <- NA
metadata <- na.omit(metadata)

# Log-transform non-normally distributed variables
metadata <- metadata %>%
  mutate(
    total_c = log10(total_c + 1),
    maom_c = log10(maom_c + 1),
    pom_c = log10(pom_c + 1),
    maom_n = log10(maom_n + 1),
    pom_n = log10(pom_n + 1)
  )

# Perform PCA, ensuring only numeric variables are used
pca_results <- PCA(metadata %>% select(-Sample, -Location), scale.unit = TRUE, graph = FALSE)

# Scree plot to visualize variance explained
fviz_eig(pca_results)

# Extract PCA results for variables
var_results <- get_pca_var(pca_results)
print(var_results)  # Check variable loadings

# Ensure 'Location' matches the filtered metadata
metadata$Location <- as.factor(metadata$Location)

# Biplot colored by location (Site) for PCA visualization
plot2 <- fviz_pca_biplot(pca_results,
                         col.ind = metadata$Location,   # Color points by location (filtered)
                         palette = "Set4",             # Custom color palette
                         col.var = "#4D4D4D",          # Color for variable vectors
                         repel = TRUE,                 # Prevent overlapping labels
                         label = "var",                # Label variables only
                         arrowsize = 0.5,              # Adjust arrow size
                         pointsize = 2,                # Adjust point size
                         pointshape = 19,              # Maintain consistent shape
                         labelsize = 5) +              # Adjust label size for variables
  labs(color = "Location")  # Add title to legend

print(plot2)

# Extract individual PCA scores (sample points)
ind_scores <- as.data.frame(pca_results$ind$coord)  # Individual scores

# Ensure `ind_scores` matches `metadata` after filtering
ind_scores$Location <- metadata$Location  # Use filtered metadata for Location

# Extract `Practice` from `data` ensuring it matches `metadata`
ind_scores$Practice <- data$Practice[match(metadata$Sample, data$Sample)] 

# Extract PCA variable loadings (arrows)
var_scores <- as.data.frame(pca_results$var$coord)  # Variable loadings
var_scores$Variable <- rownames(var_scores)        # Assign variable names

# Rescale variable scores to match sample score space
scale_factor <- max(abs(ind_scores[, 1:2])) / max(abs(var_scores[, 1:2]))
var_scores[, 1:2] <- var_scores[, 1:2] * scale_factor  

# Apply reduction factor to control arrow length in the plot
reduction_factor <- 0.8  
var_scores[, 1:2] <- var_scores[, 1:2] * reduction_factor

# PCA plot using ggplot, coloring samples by Location and shaping by Practice
pca_ggplot <- ggplot() +
  # Sample points (colored by location, shaped by practice)
  geom_point(data = ind_scores,  
             aes(x = Dim.1, y = Dim.2, color = Location, shape = Practice),
             size = 4, alpha = 0.8) +
  # Arrows for PCA variables
  geom_segment(data = var_scores,
               aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  # Variable labels
  geom_text(data = var_scores, aes(x = Dim.1, y = Dim.2, label = Variable),
            color = "black", hjust = 1.2, size = 4) +
  # Dotted lines for PC1 and PC2 axes
  geom_hline(yintercept = 0, linetype = "dotted", size = 1.2, color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", size = 1.2, color = "black") +
  # Manually define colors for each site
  scale_color_manual(name = "Location (Site)", 
                     values = c("Site 13" = "#be7245",  # Site A2
                                "Site 15" = "#46211c",  # Site A3
                                "Site 16" = "#d2ad7c",  # Site A1
                                "Site 19" = "#D8DAEB",  # Site C1
                                "Site 5"  = "#998ec3",  # Site C2
                                "Site 6"  = "#542788"   # Site C3
                     )) +
  # Axis labels and legends
  labs(x = "PC1", y = "PC2", shape = "Practice") +
  # Theme adjustments
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, size = 1.2),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black")
  )

print(pca_ggplot)


