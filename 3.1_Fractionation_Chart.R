library(dplyr)
library(ggplot2)
library(tidyr)
library(ggthemes)
library(ggpubr)

# Load your metadata
metadata <- read.csv("1.1_Metadata.csv")

### ==== BAR CHART SECTION ====

# Calculate total SOC
metadata <- metadata %>%
  select(Sample, total_c, maom_c, pom_c, Location, depth, Practice) %>% 
  drop_na()

# Prepare long format for bar chart
summary_df <- metadata %>%
  select(Location, pom_c, maom_c, total_c) %>%
  drop_na() %>%
  pivot_longer(cols = c(pom_c, maom_c, total_c),
               names_to = "Fraction",
               values_to = "Value") %>%
  mutate(
    Fraction = recode(Fraction,
                      "pom_c" = "POC",
                      "maom_c" = "MAOC",
                      "total_c" = "SOC"),
    Location = factor(Location,
                      levels = c("Site 19", "Site 5", "Site 6", "Site 16", "Site 13", "Site 15")),
    Fraction = factor(Fraction, levels = c("SOC", "MAOC", "POC"))
  )

# Summary stats for bar chart
bar_summary <- summary_df %>%
  group_by(Location, Fraction) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    se = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Bar chart
bar_plot <- ggplot(bar_summary, aes(x = Fraction, y = mean_value, fill = Fraction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se),
                width = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(~ Location, nrow = 1) +
  labs(x = "Carbon Fraction", y = "Carbon Concentration (mg C/g soil)",
       title = "Distribution of SOC Fractions Within Each Site") +
  scale_fill_manual(values = c("POC" = "gray60", "MAOC" = "#A3A990", "SOC" = "#BB946C")) +
  theme_tufte(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(color = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

print(bar_plot)

### ==== BOX PLOT SECTION (Supplemental) ====

# Subset and pivot longer for box plot
data_long <- metadata %>%
  select(depth, Practice, pom_c, maom_c, total_c) %>%
  pivot_longer(cols = c(pom_c, maom_c, total_c),
               names_to = "Fraction", values_to = "Value") %>%
  mutate(
    Fraction = recode(Fraction,
                      "pom_c" = "POC",
                      "maom_c" = "MAOC",
                      "total_c" = "SOC"),
    Fraction = factor(Fraction, levels = c("SOC", "MAOC", "POC")),
    Practice = factor(Practice, levels = c("Prescriptive", "Adaptive"))
  )

# Calculate p-values (but don't plot them)
p_values <- data_long %>%
  group_by(depth, Fraction) %>%
  summarise(p_value = t.test(Value ~ Practice)$p.value,
            .groups = "drop")

# Reorder the fraction levels
data_long$Fraction <- factor(data_long$Fraction, levels = c("MAOC", "POC", "SOC"))

# Updated box plot (no p-values shown)
box_plot <- ggplot(data_long, aes(x = Practice, y = Value, fill = Fraction)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 0.6) +
  # Optional jitter line is left commented out
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  #             size = 1.5, alpha = 0.6, color = "black") +
  facet_grid(Fraction ~ depth, scales = "free_y", switch = "y") +
  labs(x = "Practice", y = "SOC Concentration (mg C/g soil)", fill = "Carbon Fraction") +
  theme_tufte(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(color = "black", size = 1),
    strip.placement = "outside"
  ) +
  scale_fill_manual(values = c("POC" = "gray60", "MAOC" = "#A3B18A", "SOC" = "#B35806"))

# Print the plot
print(box_plot)


