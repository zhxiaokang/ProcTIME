
rm(list=ls())

library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(ggpubr)

# Load color palette
source("color_palette.R")

# Create output directory
out_dir <- "../output/stromal_validation"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Load data
message("Loading data...")
load("../data/clustering_tcga_removing_badguys.RData") # Contains df.clinical.bcr.cluster
load("../data/prad_data_tcga.RData") # Contains prad.gene
load("../data/ESTIMATE_TCGA.RData") # Contains res.deconv

# Align samples
common_samples <- intersect(rownames(df.clinical.bcr.cluster), colnames(prad.gene))
clinical_data <- df.clinical.bcr.cluster[common_samples, ]
expression_data <- prad.gene[, common_samples]

# Log transform if necessary (same check as before)
if (max(expression_data, na.rm = TRUE) > 100) {
  message("Applying log2(x+1) transformation...")
  expression_data <- log2(expression_data + 1)
}

# 1. Visualize CAF-specific markers
# List of common CAF/Stromal markers
caf_markers <- c("ACTA2", "FAP", "PDGFRB", "POSTN", "COL1A1", "FN1", "VIM", "DCN", "LUM")

# Check which markers are in the dataset
valid_markers <- intersect(caf_markers, rownames(expression_data))
message(paste("Found", length(valid_markers), "CAF markers out of", length(caf_markers)))

# Prepare data for plotting
marker_expr <- t(expression_data[valid_markers, ])
plot_data <- as.data.frame(marker_expr)
plot_data$sample <- rownames(plot_data)
plot_data$cluster <- clinical_data[plot_data$sample, "cluster"]

# Convert to long format for faceting
plot_data_long <- plot_data %>%
  pivot_longer(cols = all_of(valid_markers), names_to = "Gene", values_to = "Expression")

# Set cluster levels
plot_data_long$cluster <- factor(plot_data_long$cluster, levels = c("TCE", "EPCE", "TASCE"))

# Plot Boxplots for Markers
p_markers <- ggplot(plot_data_long, aes(x = cluster, y = Expression, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = color.three.clusters) +
  labs(title = "Expression of CAF-specific Markers", y = "Log2 Expression", x = "Subtype") +
  stat_compare_means(comparisons = list(c("TCE", "TASCE"), c("EPCE", "TASCE")), label = "p.signif")

ggsave(file.path(out_dir, "CAF_markers_boxplot.pdf"), p_markers, width = 12, height = 10)
message("Saved CAF markers boxplot.")

# 2. Calculate and Plot a General Stromal Signature Score

# Integrate ESTIMATE results
# Align res.deconv with common_samples
res.deconv_aligned <- res.deconv[, common_samples]

clinical_data$ESTIMATE_Stromal_Score <- as.numeric(res.deconv_aligned["StromalScore", ])
clinical_data$ESTIMATE_Immune_Score <- as.numeric(res.deconv_aligned["ImmuneScore", ])
clinical_data$ESTIMATE_Tumor_Purity <- as.numeric(res.deconv_aligned["TumorPurity", ])

# Visualize ESTIMATE Stromal scores
p_estimate <- ggplot(clinical_data, aes(x = cluster, y = ESTIMATE_Stromal_Score, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  theme_bw() +
  scale_fill_manual(values = color.three.clusters) +
  labs(title = "ESTIMATE Stromal Score", y = "ESTIMATE Stromal Score", x = "Subtype") +
  stat_compare_means(comparisons = list(c("TCE", "TASCE"), c("EPCE", "TASCE")), label = "p.signif")

ggsave(file.path(out_dir, "ESTIMATE_Stromal_Score_boxplot.pdf"), p_estimate, width = 6, height = 6)
message("Saved ESTIMATE Stromal Score boxplot.")

# Visualize ESTIMATE Tumor Purity scores
p_tumor_purity <- ggplot(clinical_data, aes(x = cluster, y = ESTIMATE_Tumor_Purity, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
  theme_bw() +
  scale_fill_manual(values = color.three.clusters) +
  labs(title = "ESTIMATE Tumor Purity", y = "ESTIMATE Tumor Purity", x = "Subtype") +
  stat_compare_means(comparisons = list(c("TCE", "TASCE"), c("EPCE", "TASCE")), label = "p.signif")

ggsave(file.path(out_dir, "ESTIMATE_Tumor_Purity_boxplot.pdf"), p_tumor_purity, width = 6, height = 6)
message("Saved ESTIMATE Tumor Purity boxplot.")

message("Validation analysis complete. Results saved to ", out_dir)
