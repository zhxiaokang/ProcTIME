# study the cell types composition between/among different clusters using heatmap

rm(list = ls())

# library
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)
library(cowplot)
library(gridGraphics)
library(grid)
library(circlize)

# ==== TCGA ====
# load color palette
source("color_palette.R")

load("../data/clustering_tcga_original.RData")

cell.types <- rownames(res.deconv.scale.tcga)
print(cell.types)

cell.types.short <- c("B cell", "Cancer associated fibroblast", "T cell CD4+", "T cell CD8+", "Endothelial cell", "Macrophage", "NK cell", "Epithelial cell")
rownames(res.deconv.scale.tcga) <- cell.types.short

pdf("../output/heatmap/heatmap_tcga_deconv.pdf", width = 8, height = 3)
p <- Heatmap(res.deconv.scale.tcga, column_title = NULL, show_heatmap_legend = FALSE,
             row_order = c("T cell CD4+", "T cell CD8+", "Epithelial cell", "Cancer associated fibroblast", "Endothelial cell", "Macrophage", "B cell", "NK cell"),
             cluster_columns = FALSE, cluster_rows = FALSE,
             show_column_names = FALSE)
dev.off()

pdf("../output/heatmap/heatmap_tcga_deconv_cluster.pdf", width = 8, height = 3.2)
p <- Heatmap(res.deconv.scale.tcga, column_title = "Cell type proportions across clusters on the TCGA dataset",
             show_heatmap_legend = TRUE, name = "Scaled\nproportion",
             column_split = res.deconv.sample.tcga$cluster,
             row_order = c("T cell CD4+", "T cell CD8+", "Epithelial cell", "Cancer associated fibroblast", "Endothelial cell", "Macrophage", "B cell", "NK cell"),
             cluster_columns = FALSE, cluster_rows = FALSE,
             top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color.three.clusters), labels = c("TCE", "EPCE", "TASCE"))),
             show_column_names = FALSE)
print(p)
dev.off()

# ==== DKFZ ====
rm(list = ls())

# load color palette
source("color_palette.R")

load("../data/clustering_dkfz.RData")

cell.types <- rownames(res.deconv.dkfz.scale)
print(cell.types)

cell.types.short <- c("B cell", "Cancer associated fibroblast", "T cell CD4+", "T cell CD8+", "Endothelial cell", "Macrophage", "NK cell", "Epithelial cell")
rownames(res.deconv.dkfz.scale) <- cell.types.short

pdf("../output/heatmap/heatmap_dkfz_deconv.pdf", width = 5, height = 3)
p <- Heatmap(res.deconv.dkfz.scale, column_title = NULL, show_heatmap_legend = FALSE,
             row_order = c("T cell CD4+", "T cell CD8+", "Epithelial cell", "Cancer associated fibroblast", "Endothelial cell", "Macrophage", "B cell", "NK cell"),
             cluster_columns = FALSE, cluster_rows = FALSE,
             show_column_names = FALSE)
print(p)
dev.off()

# ==== Long ====
rm(list = ls())

# load color palette
source("color_palette.R")

load("../data/clustering_long.RData")

cell.types <- rownames(res.deconv.long.scale)
print(cell.types)

cell.types.short <- c("B cell", "Cancer associated fibroblast", "T cell CD4+", "T cell CD8+", "Endothelial cell", "Macrophage", "NK cell", "Epithelial cell")
rownames(res.deconv.long.scale) <- cell.types.short

pdf("../output/heatmap/heatmap_long_deconv.pdf", width = 5, height = 3)
p <- Heatmap(res.deconv.long.scale, column_title = NULL, show_heatmap_legend = FALSE,
             row_order = c("T cell CD4+", "T cell CD8+", "Epithelial cell", "Cancer associated fibroblast", "Endothelial cell", "Macrophage", "B cell", "NK cell"),
             cluster_columns = FALSE, cluster_rows = FALSE,
             show_column_names = FALSE)
print(p)
dev.off()


