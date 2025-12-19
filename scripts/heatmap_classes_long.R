# study the cell types composition between/among different predicted classes using heatmap

rm(list = ls())

# library
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)
# library(cowplot)
library(gridGraphics)
library(grid)

# load color palette
source("color_palette.R")

# ==== load datasets ====
load("../data/prediction_3_long.RData")
cluster <- df.pred$pred

# ==== draw Heatmap ====

out.dir <- "../output/heatmap/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

# # function scale to [0, 1]
# scale_zero_one <- function(vec) {
#   min.vec <- min(vec)
#   max.vec <- max(vec)
#   
#   vec.scale <- (vec - min.vec) / (max.vec - min.vec)
#   return(vec.scale)
# }
# 
# res.deconv.tcga.zero.one <- t(apply(res.deconv.tcga, 1, scale_zero_one))

cell.types <- row.names(res.deconv.long.scale)
cell.types[8] <- "Epithelial cell"
row.names(res.deconv.long.scale) <- cell.types

pdf("../output/heatmap/heatmap_long_3_subtypes.pdf", width = 8, height = 3)
p <- Heatmap(res.deconv.long.scale, column_title = "Cell type proportions across subtypes on the Atlanta dataset", name = "Scaled\nproportion",
        column_split = cluster, 
        row_order = c("T cell CD4+", "T cell CD8+", "Epithelial cell", "Cancer associated fibroblast", "Endothelial cell", "Macrophage", "B cell", "NK cell"),
        cluster_columns = FALSE, cluster_rows = FALSE,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color.three.clusters), labels = levels(cluster))), show_column_names = FALSE)
print(p)
dev.off()
# 
# # 2 clusters
# load("../data/prediction_2_long.RData")
# cluster <- df.pred$pred
# 
# pdf("../output/heatmap/heatmap_long_2clusters.pdf", width = 7, height = 3)
# p <- Heatmap(res.deconv.long.scale, column_title = "Cell type proportions across clusters on Dataset Long", name = "Scaled \nproportion",
#         column_split = cluster, 
#         cluster_columns = FALSE, cluster_rows = TRUE,
#         top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color.two.clusters), labels = levels(cluster))), show_column_names = FALSE)
# print(p)
# dev.off()
# 
# pdf("../output/heatmap/heatmap_long_2clusters_with_dendrogram.pdf", width = 7, height = 3)
# p <- Heatmap(res.deconv.long.scale, column_title = "Cell type proportions across clusters on Dataset Long", name = "Scaled \nproportion",
#         column_split = cluster, 
#         cluster_columns = TRUE, cluster_rows = TRUE,
#         top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color.two.clusters), labels = levels(cluster))), show_column_names = FALSE)
# print(p)
# dev.off()
