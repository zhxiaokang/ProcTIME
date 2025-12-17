# study the cell types composition between/among different clusters using heatmap

rm(list = ls())

# library
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)
# library(ggpubr)
library(gridGraphics)
library(grid)

# load color palette
source("color_palette.R")

# ==== load datasets ====
load("../data/prediction_3_ouh.RData")
cluster <- df.pred$pred

load("../data/ouh_time_deconv.RData")

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

gs <- df.bcr.ouh$gleason_score_updated
gs.cat <- gs

# `%ni%` <- Negate(`%in%`)

gs.cat[gs.cat %in% c("0+3", "3+0", "3+3", "0+3+3", "0+3+4", "3+4", "3+3+0", "3+4+0", "3+4+5", "3+5+0")] <- "low"
gs.cat[gs.cat %in% c("4+3")] <- "intermediate"
gs.cat[gs.cat %in% c("0+4", "0+4+4", "0+4+5", "4+4", "4+5")] <- "high"

# gs.cat[gs.cat %in% c("0+3", "3+0", "3+3", "0+3+4", "3+4", "3+4+0", "3+4+5")] <- "low"
# gs.cat[gs.cat %in% c("4+3", "0+4", "0+4+4", "4+4", "4+5")] <- "high"

cell.types <- row.names(res.deconv.ouh.scale)
cell.types[8] <- "Epithelial cell"
row.names(res.deconv.ouh.scale) <- cell.types

pdf(file.path(out.dir, "heatmap_ouh_3clusters.pdf"), width = 8, height = 3)
p <- Heatmap(res.deconv.ouh.scale, column_title = "Cell type proportions across clusters on the Oslo dataset", name = "Scaled \nproportion",
             column_split = cluster, 
             row_order = c("T cell CD4+", "T cell CD8+", "Epithelial cell", "Cancer associated fibroblast", "Endothelial cell", "Macrophage", "B cell", "NK cell"),
             cluster_columns = FALSE, cluster_rows = FALSE,
             top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color.three.clusters), labels = levels(cluster))),
             show_column_names = FALSE)
print(p)
dev.off()

# pdf("../output/heatmap/heatmap_ouh_3clusters.pdf", width = 9, height = 3)
# p <- Heatmap(res.deconv.ouh.scale, column_title = "Cell type proportions across clusters on dataset Oslo", name = "Scaled \nproportion",
#              column_split = cluster, 
#              row_order = c("T cell CD4+", "T cell CD8+", "uncharacterized cell", "Cancer associated fibroblast", "Endothelial cell", "Macrophage", "B cell", "NK cell"),
#              cluster_columns = FALSE, cluster_rows = FALSE,
#              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color.three.clusters), labels = levels(cluster))), 
#              bottom_annotation = HeatmapAnnotation(EPE = df.bcr.ouh$epe, GleasonScore = gs.cat, col = list(GleasonScore=c("high"=color.gleason.score[1], "intermediate"=color.gleason.score[2], "low"=color.gleason.score[3]))),
#              show_column_names = FALSE)
# 
# print(p)
# dev.off()

# # 2 clusters
# load("../data/prediction_2_ouh.RData")
# cluster <- df.pred$pred
# 
# pdf("../output/heatmap/heatmap_ouh_2clusters.pdf", width = 7, height = 8)
# p <- Heatmap(res.deconv.ouh.scale, column_title = "Cell type proportions across clusters on Dataset OUH", name = "Scaled \nproportion",
#              column_split = cluster, 
#              cluster_columns = FALSE, cluster_rows = TRUE,
#              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color.two.clusters), labels = levels(cluster))), show_column_names = FALSE)
# print(p)
# dev.off()
# 
# pdf("../output/heatmap/heatmap_ouh_2clusters_with_dendrogram.pdf", width = 7, height = 8)
# p <- Heatmap(res.deconv.ouh.scale, column_title = "Cell type proportions across clusters on Dataset OUH", name = "Scaled \nproportion",
#              column_split = cluster, 
#              cluster_columns = TRUE, cluster_rows = TRUE,
#              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color.two.clusters), labels = levels(cluster))), show_column_names = FALSE)
# print(p)
# dev.off()
