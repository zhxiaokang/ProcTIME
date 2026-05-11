# study the cell types composition between/among different predicted classes using heatmap

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
# prediction labels define the heatmap column split
load("../data/prediction_3_ouh.RData")
cluster <- df.pred$pred

# deconvolution results provide the heatmap matrix
load("../data/ouh_time_deconv.RData")

# ==== draw Heatmap ====

out.dir <- "../output/heatmap/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

cell.types <- row.names(res.deconv.ouh.scale)
cell.types[8] <- "Epithelial cell"
row.names(res.deconv.ouh.scale) <- cell.types

pdf(file.path(out.dir, "heatmap_ouh_3_subtypes.pdf"), width = 8, height = 3)
p <- Heatmap(res.deconv.ouh.scale, column_title = "Cell type proportions across subtypes on the Oslo dataset", name = "Scaled \nproportion",
             column_split = cluster, 
             row_order = c("T cell CD4+", "T cell CD8+", "Epithelial cell", "Cancer associated fibroblast", "Endothelial cell", "Macrophage", "B cell", "NK cell"),
             cluster_columns = FALSE, cluster_rows = FALSE,
             top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color.three.clusters), labels = levels(cluster))),
             show_column_names = FALSE)
print(p)
dev.off()
