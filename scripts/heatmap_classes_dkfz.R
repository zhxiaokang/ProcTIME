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
load("../data/prediction_3_dkfz.RData")
cluster <- df.pred$pred

# ==== draw Heatmap ====

out.dir <- "../output/heatmap/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

cell.types <- row.names(res.deconv.dkfz.scale)
cell.types[8] <- "Epithelial cell"
row.names(res.deconv.dkfz.scale) <- cell.types

pdf("../output/heatmap/heatmap_dkfz_3_subtypes.pdf", width = 8, height = 3)
p <- Heatmap(res.deconv.dkfz.scale, column_title = "Cell type proportions across subtypes on the Hamburg dataset", name = "Scaled\nproportion",
        column_split = cluster,
        row_order = c("T cell CD4+", "T cell CD8+", "Epithelial cell", "Cancer associated fibroblast", "Endothelial cell", "Macrophage", "B cell", "NK cell"),
        cluster_columns = FALSE, cluster_rows = FALSE,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color.three.clusters), labels = levels(cluster))), show_column_names = FALSE)
print(p)
dev.off()
