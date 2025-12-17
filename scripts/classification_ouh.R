# classification using TIME composition from EPIC. Train on TCGA and test on long

library(tibble)
library(dplyr)
library(survminer)
library(ComplexHeatmap)
library(MASS)
library(reshape)

rm(list = ls())

load("../data/ouh_time_deconv.RData")

# load color palette
source("color_palette.R")

new.levels <- c("TCE", "EPCE", "TASCE")

# scale the cell types' proportions
## z-score transformation --- 0 mean 1 deviation
res.deconv.ouh.scale <- t(scale(t(res.deconv.ouh)))

# # ======= make predictions on the lda model of 2 clusters =======
# load("../data/model_2classes_tcga.RData")
# 
# model <- model.lda
# 
# pred <- predict(model, t(res.deconv.ouh.scale))$class
# pred <- factor(pred, levels = c("good", "bad"))
# 
# df.pred <- data.frame("pred" = pred)
# row.names(df.pred) <- colnames(res.deconv.ouh.scale)
# 
# df.patient.pred <- merge(df.pred, dplyr::select(df.bcr.ouh, c(patient_number, sample_name)), by.x = "row.names", by.y = "sample_name")
# 
# pdf("../output/pred2_of_sample_vs_patient.pdf", width = 16, height = 4.3)
# p <- ggplot(df.patient.pred, aes(x = as.factor(patient_number), y = pred)) + geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.08)
# print(p)
# dev.off()
# # 
# # df.pred.bcr <- merge(df.pred, dplyr::select(df.bcr.ouh, c(patient_number, sample_name, Event, DaysToEvent)), 
# #                      by.x = "row.names", by.y = "sample_name")
# # df.pred.bcr <- df.pred.bcr[order(df.pred.bcr$patient_number), ]
# # 
# # save(model, res.deconv.ouh.scale, df.pred, file = "../data/prediction_2_ouh.RData")

# ======= make predictions on the lda model of 3 clusters =======
load("../data/model_3classes_tcga.RData")

model <- model.lda

pred <- predict(model, t(res.deconv.ouh.scale))$class
pred <- factor(pred, levels = new.levels)

df.pred <- data.frame("pred" = pred)
row.names(df.pred) <- colnames(res.deconv.ouh.scale)

df.patient.pred <- merge(df.pred, dplyr::select(df.bcr.ouh, c(patient_number, localized.in.focus.number)), by = "row.names")
df.patient.pred$patient_number <- as.factor(df.patient.pred$patient_number)

# df.patient.pred <- df.patient.pred %>% dplyr::filter(focus != "F3")

df.rect <- data.frame(xmin = seq(1, 23) - 0.4,
                      xmax = seq(1, 23) + 0.4,
                      ymin = rep(0, 23), 
                      ymax = rep(3, 23) + 0.6)

pdf("../output/pred3_of_sample_vs_patient.pdf", width = 8, height = 2.2)
p <- ggplot() +
  geom_dotplot(data = df.patient.pred, aes(x = patient_number, y = pred, fill = localized.in.focus.number), 
               stackgroups=TRUE, binpositions="all", binaxis = "y", stackdir = "center", binwidth = 0.17) +
  geom_rect(data = df.rect, aes(xmin = xmin, xmax = xmax, 
                                ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.2) +
  theme_classic() +
  xlab("Patient ID") + ylab("Predicted cell type\ncomposition subtype") + labs(fill = "localized.in.focus.number")
print(p)
dev.off()

# ====== use clustering result instead of predicted results ======
# load("../data/clustering_ouh.RData")
# 
# df.bcr.ouh.class$patient_number <- as.factor(df.bcr.ouh.class$patient_number)
# 
# # df.bcr.ouh.class <- df.bcr.ouh.class %>% dplyr::filter(focus != "F3")
# 
# df.rect <- data.frame(xmin = seq(1, 23) - 0.4,
#                       xmax = seq(1, 23) + 0.4,
#                       ymin = rep(0, 23), 
#                       ymax = rep(3, 23) + 0.6)
# 
# pdf("../output/clustering_of_sample_vs_patient.pdf", width = 8, height = 2.2)
# p <- ggplot() +
#   geom_dotplot(data = df.bcr.ouh.class, aes(x = patient_number, y = class, fill = focus), 
#                stackgroups=TRUE, binpositions="all", binaxis = "y", stackdir = "center", binwidth = 0.17) +
#   geom_rect(data = df.rect, aes(xmin = xmin, xmax = xmax, 
#                                 ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.2) +
#   theme_classic() +
#   xlab("Patient ID") + ylab("Predicted TIME subtype") + labs(fill = "Focus")
# print(p)
# dev.off()

df.pred.bcr <- merge(df.pred, dplyr::select(df.bcr.ouh, c(patient_number, Event, DaysToEvent)),
                     by = "row.names")
df.pred.bcr <- df.pred.bcr[order(df.pred.bcr$patient_number), ]

save(model, res.deconv.ouh.scale, df.pred, file = "../data/prediction_3_ouh.RData")
# 
# cluster <- df.pred$pred
# col.ha <- HeatmapAnnotation(patient = as.character(df.bcr.ouh$patient_number))
# Heatmap(res.deconv.ouh.scale, column_title = "Cell type proportions across clusters on OUH", name = "Scaled \nproportion",
#                       column_split = cluster, 
#                       cluster_columns = TRUE, cluster_rows = TRUE,
#                       bottom_annotation = col.ha,
#                       top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color.three.clusters), labels = levels(cluster))), show_column_names = FALSE)


