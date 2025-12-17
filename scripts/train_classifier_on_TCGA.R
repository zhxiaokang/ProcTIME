# Train a classifier of 3 classes from TCGA's clustering results

library(dplyr)
library(MASS)  # for classification method LDA: Linear discriminant analysis

rm(list = ls())

load("../data/clustering_tcga_removing_badguys.RData")
# incluse: "res.deconv.sample.tcga" "res.deconv.scale.tcga"  "res.deconv.tcga"

# ======= lda on 3 classes =======
class.label <- res.deconv.sample.tcga$cluster
x <- as.matrix(dplyr::select(res.deconv.sample.tcga, -c(cluster)))
model.lda <- lda(x, grouping = class.label)

save(model.lda, file = "../data/model_3classes_tcga.RData")

# # ====== merge "mediate" and "bad" into "bad" ======
# cluster <- res.deconv.sample.tcga$cluster
# cluster[which(cluster == "mediate" | cluster == "bad")] <- "bad"
# res.deconv.sample.tcga$cluster <- cluster
# 
# class.label <- res.deconv.sample.tcga$cluster
# x <- as.matrix(dplyr::select(res.deconv.sample.tcga, -c(cluster)))
# model.lda <- lda(x, grouping = class.label)
# 
# save(model.lda, file = "../data/model_2classes_tcga.RData")

