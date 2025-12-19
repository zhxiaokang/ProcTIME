# classification using TIME composition from EPIC. Train on TCGA and test on long

library(tibble)
library(dplyr)
library(survival)
library(survminer)
library(MASS)
library(ggpubr)
library(sva)

rm(list = ls())

load("../data/long_time_deconv.RData")  # includes res.deconv.long, df.gene.long.tpm, df.bcr.long, time.epic.bcr
load("../data/model_3classes_tcga.RData")  # includes model.lda
load("../data/tcga_scaling_params.RData")  # includes res.mean, res.sd which were calculated from TCGA data
load("../data/clustering_tcga_removing_badguys.RData")  # includes res.deconv.tcga, res.deconv.scale.tcga, res.deconv.sample.tcga, df.clinical.bcr.cluster

# load color palette
source("color_palette.R")

new.levels <- c("TCE", "EPCE", "TASCE")

surv_long <- function(df.pred, df.bcr) {
  # Survival analysis on the predicted classes
  
  # relevel the 3 clusters
  df.pred$pred <- factor(df.pred$pred, levels = new.levels)
  
  df.bcr.class <- merge(df.bcr, df.pred, by = "row.names")
  df.bcr.class <- column_to_rownames(df.bcr.class, colnames(df.bcr.class)[1])
  
  df.bcr.class$bcr_days <- as.numeric(df.bcr.class$BCR_Months)
  df.bcr.class$bcr_event <- 1*(df.bcr.class$BCR=="1")
  
  df.surv <- dplyr::select(df.bcr.class, c(bcr_days, bcr_event, pred))
  
  options(repr.plot.width=6, repr.plot.height=5.5)
  
  # color.three.clusters <- c("TCE" = "#E64B35", "EPCE" = "#4BBAD3", "TASCE" = "#B3B3B3")
  
  km.fit <- survminer::surv_fit(Surv(bcr_days, bcr_event) ~ pred, data = df.surv)
  
  # # Filter the data to exclude time points with numbers at risk below 5
  # cut.fit <- filter(surv_summary(km.fit), n.risk>=5)
  
  # p.surv <- ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE,
  #                                 title = "The Atlanta validation dataset",
  #                                 palette = color.three.clusters,
  #                                 legend.title = "", legend.labs = c("TCE", "EPCE", "TASCE"),
  #                                 xlab = "Time in months", ylab = "BCR free probability")
  
  # to keep only 2 digits of p-value
  p.surv <- ggsurvplot(km.fit, pval = 0.04, risk.table = TRUE, ncensor.plot = FALSE,
                       title = "The Atlanta validation dataset",
                       palette = color.three.clusters,
                       legend.title = "", legend.labs = c("TCE", "EPCE", "TASCE"),
                       xlab = "Time in months", ylab = "BCR free probability")
  return(p.surv)
}

# Scale the long dataset using TCGA scaling parameters

# Use ComBat to remove batch effect
# TCGA as reference batch
combined.dat <- cbind(res.deconv.tcga, res.deconv.long)
batch <- c(rep("TCGA", ncol(res.deconv.tcga)), rep("Long", ncol(res.deconv.long)))

# ComBat with reference batch
adjusted.dat <- ComBat(dat = combined.dat, batch = batch, mod = NULL, par.prior = TRUE, ref.batch = "TCGA")

# Extract Long
res.deconv.long.adj <- adjusted.dat[, (ncol(res.deconv.tcga) + 1):ncol(adjusted.dat)]

# Apply TCGA scaling parameters

# Scale
res.deconv.long.adj.t <- t(res.deconv.long.adj)
res.deconv.long.scale <- t(scale(res.deconv.long.adj.t, center = res.mean, scale = res.sd))

# ======= make predictions on the lda model and draw surv curv based on the pred =======
model <- model.lda

pred <- predict(model, t(res.deconv.long.scale))$class
pred <- factor(pred, levels = new.levels)

df.pred <- data.frame("pred" = pred)
row.names(df.pred) <- colnames(res.deconv.long.scale)

save(model, res.deconv.long.scale, df.pred, file = "../data/prediction_3_long.RData")

pdf("../output/classification_3classes_long.pdf")
p <- surv_long(df.pred, df.bcr.long)
print(p)
dev.off()

