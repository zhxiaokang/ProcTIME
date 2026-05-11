# classification using TIME composition from EPIC. Train on TCGA and test on long

library(tibble)
library(dplyr)
library(survminer)
library(ComplexHeatmap)
library(MASS)
library(reshape)
library(sva)

rm(list = ls())
set.seed(12345)

# Load the OUH deconvolution matrix and clinical mapping used as the prediction target cohort.
load("../data/ouh_time_deconv.RData")
# Load the TCGA deconvolution reference because its subtype labels and feature layout define the training space.
load("../data/clustering_tcga_removing_badguys.RData")

# load color palette
source("color_palette.R")

new.levels <- c("TCE", "EPCE", "TASCE")

# Use ComBat to remove batch effect
# TCGA as reference batch
combined.dat <- cbind(res.deconv.tcga, res.deconv.ouh)
batch <- c(rep("TCGA", ncol(res.deconv.tcga)), rep("OUH", ncol(res.deconv.ouh)))

# ComBat with reference batch
adjusted.dat <- ComBat(dat = combined.dat, batch = batch, mod = NULL, par.prior = TRUE, ref.batch = "TCGA")

# Extract OUH
res.deconv.ouh.adj <- adjusted.dat[, (ncol(res.deconv.tcga) + 1):ncol(adjusted.dat)]

# Apply the mean and standard deviation learned from TCGA so OUH samples are transformed with the same scaling as the training data.
load("../data/tcga_scaling_params.RData")
# includes res.mean, res.sd which were calculated from TCGA data

# Scale
res.deconv.ouh.adj.t <- t(res.deconv.ouh.adj)
res.deconv.ouh.scale <- t(scale(res.deconv.ouh.adj.t, center = res.mean, scale = res.sd))

# ======= make predictions on the lda model of 3 clusters =======
load("../data/model_3classes_tcga.RData")

model <- model.lda

pred <- predict(model, t(res.deconv.ouh.scale))$class
pred <- factor(pred, levels = new.levels)

df.pred <- data.frame("pred" = pred)
row.names(df.pred) <- colnames(res.deconv.ouh.scale)

df.patient.pred <- merge(df.pred,
                         dplyr::select(df.bcr.ouh,
                                       c(patient_number, localized.in.focus.number,
                                         gleason.for.biopsy)),
                         by = "row.names")
df.patient.pred$patient_number <- as.factor(df.patient.pred$patient_number)

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
  xlab("Patient ID") + ylab("Predicted cell type\ncomposition subtype") + labs(fill = "Focus")
print(p)
dev.off()

df.pred.bcr <- merge(df.pred, dplyr::select(df.bcr.ouh, c(patient_number, Event, DaysToEvent)),
                     by = "row.names")
df.pred.bcr <- df.pred.bcr[order(df.pred.bcr$patient_number), ]

df.pred.csv <- df.patient.pred %>%
  dplyr::rename(sample_id = Row.names,
                focus = localized.in.focus.number,
                gleason_biopsy = gleason.for.biopsy) %>%
  dplyr::arrange(as.numeric(as.character(patient_number)), focus, sample_id) %>%
  dplyr::select(sample_id, patient_number, focus, gleason_biopsy, pred)
write.csv(df.pred.csv, "../output/prediction_3_ouh.csv", row.names = FALSE)

save(model, res.deconv.ouh.scale, df.pred, file = "../data/prediction_3_ouh.RData")
