# Pick the best clustering set

rm(list=ls())

library(dplyr)
library(ConsensusClusterPlus)
library(tibble)
library(survival)
library(survminer)

# load("../data/tcga_time_deconv.RData")
# Data included: res.deconv.epic, df.clinical.bcr, prad.gene, time.epic.bcr

load("../data/tcga_time_deconv_rm_outlier.RData")
# Data included: res.deconv.epic, df.clinical.bcr, prad.gene, time.epic.bcr, patient.outlier

# Pick the (best) deconv set
res.deconv <- res.deconv.epic

# Scale: z-score transformation --- 0 mean 1 deviation
# res.deconv.scale <- t(scale(t(res.deconv)))
# res.deconv.scale <- res.deconv

# Calculate mean and sd for each feature (row)
res.deconv.t <- t(res.deconv)
res.mean <- colMeans(res.deconv.t)
res.sd <- apply(res.deconv.t, 2, sd)

# Save parameters
save(res.mean, res.sd, file = "../data/tcga_scaling_params.RData")

# Apply scaling
res.deconv.scale <- t(scale(res.deconv.t, center = res.mean, scale = res.sd))


# Clustering
results.tcga <- ConsensusClusterPlus(res.deconv.scale, maxK=6, reps=1000, pItem=0.8, pFeature=1, title="../output/Consensus_clustering_epic_pam_pearson_average_tcga", 
                                     clusterAlg="pam", distance="pearson", seed=123456, plot="pdf")

res.deconv.sample <- as.data.frame(t(res.deconv.scale))

# add "bcr_months" column
df.clinical.bcr$bcr_months <- df.clinical.bcr$bcr_days/30

# ======== Survival analysis with 3 clusters =============
cluster <- unname(results.tcga[[3]][["consensusClass"]])

cluster[which(cluster == "1")] <- "good"
cluster[which(cluster == "2")] <- "mediate"
cluster[which(cluster == "3")] <- "bad"

cluster <- factor(cluster, levels = c("good", "mediate", "bad"))

change_cluster_name <- function(old.names, old.levels = c("good", "mediate", "bad"), new.levels = c("TCE", "TuCE", "TASCE")) {
  new.names <- as.character(old.names)
  
  for (i in seq(1, length(old.levels))){
    new.names[which(new.names == old.levels[i])] <- new.levels[i]
  }
  
  new.names <- factor(new.names, levels = new.levels)
  
  return(new.names)
}

new.levels <- c("TCE", "EPCE", "TASCE")

cluster <- change_cluster_name(cluster, old.levels = c("good", "mediate", "bad"), new.levels)

res.deconv.sample$cluster <- cluster

df.clinical.bcr.cluster <- merge(df.clinical.bcr, res.deconv.sample, by = "row.names")
df.clinical.bcr.cluster <- column_to_rownames(df.clinical.bcr.cluster, colnames(df.clinical.bcr.cluster)[1])

df.clinical.bcr.cluster$bcr_months <- as.numeric(df.clinical.bcr.cluster$bcr_months)
df.clinical.bcr.cluster$biochemical_recurrence <- 1*(df.clinical.bcr.cluster$biochemical_recurrence=="YES")


# load color palette
source("color_palette.R")

km.fit <- survminer::surv_fit(Surv(bcr_months, biochemical_recurrence) ~ cluster, data=df.clinical.bcr.cluster)
p.surv <- survminer::ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE,
                                title = "The TCGA training dataset",
                                palette = color.three.clusters,
                                legend.title = "", legend.labs = c("TCE", "EPCE", "TASCE"),
                                xlab = "Time in months", ylab = "BCR free probability")

pdf(file.path("../output/Consensus_clustering_epic_pam_pearson_average_tcga", "surv_curv_3_groups.pdf"))
print(p.surv)
dev.off()

pdf("../output/surv_curv_3_groups.pdf")
print(p.surv)
dev.off()

res.deconv.tcga <- res.deconv
res.deconv.scale.tcga <- res.deconv.scale
res.deconv.sample.tcga <- res.deconv.sample

save(res.deconv.tcga, res.deconv.scale.tcga, res.deconv.sample.tcga, df.clinical.bcr.cluster, file = "../data/clustering_tcga_original.RData")

# ====== pick out the "bad guys" from "good" cluster ======

# df.good <- dplyr::select(filter(df.clinical.bcr.cluster, cluster == "good"), c(bcr_days, biochemical_recurrence, cluster))
# 
# df.good <- df.good[order(-df.good$bcr_days), ]
# 
# filter(df.good, biochemical_recurrence == "1" & bcr_days > 1000)
# #                 bcr_days biochemical_recurrence cluster
# # TCGA.YL.A9WX     1506                      1       1
# # TCGA.VP.A87D     1194                      1       1

patients.outlier <- rownames(df.clinical.bcr.cluster %>% filter(cluster == new.levels[1] & biochemical_recurrence == "1" & bcr_days > 1000))
print(patients.outlier)
# "TCGA.VP.A87D" "TCGA.YL.A9WX"

# # Also another outlier with extremely different cell types proportions
# patient.extreme <- rownames(filter(time.epic.bcr, `uncharacterized cell` < 0.1))
# print(patient.extreme)

# patients.outlier <- c(patients.outlier, patient.extreme)

df.clinical.bcr.cluster <- df.clinical.bcr.cluster[!(rownames(df.clinical.bcr.cluster) %in% patients.outlier), ]

patients.keep <- rownames(df.clinical.bcr.cluster)

# update the other variables' patients
res.deconv.tcga <- res.deconv.tcga[, patients.keep]
res.deconv.scale.tcga <- res.deconv.scale.tcga[, patients.keep]
res.deconv.sample.tcga <- res.deconv.sample.tcga[patients.keep, ]

km.fit <- survminer::surv_fit(Surv(bcr_months, biochemical_recurrence) ~ cluster, data=df.clinical.bcr.cluster)
p.surv <- survminer::ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE,
                                title = "The TCGA training dataset",
                                palette = color.three.clusters,
                                legend.title = "", legend.labs = c("TCE", "EPCE", "TASCE"),
                                xlab = "Time in months", ylab = "BCR free probability")
                                
pdf("../output/surv_curv_3_groups.pdf")
print(p.surv)
dev.off()

# save the data
save(res.deconv.tcga, res.deconv.scale.tcga, res.deconv.sample.tcga, df.clinical.bcr.cluster, file = "../data/clustering_tcga_removing_badguys.RData")

