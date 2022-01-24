# Clustering on MSKCC

library(dplyr)
library(ConsensusClusterPlus)
library(tibble)
library(survival)
library(survminer)

rm(list = ls())

load("../data/mskcc_time_deconv.RData")
# include: res.deconv.mskcc, df.bcr.mskcc, df.gene.mskcc.tpm

# Scale: z-score transformation --- 0 mean 1 deviation
res.deconv.mskcc.scale <- t(scale(t(res.deconv.mskcc)))

results.mskcc <- ConsensusClusterPlus(as.matrix(res.deconv.mskcc.scale), maxK=6, reps=1000, pItem=0.8, pFeature=1, title="../output/Consensus_clustering_epic_pam_pearson_average_mskcc", 
                                clusterAlg="pam", distance="pearson", innerLinkage="average", seed=123456, plot="pdf")

res.deconv.mskcc.scale.sample <- as.data.frame(t(res.deconv.mskcc.scale))

# ===== draw KM plot of clusters =====
cluster <- unname(results.mskcc[[3]][["consensusClass"]])

res.deconv.mskcc.scale.sample$class <- cluster

df.bcr.mskcc.class <- merge(res.deconv.mskcc.scale.sample, df.bcr.mskcc, by = "row.names")
df.bcr.mskcc.class <- column_to_rownames(df.bcr.mskcc.class, colnames(df.bcr.mskcc.class)[1])

df.bcr.mskcc.class$bcr_days <- as.numeric(df.bcr.mskcc.class$BCR_FreeTime)
df.bcr.mskcc.class$biochemical_recurrence <- 1*(df.bcr.mskcc.class$BCR_Event=="BCR_Algorithm")

df.surv <- dplyr::select(df.bcr.mskcc.class, c(bcr_days, biochemical_recurrence, class))

df.surv <- na.omit(df.surv)

options(repr.plot.width=6, repr.plot.height=5.5)

km.fit <- survfit(Surv(bcr_days, biochemical_recurrence) ~ class, data=df.surv)
p.surv <- ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE)

pdf(file.path("../output/Consensus_clustering_epic_pam_pearson_average_mskcc", "surv_curv_3_groups.pdf"))
print(p.surv)
dev.off()

# Merge cluster 2 and 3
cluster <- unname(results.mskcc[[3]][["consensusClass"]])
cluster[which(cluster == "2" | cluster == "3")] <- "good"
cluster[which(cluster == "1")] <- "bad"

res.deconv.mskcc.scale.sample$class <- cluster

df.bcr.mskcc.class <- merge(res.deconv.mskcc.scale.sample, df.bcr.mskcc, by = "row.names")
df.bcr.mskcc.class <- column_to_rownames(df.bcr.mskcc.class, colnames(df.bcr.mskcc.class)[1])

df.bcr.mskcc.class$bcr_days <- as.numeric(df.bcr.mskcc.class$BCR_FreeTime)
df.bcr.mskcc.class$biochemical_recurrence <- 1*(df.bcr.mskcc.class$BCR_Event=="BCR_Algorithm")

df.surv <- dplyr::select(df.bcr.mskcc.class, c(bcr_days, biochemical_recurrence, class))

df.surv <- na.omit(df.surv)

options(repr.plot.width=6, repr.plot.height=5.5)

km.fit <- survfit(Surv(bcr_days, biochemical_recurrence) ~ class, data=df.surv)
p.surv <- ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE)

pdf(file.path("../output/Consensus_clustering_epic_pam_pearson_average_mskcc", "surv_curv_2_groups.pdf"))
print(p.surv)
dev.off()

save(results.mskcc, res.deconv.mskcc, res.deconv.mskcc.scale, res.deconv.mskcc.scale.sample, file = "../data/clustering_mskcc.RData")
