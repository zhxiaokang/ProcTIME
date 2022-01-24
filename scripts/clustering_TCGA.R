# Pick the best clustering set

rm(list=ls())

library(dplyr)
library(ConsensusClusterPlus)
library(tibble)
library(survival)
library(survminer)

load("../data/tcga_time_deconv.RData")
# Data included: res.deconv.epic, df.clinical.bcr, prad.gene

# Pick the (best) deconv set
res.deconv <- res.deconv.epic

# Scale: z-score transformation --- 0 mean 1 deviation
res.deconv.scale <- t(scale(t(res.deconv)))

# Clustering
results <- ConsensusClusterPlus(as.matrix(res.deconv.scale), maxK=6, reps=1000, pItem=0.8, pFeature=1, title="../output/Consensus_clustering_epic_pam_pearson_average_tcga", 
                                clusterAlg="pam", distance="pearson", seed=123456, plot="pdf")

# ======== Survival analysis with 3 clusters (but merge 2 clusters) =============
res.deconv.sample <- as.data.frame(t(res.deconv.scale))
res.deconv.sample$cluster <- unname(results[[3]][["consensusClass"]])

# Merge cluster 1 and 2
cluster <- unname(results[[3]][["consensusClass"]])
cluster[which(cluster == "1" | cluster == "2")] <- "good"
cluster[which(cluster == "3")] <- "bad"

res.deconv.sample$cluster <- cluster

df.clinical.bcr.cluster <- merge(df.clinical.bcr, res.deconv.sample, by = "row.names")
df.clinical.bcr.cluster <- column_to_rownames(df.clinical.bcr.cluster, colnames(df.clinical.bcr.cluster)[1])

df.clinical.bcr.cluster$bcr_days <- as.numeric(df.clinical.bcr.cluster$bcr_days)
df.clinical.bcr.cluster$biochemical_recurrence <- 1*(df.clinical.bcr.cluster$biochemical_recurrence=="YES")

# options(repr.plot.width=6, repr.plot.height=5.5)

km.fit <- survfit(Surv(bcr_days, biochemical_recurrence) ~ cluster, data=df.clinical.bcr.cluster)
p.surv <- ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE)

pdf(file.path("../output/Consensus_clustering_epic_pam_pearson_average_tcga", "surv_curv_2_groups.pdf"))
print(p.surv)
dev.off()

save(results, res.deconv, res.deconv.scale, res.deconv.sample, file = "../data/clustering_tcga.RData")
