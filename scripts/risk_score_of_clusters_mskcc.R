# combine the clustering and risk score of MSKCC

library(ggplot2)
library(tibble)
library(dplyr)
library(ggpubr)

rm(list = ls())

load("../data/clustering_mskcc.RData")

df.clustering <- select(res.deconv.mskcc.scale.sample, c(class))

rs <- read.csv("../data/risk_score_from_DeepSurv_mskcc.csv")

rs <- column_to_rownames(rs, names(rs)[1])

df.class.rs <- merge(df.clustering, rs, by = "row.names")
df.class.rs <- column_to_rownames(df.class.rs, names(df.class.rs)[1])
df.class.rs$class <- as.factor(df.class.rs$class)

p.violin <- ggplot(df.class.rs, aes(x = class, y = risk.score, fill = class)) +
  geom_violin(trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6, aes(fill = class)) +
  stat_compare_means(method = "wilcox.test",
                     label.x = 1.2, 
                     label.y = 0.055) +
  labs(x = "Cluster", y = "Risk Score")

pdf("../output/risk_score_of_clusters_mskcc.pdf")
print(p.violin)
dev.off()

