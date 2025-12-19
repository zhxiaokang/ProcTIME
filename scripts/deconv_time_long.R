# deconvolute the TIME for long dataset

rm(list = ls())

library(tibble)
library(immunedeconv)
library(dplyr)

load("../data/prad_data_long.RData")

res.deconv.long <- deconvolute(df.gene.long.tpm, "epic", tumor = TRUE)

res.deconv.long <- column_to_rownames(res.deconv.long, colnames(res.deconv.long)[1])

res.deconv.sample <- as.data.frame(t(res.deconv.long))

time.epic.bcr <- merge(res.deconv.sample, df.bcr.long[, c("BCR_Months", "BCR")], by = "row.names")
time.epic.bcr <- column_to_rownames(time.epic.bcr, colnames(time.epic.bcr)[1])

write.csv(time.epic.bcr, "../data/time_epic_bcr_long.csv")

patients.bcr <- rownames(time.epic.bcr)

res.deconv.long <- res.deconv.long[, patients.bcr]
df.gene.long.tpm <- df.gene.long.tpm[, patients.bcr]
df.bcr.long <- df.bcr.long[patients.bcr, ]

save(res.deconv.long, df.gene.long.tpm, df.bcr.long, time.epic.bcr, file = "../data/long_time_deconv.RData")

