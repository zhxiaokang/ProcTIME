# deconvolute the TIME for dkfz dataset

rm(list = ls())

library(tibble)
library(immunedeconv)
library(dplyr)

load("../data/prad_data_dkfz.RData")

res.deconv.dkfz <- deconvolute(df.gene.dkfz.tpm, "epic", tumor = TRUE)

res.deconv.dkfz <- column_to_rownames(res.deconv.dkfz, colnames(res.deconv.dkfz)[1])

res.deconv.sample <- as.data.frame(t(res.deconv.dkfz))

time.epic.bcr <- merge(res.deconv.sample, df.bcr.dkfz[, c("TIME_FROM_SURGERY_TO_BCR_LASTFU", "BCR_STATUS")], by = "row.names")
time.epic.bcr <- column_to_rownames(time.epic.bcr, colnames(time.epic.bcr)[1])

write.csv(time.epic.bcr, "../data/time_epic_bcr_dkfz.csv")

patients.bcr <- rownames(time.epic.bcr)

res.deconv.dkfz <- res.deconv.dkfz[, patients.bcr]
df.gene.dkfz.tpm <- df.gene.dkfz.tpm[, patients.bcr]
df.bcr.dkfz <- df.bcr.dkfz[patients.bcr, ]

save(res.deconv.dkfz, df.gene.dkfz.tpm, df.bcr.dkfz, time.epic.bcr, file = "../data/dkfz_time_deconv.RData")

