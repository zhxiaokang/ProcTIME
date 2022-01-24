# deconvolute the TIME for MSKCC dataset

rm(list = ls())

library(tibble)
library(immunedeconv)
library(dplyr)

load("../data/prad_data_mskcc.RData")

# convert FPKM to TPM
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

df.gene.mskcc.tpm <- t(apply(df.gene.mskcc, 1, fpkmToTpm))

res.deconv.mskcc <- deconvolute(as.matrix(df.gene.mskcc.tpm), "epic", tumor = TRUE)

res.deconv.mskcc <- column_to_rownames(res.deconv.mskcc, colnames(res.deconv.mskcc)[1])

res.deconv.sample <- as.data.frame(t(res.deconv.mskcc))

time.epic.bcr <- merge(res.deconv.sample, df.bcr.mskcc[, c("BCR_FreeTime", "BCR_Event")], by = "row.names")
time.epic.bcr <- column_to_rownames(time.epic.bcr, colnames(time.epic.bcr)[1])

time.epic.bcr <- filter(time.epic.bcr, BCR_Event != "NA")

write.csv(time.epic.bcr, "../data/time_epic_bcr_mskcc.csv")

patients.bcr <- rownames(time.epic.bcr)

res.deconv.mskcc <- res.deconv.mskcc[, patients.bcr]
df.gene.mskcc.tpm <- df.gene.mskcc.tpm[, patients.bcr]
df.bcr.mskcc <- df.bcr.mskcc[patients.bcr, ]

save(res.deconv.mskcc, df.gene.mskcc.tpm, df.bcr.mskcc, time.epic.bcr, file = "../data/mskcc_time_deconv.RData")

