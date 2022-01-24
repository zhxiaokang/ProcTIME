# deconvolute the TIME for TCGA dataset

library(tibble)
library(immunedeconv)

load("../data/prad_data_tcga.RData")

# deconvolution
res.deconv <- deconvolute(as.matrix(prad.gene), "epic", tumor = TRUE)

res.deconv.epic <- column_to_rownames(res.deconv, colnames(res.deconv)[1])

res.deconv.sample <- as.data.frame(t(res.deconv.epic))

time.epic.bcr <- merge(res.deconv.sample, df.clinical.bcr[, c("bcr_days", "biochemical_recurrence")], by = "row.names")
time.epic.bcr <- column_to_rownames(time.epic.bcr, colnames(time.epic.bcr)[1])

write.csv(time.epic.bcr, "../data/time_epic_bcr_tcga.csv")

save(res.deconv.epic, prad.gene, df.clinical.bcr, time.epic.bcr, file = "../data/tcga_time_deconv.RData")
