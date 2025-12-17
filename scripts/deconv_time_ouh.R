# deconvolute own data using EPIC

rm(list=ls())

library(immunedeconv)
library(tibble)

# load("../data/prad_data_ouh_bak.RData")
load("../data/prad_data_ouh.RData")

res.deconv.ouh <- deconvolute(df.gene.ouh.tpm, "epic", tumor = TRUE)
res.deconv.ouh <- column_to_rownames(res.deconv.ouh, names(res.deconv.ouh)[1])

df.sample.deconv <- data.frame(t(res.deconv.ouh))

# time.epic.bcr <- merge(df.sample.deconv, df.bcr.ouh[, c("clinical_outcome", "days_to_event")], by = "row.names")
time.epic.bcr <- merge(df.sample.deconv, df.bcr.ouh[, c("Event", "DaysToEvent")], by = "row.names")
time.epic.bcr <- column_to_rownames(time.epic.bcr, colnames(time.epic.bcr)[1])

write.csv(time.epic.bcr, "../data/time_epic_bcr_ouh.csv")

save(df.sample.deconv, res.deconv.ouh, df.bcr.ouh, time.epic.bcr, file = "../data/ouh_time_deconv.RData")

# # ============= OLD ==============
# df.tpm <- read.table("../data/OUH/tpm_batch_corrected_18349_genes_reshaped.tsv", header = TRUE, sep = "\t", row.names = 1)
# res.deconv.ouh <- deconvolute(df.tpm, "epic")
# res.deconv.ouh <- column_to_rownames(res.deconv.ouh, names(res.deconv.ouh)[1])
# 
# df.sample.deconv <- data.frame(t(res.deconv.ouh))
# 
# write.csv(df.sample.deconv, file = "../data/OUH/time_epic_ouh.csv", row.names = TRUE)
# 
# # clinical info
# df.patient.bcr <- read.csv("../data/OUH/bcr_23_patients.tsv", sep = "\t")
# df.sample.clinical <- read.csv("../data/OUH/88_samples_info_canonical_modify.tsv", sep = "\t")
# 
# df.sample.clinical.bcr <- merge(df.sample.clinical, df.patient.bcr, by = "patient_number")
# 
# save(df.sample.deconv, res.deconv.ouh, df.sample.clinical.bcr, file = "../data/ouh_time_deconv.RData")

