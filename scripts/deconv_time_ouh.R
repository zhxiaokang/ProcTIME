# deconvolute own data using EPIC

rm(list=ls())

library(immunedeconv)
library(tibble)

load("../data/prad_data_ouh.RData")

res.deconv.ouh <- deconvolute(df.gene.ouh.tpm, "epic", tumor = TRUE)
res.deconv.ouh <- column_to_rownames(res.deconv.ouh, names(res.deconv.ouh)[1])

df.sample.deconv <- data.frame(t(res.deconv.ouh))

time.epic.bcr <- merge(df.sample.deconv, df.bcr.ouh[, c("Event", "DaysToEvent")], by = "row.names")
time.epic.bcr <- column_to_rownames(time.epic.bcr, colnames(time.epic.bcr)[1])

write.csv(time.epic.bcr, "../data/time_epic_bcr_ouh.csv")

save(df.sample.deconv, res.deconv.ouh, df.bcr.ouh, time.epic.bcr, file = "../data/ouh_time_deconv.RData")
