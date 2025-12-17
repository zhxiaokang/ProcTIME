# pre-process OUH data

library(openxlsx)
library(tibble)
library(dplyr)

rm(list=ls())

# ==== gene expression ====
df.tpm <- read.table("../data/OUH/tpm_batch_corrected_18349_genes_reshaped.tsv", header = TRUE, sep = "\t", row.names = 1)
# remove prad_sample14 since it is duplicate of prad_sample55
df.tpm <- df.tpm[, colnames(df.tpm) != "prad_sample14"]

# ==== clinical info ====
# samples' clinical info
df.sample.clinical <- read.xlsx("../data/OUH/til_xiaokang_data.xlsx", sheet = "updated")
# remove "Benign" samples
df.sample.clinical <- dplyr::filter(df.sample.clinical, localized.in.focus.number != "Benign")
# remove prad_sample14 since it is duplicate of prad_sample55
df.sample.clinical <- dplyr::filter(df.sample.clinical, sample_name != "prad_sample14")

# patients' BCR info
df.patient.bcr <- read.xlsx("../data/OUH/til_xiaokang_bcr.xlsx", sheet = "updated")
# patients' age info
df.patient.age <- read.xlsx("../data/OUH/til_xiaokang_age.xlsx")

# merge BCR and age
df.patient.bcr.age <- merge(df.patient.bcr, df.patient.age, by = "patient_number")

# merge samples and patients
df.sample.clinical.bcr <- merge(df.sample.clinical, df.patient.bcr.age, by = "patient_number")
df.sample.clinical.bcr <- column_to_rownames(df.sample.clinical.bcr, var = "sample_name")

# ==== intersect gene expression and clinical patients/samples ====
samples.gene <- colnames(df.tpm)
samples.clinical <- rownames(df.sample.clinical.bcr)

samples.inter <- intersect(samples.gene, samples.clinical)

# Only keep the patients that have both gene exp and BCR info
df.gene.ouh.tpm <- df.tpm[, samples.inter]
df.bcr.ouh <- df.sample.clinical.bcr[samples.inter, ]

# ==== save the data ====
save(df.gene.ouh.tpm, df.bcr.ouh, file = "../data/prad_data_ouh.RData")

