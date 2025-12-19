# pre-process DKFZ2018 data

rm(list=ls())

library(tidyverse)
library(tibble)
library(dplyr)

# ====== preprocess gene expression info ======
dkfz.gene <- read.csv("../data/prostate_dkfz_2018/data_mrna_seq_rpkm.txt", sep = "\t", na.strings = "")

pick_dup <- function(vec) {
  # pick out the duplicates from vector
  tab.vec <- as.data.frame(table(vec))
  tab.dup <- filter(tab.vec, Freq > 1)
  vec.dup <- as.character(tab.dup[, 1])
  
  return(list(vec.dup, tab.dup))
}

remove_dup <- function(vec) {
  # remove duplicates from vector (the duplicates are all removed, none of them is kept)
  tab.vec <- as.data.frame(table(vec))
  tab.uniq <- filter(tab.vec, Freq == 1)
  vec.uniq <- as.character(tab.uniq[, 1])
  
  return(vec.uniq)
}

list.dup <- pick_dup(dkfz.gene$Hugo_Symbol)

dkfz.gene.dup <- filter(dkfz.gene, Hugo_Symbol %in% list.dup[[1]])
# mostly with very low expression, so just delete them

gene.uniq <- remove_dup(dkfz.gene$Hugo_Symbol)

dkfz.gene.uniq <- filter(dkfz.gene, Hugo_Symbol %in% gene.uniq)
dkfz.gene.uniq <- column_to_rownames(dkfz.gene.uniq, "Hugo_Symbol")
dkfz.gene.uniq <- dplyr::select(dkfz.gene.uniq, -c(Entrez_Gene_Id))

# when there are multiple samples from one patient, pick the first sample from that patient as representative

pick_first_from_multi <- function(str.vec) {
  # pick the first (by index) when there are multiple samples from the same patient
  # str.vec: a vector of strings, where the strings are the sample IDs like "ICGC_PCA192_T01"
  patient.id <- unique(substr(str.vec, 1, 11))
  
  first.samples.patient <- c()
  for(id in patient.id){
    id.of.patient <- grep(id, str.vec, value = TRUE)
    first.sample.patient <- sort(id.of.patient)[1]
    first.samples.patient <- c(first.samples.patient, first.sample.patient)
  }
  
  return(first.samples.patient)
}

sample.id <- colnames(dkfz.gene.uniq)
first.samples.patient <- pick_first_from_multi(sample.id)

dkfz.gene.unifocal <- dkfz.gene.uniq[, first.samples.patient]
# remove the "_T0*" part from the sample/patient IDs

first.samples <- substr(first.samples.patient, 1, 11)
colnames(dkfz.gene.unifocal) <- first.samples

# function to convert FPKM to TPM
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

df.gene.dkfz.tpm <- apply(dkfz.gene.unifocal, 2, fpkmToTpm)

# ====== preprocess clinical info ======
dkfz.clinical <- read.csv("../data/prostate_dkfz_2018/data_clinical_patient.txt", skip = 4, header = TRUE, sep = "\t", na.strings = "")

dkfz.clinical <- column_to_rownames(dkfz.clinical, colnames(dkfz.clinical)[1])

df.bcr.dkfz <- dplyr::filter(dkfz.clinical, BCR_STATUS != "NA")

# the above file misses some info, such as gleason score
df.clinical.complete <- read.csv("../data/prostate_dkfz_2018/prostate_dkfz_2018_clinical_data.tsv", sep = "\t")
# for the patients with several samples, only pick out the 1st one
temp <- df.clinical.complete %>% dplyr::filter(str_detect(Sample.ID, "CPCG") | str_detect(Sample.ID, "T01"))
df.clinical.complete <- column_to_rownames(temp, "Patient.ID")
# remove the patients whose BCR status is unavailable
df.clinical.complete <- dplyr::filter(df.clinical.complete, BCR.Status != "NA")

# since all the downstream analyses were done based on the original df.bcr.dkfz, so still keep it but the other cols from df.clinical.complete shall be added
df.merge <- merge(df.clinical.complete, df.bcr.dkfz, by = "row.names")
df.bcr.dkfz <- column_to_rownames(df.merge, colnames(df.merge)[1])

# ====== intersect the gene expression and clinical patients ======
patient.gene <- colnames(df.gene.dkfz.tpm)

patient.clinical <- rownames(df.bcr.dkfz)

id.inter <- intersect(patient.gene, patient.clinical)

# Only keep the patients that have both gene exp and BCR info
df.gene.dkfz.tpm <- df.gene.dkfz.tpm[, id.inter]
df.bcr.dkfz <- df.bcr.dkfz[id.inter, ]

save(df.gene.dkfz.tpm, df.bcr.dkfz, file = "../data/prad_data_dkfz.RData")

