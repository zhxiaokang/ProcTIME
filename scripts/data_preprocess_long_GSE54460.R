# preprocess GSE54460

library(openxlsx)
library(tibble)
library(dplyr)

rm(list=ls())

# ====== gene expression data preprocess ======
long.fpkm <- read.xlsx("../data/GSE54460/GSE54460_FPKM.xlsx", rowNames = TRUE)

# replace NA with 0
long.fpkm[is.na(long.fpkm)] <- 0

# convert some string cols into numeric
m.fpkm <- apply(long.fpkm, 2, as.numeric)
row.names(m.fpkm) <- row.names(long.fpkm)
long.fpkm <- as.data.frame(m.fpkm)

# function to convert FPKM to TPM
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

df.gene.long.tpm <- apply(long.fpkm, 2, fpkmToTpm)

# ====== clinical data preprocess ======
long.clinical <- read.xlsx("../data/GSE54460/GSE54460_clinical_info.xlsx", rowNames = TRUE)

df.clinical.long <- as.data.frame(t(long.clinical))

# function to merge "months to BCR" and "last Follow-up"
merge_bcr <- function(months.to.bcr, last.follow.up, bcr.status) {
  bcr.months <- months.to.bcr
  
  for (i in seq(1, length(bcr.status))){
    if (bcr.status[i] == "0"){
      bcr.months[i] <- last.follow.up[i]
    }
  }
  
  return(as.numeric(bcr.months))
}

bcr.months <- merge_bcr(df.clinical.long$`Months to BCR`, df.clinical.long$`Months total F/U`, df.clinical.long$BCR)

df.bcr.long <- df.clinical.long
df.bcr.long$BCR_Months <- bcr.months

sum(colnames(df.gene.long.tpm) != rownames(df.bcr.long))

# # remove patients of 0 BCR_Months
# df.bcr.long <- filter(df.bcr.long, BCR_Months != 0)

# then intersect the patients
patients.bcr <- rownames(df.bcr.long)
patients.gene <- colnames(df.gene.long.tpm)

patients.inter <- intersect(patients.bcr, patients.gene)

df.bcr.long <- df.bcr.long[patients.inter, ]
df.gene.long.tpm <- df.gene.long.tpm[, patients.inter]

save(df.gene.long.tpm, df.bcr.long, file = "../data/prad_data_long.RData")

