# pre-process MSKCC data

library(readxl)
library(tibble)
library(dplyr)

mskcc.gene <- read.csv("../data/prad_mskcc/MSKCC_PCa_mRNA_data.txt", sep = "\t", row.names = 1, na.strings = "")

mskcc.clinical <- read_excel("../data/prad_mskcc/MSKCC_PCa_Clinical_Annotation.xls")

mskcc.clinical <- column_to_rownames(mskcc.clinical, colnames(mskcc.clinical)[1])

# only pick out the intersecting patients
id.inter <- intersect(colnames(mskcc.gene), rownames(mskcc.clinical))

mskcc.gene.inter <- mskcc.gene[, id.inter]

mskcc.clinical.inter <- mskcc.clinical[id.inter, ]

# =========== preprocess the clinical data ===========
# remove the patients without BCR info
df.bcr.mskcc <- filter(mskcc.clinical.inter, !(BCR_Event == "NA"))

# filter those people out from the gene expression matrix
df.gene.mskcc <- mskcc.gene.inter[, rownames(df.bcr.mskcc)]

save(df.bcr.mskcc, df.gene.mskcc, file = "../data/prad_data_mskcc.RData")