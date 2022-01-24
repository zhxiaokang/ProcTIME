# To download and pre-process TCGA data

# load required packages
libs <- c("TCGAbiolinks", "RTCGAToolbox", "tibble", "dplyr")

invisible(lapply(libs, library, character.only = TRUE))

# clinical data
# Date: 20.1.2022
query <- GDCquery(project = "TCGA-PRAD",
                  data.category = "Clinical",
                  file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query,clinical.info = "patient")

# =========== pre-process the clinical data ===========
# remove duplicate rows
df.clinical <- distinct(clinical)
df.clinical <- column_to_rownames(df.clinical, colnames(df.clinical)[1])
write.table(df.clinical, "../data/prad_clinical_tcga.txt", sep = "\t", col.names = NA, row.names = TRUE)

# pick out the patients with biochemical_recurrence info
df.clinical.with.bcr.info <- filter(df.clinical, df.clinical$biochemical_recurrence %in% c("NO", "YES"))
print(paste("Num of patients with BCR info:", nrow(df.clinical.with.bcr.info)))

# remove the patients with biochemical_recurrence but no days to first biochemical recurrence
df.clinical.bcr <- filter(df.clinical.with.bcr.info, !(biochemical_recurrence == "YES" & is.na(days_to_first_biochemical_recurrence)))
print(paste("Num of patients with correct BCR 'YES' info:", nrow(df.clinical.bcr)))

# remove the patients that "biochemical_recurrence" == "NO" but has a "days_to_first_biochemical_recurrence"
df.clinical.bcr <- filter(df.clinical.bcr, !(biochemical_recurrence == "NO" & !is.na(days_to_first_biochemical_recurrence)))
print(paste("Num of patients with correct BCR 'NO' info:", nrow(df.clinical.bcr)))

# if biochemical_recurrence == "NO" then use days_to_last_followup
# merge the 2 cols (days_to_first_biochemical_recurrence, days_to_last_followup) into one col: bcr_days
## - use days_to_first_biochemical_recurrence if biochemical_recurrence == YES, 
## - use days_to_last_followup if biochemical_recurrence == NO
bcr.days <- coalesce(df.clinical.bcr$days_to_last_followup, df.clinical.bcr$days_to_death)
for (i in seq(1, nrow(df.clinical.bcr))) {
  if (df.clinical.bcr$biochemical_recurrence[i] == "YES") {
    bcr.days[i] <- df.clinical.bcr$days_to_first_biochemical_recurrence[i]
  }
}
df.clinical.bcr$bcr_days <- bcr.days

# modify the patient ID
patient.id <- row.names(df.clinical.bcr)
row.names(df.clinical.bcr) <- gsub("-", ".", patient.id)

# =========== gene expression data ============
# Date: 19.1.2022
pradData <- getFirehoseData(dataset="PRAD", runDate="20160128",
                            forceDownload=TRUE, clinical=FALSE, RNASeq2GeneNorm = TRUE)

prad.gene.intermediate <- getData(pradData, type = "RNASeq2GeneNorm")
prad.gene <- getElement(prad.gene.intermediate[[1]], "DataMatrix")

# ========== pre-process the gene expression data =========
# modify the patient ID
rid <- substr(colnames(prad.gene),1,12)
rid <- gsub("-", ".", rid)
colnames(prad.gene) <- rid

# ====== intersect patients of gene and clinical info =======
id.inter <- intersect(colnames(prad.gene), rownames(df.clinical.bcr))

df.clinical.bcr <- df.clinical.bcr[id.inter, ]
prad.gene <- prad.gene[, id.inter]

# write to file
write.table(df.clinical.bcr, "../data/prad_clinical_bcr_tcga.txt", sep = "\t", col.names = NA, row.names = TRUE)
write.table(prad.gene, "../data/prad_gene_tcga.txt", sep = "\t", col.names = NA, row.names = TRUE)

save(df.clinical.bcr, prad.gene, file = "../data/prad_data_tcga.RData")

