# Draw PCA of TIME of TCGA

rm(list = ls())

library(dplyr)
library(factoextra)
library(ggfortify)
# library(pca3d)
library(plotly)

# load the data
load("../data/tcga_time_deconv.RData")

# pick out the df
df.time <- dplyr::select(time.epic.bcr, -c(bcr_days, biochemical_recurrence))

# draw PCA
pca <- prcomp(df.time, scale. = TRUE)

p.pca <- fviz_pca_ind(pca) + xlim(-5, 15) + ylim (-25, 10)
print(p.pca)

# outlier: TCGA.4L.AA1F

df.time <- df.time[!(rownames(df.time) %in% "TCGA.4L.AA1F"), ]

# draw PCA
pca <- prcomp(df.time, scale. = TRUE)

p.pca <- fviz_pca_ind(pca, geom = c("point", "text")) + xlim(-5, 6) + ylim (-4, 5)
print(p.pca)

# ====== 3D PCA with lib 'pca3d' ======
# pca3d(pca, show.labels = TRUE)

# outliers: TCGA.J4.A6G1, TCGA.J9.A8CK, TCGA.YL.A9WY

# patient.outlier <- c("TCGA.4L.AA1F", "TCGA.J4.A6G1", "TCGA.J9.A8CK", "TCGA.YL.A9WY")

# ====== 3D PCA with lib 'plotly' ======
components <- as.data.frame(pca[['x']])
plot_ly(components, x = ~PC1, y = ~PC2, z = ~PC3, text = rownames(df.time))

# outliers: TCGA.J4.A6G1, TCGA.J9.A8CK

patient.outlier <- c("TCGA.4L.AA1F", "TCGA.J4.A6G1", "TCGA.J9.A8CK")  # "TCGA.4L.AA1F" detected from 2D above

# ==== remove the 3 patients for downstream analysis (clustering and classification) ======

df.clinical.bcr <- df.clinical.bcr[!(rownames(df.clinical.bcr) %in% patient.outlier), ]

prad.gene <- prad.gene[, !(colnames(prad.gene) %in% patient.outlier)]

res.deconv.epic <- res.deconv.epic[, !(colnames(res.deconv.epic) %in% patient.outlier)]

time.epic.bcr <- time.epic.bcr[!(rownames(time.epic.bcr) %in% patient.outlier), ]

save(patient.outlier, df.clinical.bcr, prad.gene, res.deconv.epic, time.epic.bcr, file = "../data/tcga_time_deconv_rm_outlier.RData")
