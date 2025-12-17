# deconvolute the TIME for TCGA dataset

library(tibble)
library(immunedeconv)

load("../data/prad_data_tcga.RData")

res.deconv <- deconvolute_estimate(as.matrix(prad.gene))

save(res.deconv, file = "../data/ESTIMATE_TCGA.RData")