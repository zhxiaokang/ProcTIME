# association between "uncharacterized cell" and "tumor purity"

library(openxlsx)
library(ggpubr)

# load TIME OUH data
load("../data/ouh_time_deconv.RData")

# load cancer purity data
tumor.purity <- read.xlsx('../data/OUH/SupplementaryTables_NOT_FINAL_VERSION.xlsx', sheet = 7, startRow = 4)

df.merge <- merge(df.bcr.ouh, df.sample.deconv, by = "row.names")

df.merge.purity <- merge(df.merge, tumor.purity, by.x = "sample_name_lovf_et_al", by.y = "Sample")

pdf("../output/association_epithelial_tumor_purity_oslo.pdf", width = 6, height = 4)
p <- ggscatter(df.merge.purity, x = "uncharacterized.cell", y = "Tumor_purity", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Proportions of epithelial cell", ylab = "Tumor purity")
print(p)
dev.off()

# their distributions
gghistogram(df.merge.purity, x = "uncharacterized.cell")

gghistogram(df.merge.purity, x = "Tumor_purity")