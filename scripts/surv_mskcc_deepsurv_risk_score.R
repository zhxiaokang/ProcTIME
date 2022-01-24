# The KM plot of MSKCC based on the predicted risk score from DeepSurv

library(dplyr)
library(tibble)

rm(list=ls())

load("../data/prad_data_mskcc.RData")

rs <- read.csv("../data/risk_score_from_DeepSurv_mskcc.csv")

rs <- column_to_rownames(rs, names(rs)[1])

df.bcr <- select(df.bcr.mskcc, c(BCR_Event, BCR_FreeTime))

df.surv <- merge(rs, df.bcr, by = "row.names")
df.surv <- column_to_rownames(df.surv, names(df.surv)[1])

# draw survival curv
# change data type
df.surv$BCR_FreeTime <- as.numeric(df.surv$BCR_FreeTime)
df.surv$BCR_Event <- 1*(df.surv$BCR_Event=="BCR_Algorithm")

# find the optimal cutpoint
res.cut <- surv_cutpoint(data = df.surv, time = "BCR_FreeTime", event = "BCR_Event", variables = "risk.score")

res.cat <- surv_categorize(res.cut)

km.fit <- survfit(Surv(BCR_FreeTime, BCR_Event) ~ risk.score, data = res.cat)
p.surv <- ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE, xlab = "Day")

pdf("../output/km_mskcc_on_risk_score.pdf")
print(p.surv)
dev.off()

