# associations between 3 subtypes and some variabels, such as gleason score

library(dplyr)
library(ggstatsplot)
library(stringr)

# ====== on TCGA dataset ======
rm(list = ls())

load('../data/clustering_tcga_original.RData')

gleason.score <- df.clinical.bcr.cluster$stage_event_gleason_grading
index.intermediate.gs <- c(grep("^743", gleason.score))
index.high.gs <- c(grep("^8", gleason.score), grep("^9", gleason.score), grep("^10", gleason.score))

gleason.score.cat <- gleason.score
gleason.score.cat[index.high.gs] <- "high"
gleason.score.cat[index.intermediate.gs] <- "intermediate"
gleason.score.cat[-c(index.high.gs, index.intermediate.gs)] <- "low"

# `%ni%` <- Negate(`%in%`)

df.assoc <- data.frame(GleasonScore = factor(gleason.score.cat, levels = c("low", "intermediate", "high")),
                       cluster = df.clinical.bcr.cluster$cluster)

dat <- table(df.assoc$GleasonScore, df.assoc$cluster)

chisq <- chisq.test(dat)

p.chisq <- ggbarstats(
  df.assoc, GleasonScore, cluster,
  xlab = "Clustered cell type composition subtype",
  results.subtitle = FALSE,
  package = "RColorBrewer",
  palette = "Set2",
  subtitle = paste0(
    "TCGA, Pearson's Chi-squared test", ", p-value = ",
    ifelse(chisq$p.value < 0.001, "< 0.001", round(chisq$p.value, 3))
  )
)

pdf("../output/association_subtypes_gleasonscore_tcga.pdf", width = 6, height = 4)
print(p.chisq)
dev.off()

pt <- df.clinical.bcr.cluster$ajcc_pathologic_t

# table(pt)
# pt
# pt
# '-- T2a T2b T2c T3a T3b  T4 
#   5   7  10 132 131 110   7 

# pick out only T2 and T3
df.pt23 <- df.clinical.bcr.cluster %>% dplyr::filter(str_detect(ajcc_pathologic_t, '^T2|^T3'))

# merge T2* into T2
pt23 <- df.pt23$ajcc_pathologic_t
pt23 <- str_replace(pt23, 'T2[abc]', "T2")
# table(pt23)
# pt23
# T2 T3a T3b 
# 149 131 110 

df.assoc <- data.frame(pT_stage = factor(pt23),
                       cluster = df.pt23$cluster)
df.assoc <- na.omit(df.assoc)
df.assoc$pT_stage <- factor(df.assoc$pT_stage, levels = c("T2", "T3a", "T3b"))

dat <- table(df.assoc$pT_stage, df.assoc$cluster)
chisq <- chisq.test(dat)

p.chisq <- ggbarstats(
  df.assoc, pT_stage, cluster,
  xlab = "Clustered cell type composition subtype",
  results.subtitle = FALSE,
  package = "RColorBrewer",
  palette = "Set3",
  subtitle = paste0(
    "TCGA, Pearson's Chi-squared test", ", p-value = ",
    ifelse(chisq$p.value < 0.001, "< 0.001", round(chisq$p.value, 3))
  )
)

pdf("../output/association_subtypes_ptstage_tcga.pdf", width = 6, height = 4)
print(p.chisq)
dev.off()

# ==== load OUH datasets ====
rm(list = ls())

load("../data/prediction_3_ouh.RData")
cluster <- df.pred$pred

load("../data/ouh_time_deconv.RData")

gs <- df.bcr.ouh$gleason_score_updated
gs.cat <- gs

# `%ni%` <- Negate(`%in%`)

gs.cat[gs.cat %in% c("0+3", "3+0", "3+3", "0+3+4", "3+4", "3+4+0", "3+4+5")] <- "low"
gs.cat[gs.cat %in% c("4+3")] <- "intermediate"
gs.cat[gs.cat %in% c("0+4", "0+4+4", "4+4", "4+5")] <- "high"

df.pred$GleasonScore <- factor(gs.cat, levels = c("low", "intermediate", "high"))

dat <- table(df.pred$GleasonScore, df.pred$pred)

fisher <- fisher.test(dat)

p.fisher <- ggbarstats(
  df.pred, GleasonScore, pred,
  xlab = "Predicted cell type composition subtype",
  results.subtitle = FALSE,
  package = "RColorBrewer",
  palette = "Set2",
  subtitle = paste0(
    "Oslo, Fisher's exact test", ", p-value = ",
    ifelse(fisher$p.value < 0.001, "< 0.001", round(fisher$p.value, 3))
  )
)

pdf("../output/association_subtypes_gleasonscore_ouh.pdf", width = 6, height = 4)
print(p.fisher)
dev.off()

# 
# ptnm <- df.bcr.ouh$clinical_outcome
# 
# # pt <- str_match(ptnm, "T[0-9][a-z]?")
# # table(pt)
# # pt
# # T1a T1b T1c  T2 T2a T2b T2c T3a T3b  T4 
# # 1   2 147   9  49  42  59  50  39   2 
# 
# pt <- str_match(ptnm, "T[0-9]")
# table(pt)
# # pt
# # T1  T2  T3  T4 
# # 150 159  89   2 
# 
# df.assoc <- data.frame(pT_stage = factor(pt),
#                        cluster = df.clinical.bcr.cluster$cluster)
# df.assoc <- na.omit(df.assoc)
# df.assoc <- dplyr::filter(df.assoc, pT_stage != "T4")
# df.assoc$pT_stage <- factor(df.assoc$pT_stage, levels = c("T1", "T2", "T3"))
# 
# dat <- table(df.assoc$pT_stage, df.assoc$cluster)
# dat <- dat[1:3, ]
# chisq <- chisq.test(dat)
# 
# p.chisq <- ggbarstats(
#   df.assoc, pT_stage, cluster,
#   results.subtitle = FALSE,
#   package = "nbapalettes",
#   palette = "bucks_00s",
#   subtitle = paste0(
#     "TCGA, Pearson's Chi-squared test", ", p-value = ",
#     ifelse(chisq$p.value < 0.001, "< 0.001", round(chisq$p.value, 3))
#   )
# )
# 
# pdf("../output/association_subtypes_ptstage_tcga.pdf", width = 6, height = 4)
# print(p.chisq)
# dev.off()

# ==== load Long datasets ====
rm(list = ls())

load("../data/prediction_3_long.RData")

load("../data/long_time_deconv.RData")

gleason.score <- df.bcr.long$`Gleason Score`

# unique(gleason.score)
# [1] "549" "347" "437" "336" "459" "448" "538" "325"

index.intermediate.gs <- c(grep("437$", gleason.score))
index.high.gs <- c(grep("8$", gleason.score), grep("9$", gleason.score))

gleason.score.cat <- gleason.score
gleason.score.cat[index.high.gs] <- "high"
gleason.score.cat[index.intermediate.gs] <- "intermediate"
gleason.score.cat[-c(index.high.gs, index.intermediate.gs)] <- "low"

# `%ni%` <- Negate(`%in%`)

df.pred$GleasonScore <- factor(gleason.score.cat, levels = c("low", "intermediate", "high"))

dat <- table(df.pred$GleasonScore, df.pred$pred)

fisher <- fisher.test(dat)

p.fisher <- ggbarstats(
  df.pred, GleasonScore, pred,
  xlab = "Predicted cell type composition subtype",
  results.subtitle = FALSE,
  package = "RColorBrewer",
  palette = "Set2",
  subtitle = paste0(
    "Atlanta, Fisher's exact test", ", p-value = ",
    ifelse(fisher$p.value < 0.001, "< 0.001", round(fisher$p.value, 4))
  )
)

pdf("../output/association_subtypes_gleasonscore_long.pdf", width = 6, height = 4)
print(p.fisher)
dev.off()


pt <- df.bcr.long$stage
table(pt)
# > table(pt)
# pt
# T1C  T2 T2A T2B T2C  T3 T3A T3B  T4 
# 14  10  23  10  30   2   6   9   1 
pt <- str_replace(pt, "T[1-2][A-C]?", "T1c/2")  # merge T1C, T2
pt <- str_replace(pt, "T[3][A]?$", "T3/3a")  # merge T3, T3A
pt <- str_replace(pt, "T3B", "T3b") 
table(pt)
# > table(pt)
# pt
# T1 T2 T3 T4 
# 14 73 17  1 

df.assoc <- data.frame(pT_stage = pt, prediction = df.pred$pred)

df.assoc <- dplyr::filter(df.assoc, pT_stage != "T4")
df.assoc$pT_stage <- factor(df.assoc$pT_stage, levels = c("T1c/2", "T3/3a", "T3b"))

df.assoc$prediction <- factor(df.assoc$prediction, levels = c("TCE", "EPCE", "TASCE"))

dat <- table(df.assoc$pT_stage, df.assoc$prediction)

fisher <- fisher.test(dat)

p.fisher <- ggbarstats(
  df.assoc, pT_stage, prediction,
  xlab = "Predicted cell type composition subtype",
  results.subtitle = FALSE,
  package = "RColorBrewer",
  palette = "Set3",
  subtitle = paste0(
    "Atlanta, Fisher's exact test", ", p-value = ",
    ifelse(fisher$p.value < 0.001, "< 0.001", round(fisher$p.value, 4))
  )
)

pdf("../output/association_subtypes_ptstage_long.pdf", width = 6, height = 4)
print(p.fisher)
dev.off()

# ==== load DKFZ datasets ====
rm(list = ls())

load("../data/prediction_3_dkfz.RData")

load("../data/dkfz_time_deconv.RData")

df.bcr.dkfz.pred <- merge(df.bcr.dkfz, df.pred, by = "row.names")

pt <- df.bcr.dkfz.pred$STAGE
# > table(pt)
# pt
# pT2a pT2c pT3a pT3b  pT4
#  8    48   10   13    3

# pick out only T2 and T3
df.pt23 <- df.bcr.dkfz.pred %>% dplyr::filter(str_detect(STAGE, '^pT2|^pT3'))

# merge T2* into T2
pt23 <- df.pt23$STAGE
pt23 <- str_replace(pt23, 'T2[abc]', "T2")
pt23 <- str_replace(pt23, '^p', "")
# table(pt23)
# pt23
# T2 T3a T3b 
# 56  10  13 

df.assoc <- data.frame(pT_stage = pt23,
                     prediction = df.pt23$pred)

df.assoc$pT_stage <- factor(df.assoc$pT_stage, levels = c("T2", "T3a", "T3b"))

df.assoc$prediction <- factor(df.assoc$prediction, levels = c("TCE", "EPCE", "TASCE"))

dat <- table(df.assoc$pT_stage, df.assoc$prediction)

fisher <- fisher.test(dat)

p.fisher <- ggbarstats(
  df.assoc, pT_stage, prediction,
  results.subtitle = FALSE,
  xlab = "Predicted cell type composition subtype",
  package = "RColorBrewer",
  palette = "Set3",
  subtitle = paste0(
    "Hamburg, Fisher's exact test", ", p-value = ",
    ifelse(fisher$p.value < 0.001, "< 0.001", round(fisher$p.value, 4))
  )
)

pdf("../output/association_subtypes_ptstage_dkfz.pdf", width = 6, height = 4)
print(p.fisher)
dev.off()

gleason.score <- df.bcr.dkfz$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer

# unique(gleason.score)
# [1] "3+4" "5+4" "4+3" "3+3" "4+5" "5+5" "4+4"

index.low.gs <- c(grep("^3", gleason.score))
index.intermediate.gs <- grep("4+3", gleason.score, fixed = TRUE)

gleason.score.cat <- gleason.score
gleason.score.cat[index.low.gs] <- "low"
gleason.score.cat[index.intermediate.gs] <- "intermediate"
gleason.score.cat[-c(index.low.gs, index.intermediate.gs)] <- "high"

# `%ni%` <- Negate(`%in%`)

df.pred$GleasonScore <- factor(gleason.score.cat, levels = c("low", "intermediate", "high"))

dat <- table(df.pred$GleasonScore, df.pred$pred)

fisher <- fisher.test(dat)

p.fisher <- ggbarstats(
  df.pred, GleasonScore, pred,
  xlab = "Predicted cell type composition subtype",
  results.subtitle = FALSE,
  package = "RColorBrewer",
  palette = "Set2",
  subtitle = paste0(
    "Hamburg, Fisher's exact test", ", p-value = ",
    ifelse(fisher$p.value < 0.001, "< 0.001", round(fisher$p.value, 4))
  )
)

pdf("../output/association_subtypes_gleasonscore_dkfz.pdf", width = 6, height = 4)
print(p.fisher)
dev.off()







# tce.low <- nrow(dplyr::filter(df.pred, pred == "TCE" & gleason_score == "low"))
# tce.intermediate <- nrow(dplyr::filter(df.pred, pred == "TCE" & gleason_score == "intermediate"))
# tce.high <- nrow(dplyr::filter(df.pred, pred == "TCE" & gleason_score == "high"))
# 
# tuce.low <- nrow(dplyr::filter(df.pred, pred == "TuCE" & gleason_score == "low"))
# tuce.intermediate <- nrow(dplyr::filter(df.pred, pred == "TuCE" & gleason_score == "intermediate"))
# tuce.high <- nrow(dplyr::filter(df.pred, pred == "TuCE" & gleason_score == "high"))
# 
# tasce.low <- nrow(dplyr::filter(df.pred, pred == "TASCE" & gleason_score == "low"))
# tasce.intermediate <- nrow(dplyr::filter(df.pred, pred == "TASCE" & gleason_score == "intermediate"))
# tasce.high <- nrow(dplyr::filter(df.pred, pred == "TASCE" & gleason_score == "high"))

