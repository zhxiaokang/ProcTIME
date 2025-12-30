# adding power to grading system, such as Gleason Score

library(tibble)
library(dplyr)
library(survival)
library(survminer)
library(stringr)
library(ggplot2)

rm(list = ls())

# load color palette
source("color_palette.R")

plot_surv <- function(bcr.time, bcr.event, event.symbol, cluster, grade, corhort) {
  # Survival analysis 
  
  cluster <- factor(cluster, levels = c("TCE", "EPCE", "TASCE"))
  grade <- factor(grade, levels = c("low", "intermediate", "high"))
  
  df.surv <- data.frame(bcr.time, bcr.event, cluster, grade)
  
  df.surv$bcr.time <- as.numeric(bcr.time)
  df.surv$bcr.event <- 1*(bcr.event == event.symbol)
  
  km.fit <- survminer::surv_fit(Surv(bcr.time, bcr.event) ~ cluster, data = df.surv)
  p.surv.cluster <- survminer::ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE,
                                  title = paste("KM plot on", corhort),
                                  palette = color.three.clusters,
                                  xlab = "Time in months", ylab = "BCR free probability")
  
  km.fit <- survminer::surv_fit(Surv(bcr.time, bcr.event) ~ grade, data = df.surv)
  p.surv.gleason <- survminer::ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE,
                                          title = paste("KM plot on", corhort),
                                          palette = color.gleason.score,
                                          xlab = "Time in months", ylab = "BCR free probability")
  
  df.surv$merge <- paste(df.surv$cluster, df.surv$grade, sep = ", ")
  
  km.fit <- survminer::surv_fit(Surv(bcr.time, bcr.event) ~ merge, data = df.surv)
  p.surv.merge <- survminer::ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE,
                                          title = paste("KM plot on", corhort),
                                          xlab = "Time in months", ylab = "BCR free probability")
  
  return(list(p.surv.cluster, p.surv.gleason, p.surv.merge))
}

cox_uni_regression <- function(bcr.time, bcr.event, event.symbol, group) {
  # Survival analysis 
  # group: can be cluster or any grading stages
  
  df.surv <- data.frame(bcr.time, bcr.event, group)
  
  df.surv$bcr.time <- as.numeric(bcr.time)
  df.surv$bcr.event <- 1*(bcr.event == event.symbol)
  
  cox.reg <- coxph(Surv(bcr.time, bcr.event) ~ group, data = df.surv)
  
  return(cox.reg)
}

cox_multi_regression <- function(bcr.time, bcr.event, event.symbol, cluster, grade1, grade2) {
  # Survival analysis
  
  if(nargs() == 5){
    df.surv <- data.frame(bcr.time, bcr.event, cluster, grade1)
    
    df.surv$bcr.time <- as.numeric(bcr.time)
    df.surv$bcr.event <- 1*(bcr.event == event.symbol)
    
    cox.reg <- coxph(Surv(bcr.time, bcr.event) ~ cluster + grade1, data = df.surv)
    
    return(cox.reg)
  } else{
    df.surv <- data.frame(bcr.time, bcr.event, cluster, grade1, grade2)
    
    df.surv$bcr.time <- as.numeric(bcr.time)
    df.surv$bcr.event <- 1*(bcr.event == event.symbol)
    
    cox.reg <- coxph(Surv(bcr.time, bcr.event) ~ cluster + grade1 + grade2, data = df.surv)
    
    return(cox.reg)
  }
  
}

cox_interaction_regression <- function(bcr.time, bcr.event, event.symbol, cluster, grade) {
  # Survival analysis
  
  df.surv <- data.frame(bcr.time, bcr.event, cluster, grade)
  
  df.surv$bcr.time <- as.numeric(bcr.time)
  df.surv$bcr.event <- 1*(bcr.event == event.symbol)
  
  cox.reg <- coxph(Surv(bcr.time, bcr.event) ~ cluster * grade, data = df.surv)
  
  return(cox.reg)
}

# ====== TCGA ======
load('../data/clustering_tcga_removing_badguys.RData')

# pT stage
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

# gleason score
gleason.score <- df.pt23$stage_event_gleason_grading
index.intermediate.gs <- c(grep("^743", gleason.score))
index.high.gs <- c(grep("^8", gleason.score), grep("^9", gleason.score), grep("^10", gleason.score))

gleason.score.cat <- gleason.score
gleason.score.cat[index.high.gs] <- "high"
gleason.score.cat[index.intermediate.gs] <- "intermediate"
gleason.score.cat[-c(index.high.gs, index.intermediate.gs)] <- "low"

gleason.score.cat <- factor(gleason.score.cat, levels = c("low", "intermediate", "high"))

# do Cox regression
df.cox <- data.frame(pT_stage = pt23, GleasonScore = gleason.score.cat,
                     cluster = df.pt23$cluster,
                     bcr_months = df.pt23$bcr_months,
                     biochemical_recurrence = df.pt23$biochemical_recurrence)

df.cox$pT_stage <- factor(df.cox$pT_stage, levels = c("T2", "T3a", "T3b"))

df.cox$cluster <- factor(df.cox$cluster, levels = c("TCE", "EPCE", "TASCE"))

# plots <- plot_surv(df.clinical.bcr.cluster$bcr_months, df.clinical.bcr.cluster$biochemical_recurrence, "1",
#           df.clinical.bcr.cluster$cluster, gleason.score.cat, "TCGA")

cox.tcga.uni.cluster <- cox_uni_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", 
                                           df.cox$cluster)
summary(cox.tcga.uni.cluster)
cox.zph(cox.tcga.uni.cluster)

cox.tcga.uni.gs <- cox_uni_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", 
                                      df.cox$GleasonScore)
summary(cox.tcga.uni.gs)
cox.zph(cox.tcga.uni.gs)

cox.tcga.uni.pt <- cox_uni_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", 
                                      df.cox$pT_stage)
summary(cox.tcga.uni.pt)
cox.zph(cox.tcga.uni.pt)

cox.tcga.multi.three <- cox_multi_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", 
                                       df.cox$cluster, df.cox$GleasonScore, df.cox$pT_stage)
summary(cox.tcga.multi.three)
cox.zph(cox.tcga.multi.three)

cox.tcga.multi.three <- cox_multi_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", 
                                             df.cox$cluster, df.cox$GleasonScore, df.cox$pT_stage)
summary(cox.tcga.multi.three)
cox.zph(cox.tcga.multi.three)

cox.tcga.multi.gs.pt <- cox_multi_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", 
                                             df.cox$GleasonScore, df.cox$pT_stage)
summary(cox.tcga.multi.gs.pt)
cox.zph(cox.tcga.multi.gs.pt)

cox.tcga.multi.cluster.gs <- cox_multi_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", 
                                             df.cox$cluster, df.cox$GleasonScore)
summary(cox.tcga.multi.cluster.gs)
cox.zph(cox.tcga.multi.cluster.gs)

cox.tcga.multi.cluster.pt <- cox_multi_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", 
                                                  df.cox$cluster, df.cox$pT_stage)
summary(cox.tcga.multi.cluster.pt)
cox.zph(cox.tcga.multi.cluster.pt)

df.ci <- data.frame("Gleason Score"=0.677, "Gleason Score + pT Stage"=0.722,
                    "Gleason Score + TIME Subtype"=0.742, "Gleason Score + pT Stage + TIME Subtype"=0.76)


# draw plot for CI
rect.ci <- data.frame(
  ymin <- c(1, 3, 5, 7),
  ymax <- ymin + 0.8,
  xmin <- rep(0.65, 4),
  xmax <- c(0.76, 0.742, 0.722, 0.677)
)

pdf("../output/CI_tcga.pdf", width = 4, height = 3)
ggplot() +
  geom_rect(data = rect.ci, aes(xmin = xmin, xmax = xmax, 
                                ymin = ymin, ymax = ymax), fill = "black") +
  geom_vline(xintercept=0.65, colour = "grey50") +
  theme(axis.text.y = element_blank(), panel.background = element_blank(), 
        axis.line.x = element_line(colour = "grey50"),
        axis.ticks.y = element_blank())
dev.off()



# cox.tcga.interaction <- cox_interaction_regression(df.clinical.bcr.cluster$bcr_months, df.clinical.bcr.cluster$biochemical_recurrence, "1", 
#                                              df.clinical.bcr.cluster$cluster, gleason.score.cat)
# summary(cox.tcga.interaction)

# ====== Long/Atlanta ======

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

gleason.score.cat <- factor(gleason.score.cat, levels = c("low", "intermediate", "high"))

pt <- df.bcr.long$stage
# > table(pt)
# pt
# T1C  T2 T2A T2B T2C  T3 T3A T3B  T4 
# 14  10  23  10  30   2   6   9   1 
pt <- str_match(pt, "T[0-9]")
# > table(pt)
# > table(pt)
# pt
# T1 T2 T3 T4 
# 14 73 17  1 

df.cox <- data.frame(pT_stage = pt, GleasonScore = gleason.score.cat,
                     prediction = df.pred$pred,
                     bcr_months = df.bcr.long$BCR_Months,
                     biochemical_recurrence = df.bcr.long$BCR)

df.cox <- dplyr::filter(df.cox, pT_stage != "T4")
df.cox$pT_stage <- factor(df.cox$pT_stage, levels = c("T1", "T2", "T3"))

df.cox$prediction <- factor(df.cox$prediction, levels = c("TCE", "EPCE", "TASCE"))

# plots <- plot_surv(df.cox$bcr_months, df.cox$biochemical_recurrence, "1",
#                    df.cox$prediction, df.cox$GleasonScore, "Atlanta")

cox.uni.cluster <- cox_uni_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", df.cox$prediction)
summary(cox.uni.cluster)
cox.zph(cox.uni.cluster)

cox.uni.gs <- cox_uni_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", df.cox$GleasonScore)
summary(cox.uni.gs)
cox.zph(cox.uni.gs)

cox.uni.pt <- cox_uni_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", df.cox$pT_stage)
summary(cox.uni.pt)
cox.zph(cox.uni.pt)

# cox.long.multi <- cox_multi_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1",
#                                        df.cox$prediction, df.cox$GleasonScore, df.cox$pT_stage)
# summary(cox.long.multi)
# cox.zph(cox.long.multi)

cox.long.multi <- cox_multi_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1",
                                       df.cox$prediction, df.cox$GleasonScore)
summary(cox.long.multi)
cox.zph(cox.long.multi)

# ====== DKFZ/Hamburg ======

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

gleason.score <- df.pt23$Radical.Prostatectomy.Gleason.Score.for.Prostate.Cancer

# unique(gleason.score)
# [1] "3+4" "5+4" "4+3" "3+3" "4+5" "5+5" "4+4"

index.low.gs <- c(grep("^3", gleason.score))
index.intermediate.gs <- grep("4+3", gleason.score, fixed = TRUE)

gleason.score.cat <- gleason.score
gleason.score.cat[index.low.gs] <- "low"
gleason.score.cat[index.intermediate.gs] <- "intermediate"
gleason.score.cat[-c(index.low.gs, index.intermediate.gs)] <- "high"

# `%ni%` <- Negate(`%in%`)

GleasonScore <- factor(gleason.score.cat, levels = c("low", "intermediate", "high"))

df.cox <- data.frame(pT_stage = pt23, GleasonScore = GleasonScore,
                     prediction = df.pt23$pred,
                     bcr_months = df.pt23$TIME_FROM_SURGERY_TO_BCR_LASTFU,
                     biochemical_recurrence = df.pt23$BCR_STATUS)

df.cox$pT_stage <- factor(df.cox$pT_stage, levels = c("T2", "T3a", "T3b"))

df.cox$prediction <- factor(df.cox$prediction, levels = c("TCE", "EPCE", "TASCE"))

# plots <- plot_surv(df.cox$bcr_months, df.cox$biochemical_recurrence, "1",
#                    df.cox$prediction, df.cox$GleasonScore, "Atlanta")

cox.uni.cluster <- cox_uni_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", df.cox$prediction)
summary(cox.uni.cluster)
cox.zph(cox.uni.cluster)

cox.uni.pt <- cox_uni_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", df.cox$pT_stage)
summary(cox.uni.pt)
cox.zph(cox.uni.pt)

cox.uni.gs <- cox_uni_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1", df.cox$GleasonScore)
summary(cox.uni.gs)
cox.zph(cox.uni.gs)

# cox.long.multi <- cox_multi_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1",
#                                        df.cox$prediction, df.cox$GleasonScore, df.cox$pT_stage)
# summary(cox.long.multi)
# cox.zph(cox.long.multi)

cox.dkfz.multi <- cox_multi_regression(df.cox$bcr_months, df.cox$biochemical_recurrence, "1",
                                       df.cox$prediction, df.cox$pT_stage, df.cox$GleasonScore)
summary(cox.dkfz.multi)
cox.zph(cox.dkfz.multi)


