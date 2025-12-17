# Compare the correlation/similarity of cell types' composition of intra- and inter- foci and patients

# 3 lists:
# - samples from the same focus
# - samples from different foci but from the same patient
# - samples from different patients (randomly selecting one sample representing that patient. Do sampling several times)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(stats)
library(lsa)

rm(list = ls())

load("../data/ouh_time_deconv.RData")

sample.info <- df.bcr.ouh

# scale the cell types' proportions
## z-score transformation --- 0 mean 1 deviation
res.deconv.ouh.scale <- data.frame(t(scale(t(res.deconv.ouh))))
df.sample.deconv.scale <- data.frame(t(res.deconv.ouh.scale))

# ====== function to calculate the corr of all input pairs ======
# calculate correlations of all sample pairs
cor_pairs <- function(df, list.pair) {
  # df: df of all samples and their cell fraction
  # list.pair: list of all pairs of samples (in each pair, number of samples is 2+)
  list.cor <- c()  # vector to collect all correlations
  for (pair in list.pair) {
    df.pair <- dplyr::select(df, dplyr::all_of(pair))
    m.cor <- cor(df.pair, method = "spearman")
    list.cor <- c(list.cor, m.cor[upper.tri(m.cor)])
  }
  return(list.cor)
}

dist_pairs <- function(df, list.pair) {
  list.dist <- c()
  
  for (pair in list.pair){
    df.pair <- df[all_of(pair), ]
    m.dist <- dist(df.pair, method = "euclidean")
    list.dist <- c(list.dist, c(m.dist))
  }
  return(list.dist)
}

sim_pairs <- function(df, list.pair) {
  # df: df of all samples and their cell fraction
  # list.pair: list of all pairs of samples (in each pair, number of samples is 2+)
  list.sim <- c()  # vector to collect all correlations
  for (pair in list.pair) {
    df.pair <- dplyr::select(df, all_of(pair))
    m.sim <- cosine(as.matrix(df.pair))
    # m.sim <- cor(df.pair, method = "pearson")
    list.sim <- c(list.sim, m.sim[upper.tri(m.sim)])
  }
  return(list.sim)
}

# ====== Pick out the qualified patients ======
# multifocal patients with multi-sample focus, 
#     or to say: the patients that have samples from multi foci, and at least one of the foci has multi samples
# so that, the patients being included in the 3 (especially the first 2) lists are the same patients
patient.number <- unique(sample.info$patient_number)
patient.keep <- c()
for (i in patient.number) {
  # the df for one patient
  df.patient <- dplyr::filter(sample.info, patient_number == i)
  # exclude patients that only have samples from 1 focus
  if (length(unique(df.patient$localized.in.focus.number)) < 2) {
    print("only one focus")
    next
  }
  # the focus that has more than 1 sample
  focus.rep <- names(which(table(dplyr::select(df.patient, localized.in.focus.number))>1))
  if (length(focus.rep) == 0) {
    print("single sample for all foci")
    next
  }
  patient.keep <- c(patient.keep, i)
}

# ====== construct the pairs of samples from the same focus ======
# build a list of all sample pairs from the same focus
list.pair.same.focus <- list()

for (i in patient.keep) {
  # the df for one patient
  df.patient <- dplyr::filter(sample.info, patient_number == i)
  
  # the focus that has more than 1 sample
  focus.rep <- names(which(table(dplyr::select(df.patient, localized.in.focus.number))>1))
  
  # the sample id that belongs to that focus
  sample.id <- rownames(dplyr::filter(df.patient, localized.in.focus.number == focus.rep))
  sample.pair <- list(sample.id)
  list.pair.same.focus <- c(list.pair.same.focus, sample.pair)
}

cor.same.focus <- cor_pairs(res.deconv.ouh.scale, list.pair.same.focus)
dist.same.focus <- dist_pairs(df.sample.deconv.scale, list.pair.same.focus)
sim.same.focus <- sim_pairs(res.deconv.ouh.scale, list.pair.same.focus)

# ====== construct the pairs of samples from different foci but the same patient ======
list.pair.diff.foci <- list()

for (i in patient.keep) {
  # the df for one patient
  df.patient <- dplyr::filter(sample.info, patient_number == i)
  
  # the focus that has more than 1 sample
  focus.rep <- names(which(table(dplyr::select(df.patient, localized.in.focus.number))>1))
  
  sample.id.rep.focus <- rownames(dplyr::filter(df.patient, localized.in.focus.number == focus.rep))  # there can be either 2 or 3 samples
  sample.id.another.focus <- rownames(dplyr::filter(df.patient, localized.in.focus.number != focus.rep))  # there is always only 1 sample
  
  for (sample.rep in sample.id.rep.focus){
    sample.pair <- list(c(sample.rep, sample.id.another.focus))
    list.pair.diff.foci <- c(list.pair.diff.foci, sample.pair)
  }
}

cor.diff.foci <- cor_pairs(res.deconv.ouh.scale, list.pair.diff.foci)
dist.diff.foci <- dist_pairs(df.sample.deconv.scale, list.pair.diff.foci)
sim.diff.foci <- sim_pairs(res.deconv.ouh.scale, list.pair.diff.foci)

# ====== construct the pairs of samples from different patients ======

vec.pair.diff.patient <- c()

for (i in patient.keep) {
  df.patient <- dplyr::filter(sample.info, patient_number == i)
  
  sample.ids <- rownames(df.patient)
  
  sample.id <- sort(sample.ids)[1]
  
  vec.pair.diff.patient <- c(vec.pair.diff.patient, sample.id)
}

list.pair.diff.patient <- list(vec.pair.diff.patient)

cor.diff.patient <- cor_pairs(res.deconv.ouh.scale, list.pair.diff.patient)
dist.diff.patient <- dist_pairs(df.sample.deconv.scale, list.pair.diff.patient)
sim.diff.patient <- sim_pairs(res.deconv.ouh.scale, list.pair.diff.patient)

df.box.cor <- data.frame("cor" = c(cor.diff.patient, cor.diff.foci, cor.same.focus), 
                         "focus" = c(rep("Different patients", length(cor.diff.patient)), 
                                     rep("Different foci", length(cor.diff.foci)),
                                     rep("Same focus", length(cor.same.focus))))

df.box.dist <- data.frame("dist" = c(dist.diff.patient, dist.diff.foci, dist.same.focus), 
                          "focus" = c(rep("Different patients", length(dist.diff.patient)), 
                                      rep("Different foci", length(dist.diff.foci)),
                                      rep("Same focus", length(dist.same.focus))))

df.box.sim <- data.frame("sim" = c(sim.diff.patient, sim.diff.foci, sim.same.focus), 
                         "focus" = c(rep("Different patients", length(sim.diff.patient)), 
                                     rep("Different foci", length(sim.diff.foci)),
                                     rep("Same focus", length(sim.same.focus))))

all_comparisons <- list(c("Different patients", "Different foci"), c("Different foci", "Same focus"), c("Different patients", "Same focus"))

p <- ggboxplot(df.box.cor, x = "focus", y = "cor",
               color = "focus", palette = "jco",
               add = "jitter", short.panel.labs = TRUE)

p.3comparisons.cor <- p + stat_compare_means(comparisons = all_comparisons, label = "p.format")

p <- ggboxplot(df.box.dist, x = "focus", y = "dist",
               color = "focus", palette = "jco",
               add = "jitter", short.panel.labs = TRUE)

p.3comparisons.dist <- p + stat_compare_means(comparisons = all_comparisons, label = "p.format")

p <- ggboxplot(df.box.sim, x = "focus", y = "sim",
               color = "focus", palette = "jco",
               add = "jitter") + labs(x= "", y="Similarity")

# p.3comparisons.sim <- p + 
#                       stat_compare_means(comparisons = all_comparisons, method = 'wilcox.test',
#                                          aes(label = paste("P =", formatC(..p.format.., format = "f", digits = 2)))) +
#                       theme(legend.position = "none")

p.3comparisons.sim <- p + 
  stat_compare_means(comparisons = all_comparisons, method = 'wilcox.test',
                     label = "p.signif") +
  theme(legend.position = "none")

pdf("../output/corr_3comparisons_similarity.pdf", width = 6, height = 4)
print(p.3comparisons.sim)
dev.off()

# Load necessary libraries
# library(ggplot2)
# library(ggpubr)

# # Create a sample data frame
# data <- data.frame(
#   Group = rep(c("A", "B", "C"), each = 20),
#   Value = c(rnorm(20, mean = 0), rnorm(20, mean = 1), rnorm(20, mean = 2))
# )

# # Create a boxplot with custom p-value label
# p <- ggboxplot(data, x = "Group", y = "Value")

# # Customize p-value label
# p + stat_compare_means(comparisons = list(c("A", "B"), c("A", "C"), c("B", "C")), method = "t.test",
#                        aes(label = paste("P = ", formatC(..p.format.., format = "f", digits = 2))))

# p
# 
# df.box.cor <- data.frame("cor" = c(cor.same.focus, cor.diff.foci), "focus" = c(rep("Same focus", length(cor.same.focus)), rep("Different foci", length(cor.diff.foci))))
# 
# p <- ggboxplot(df.box.cor, x = "focus", y = "cor",
#                color = "focus", palette = "jco",
#                add = "jitter", short.panel.labs = TRUE,
#                xlab = "Focus", ylab = "Spearman Correlation")
# p.2comparisons <- p + stat_compare_means()
# 
# pdf("../output/corr_2comparisons.pdf", width = 6, height = 4)
# print(p.2comparisons)
# dev.off()


