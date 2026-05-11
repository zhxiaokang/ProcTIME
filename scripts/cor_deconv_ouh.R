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
# # calculate correlations of all sample pairs

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

sim.diff.patient <- sim_pairs(res.deconv.ouh.scale, list.pair.diff.patient)

df.box.sim <- data.frame("sim" = c(sim.diff.patient, sim.diff.foci, sim.same.focus), 
                         "focus" = c(rep("Different patients", length(sim.diff.patient)), 
                                     rep("Different foci", length(sim.diff.foci)),
                                     rep("Same focus", length(sim.same.focus))))

all_comparisons <- list(c("Different patients", "Different foci"), c("Different foci", "Same focus"), c("Different patients", "Same focus"))

p <- ggboxplot(df.box.sim, x = "focus", y = "sim",
               color = "focus", palette = "jco",
               add = "jitter") + labs(x= "", y="Similarity")

p.3comparisons.sim <- p + 
  stat_compare_means(comparisons = all_comparisons, method = 'wilcox.test',
                     label = "p.signif") +
  theme(legend.position = "none")

pdf("../output/corr_3comparisons_similarity.pdf", width = 6, height = 4)
print(p.3comparisons.sim)
dev.off()

# Plot with exact P-values
p.3comparisons.sim.pval <- p + 
  stat_compare_means(comparisons = all_comparisons, method = 'wilcox.test',
                     label = "p.format") +
  theme(legend.position = "none")

pdf("../output/corr_3comparisons_similarity_pval.pdf", width = 6, height = 4)
print(p.3comparisons.sim.pval)
dev.off()
