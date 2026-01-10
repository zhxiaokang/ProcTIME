# compare the deconvoluted results from EPIC and quanTIseq
library(tibble)
library(immunedeconv)
library(ggplot2)
library(gridExtra)

# load TCGA data
load("../data/prad_data_tcga.RData")

# deconvolution with EPIC
# res.deconv <- deconvolute(as.matrix(prad.gene), "epic", tumor = TRUE)
# res.deconv.epic <- column_to_rownames(res.deconv, colnames(res.deconv)[1])
# res.deconv.epic.sample <- as.data.frame(t(res.deconv.epic))
load("../data/tcga_time_deconv.RData")

# deconvolution with quanTIseq
# res.deconv <- deconvolute(as.matrix(prad.gene), "quantiseq", tumor = TRUE)
# res.deconv.quantiseq <- column_to_rownames(res.deconv, colnames(res.deconv)[1])
# res.deconv.quantiseq.sample <- as.data.frame(t(res.deconv.quantiseq))
load("../data/tcga_time_deconv_quantiseq.RData")

# res.deconv.quantiseq has these two rows: T cell CD4+ (non-regulatory), T cell regulatory (Tregs)
# combine them into one row: T cell CD4+
res.deconv.quantiseq["T cell CD4+", ] <- res.deconv.quantiseq["T cell CD4+ (non-regulatory)", ] + res.deconv.quantiseq["T cell regulatory (Tregs)", ]
# remove the "T cell CD4+ (non-regulatory)" and "T cell regulatory (Tregs)" rows
res.deconv.quantiseq <- res.deconv.quantiseq[!(rownames(res.deconv.quantiseq) %in% c("T cell CD4+ (non-regulatory)", "T cell regulatory (Tregs)")), ]

# Ensure patient ID columns are in the same order between EPIC and quanTIseq results
common_patients <- intersect(colnames(res.deconv.epic), colnames(res.deconv.quantiseq))
res.deconv.epic.ordered <- res.deconv.epic[, common_patients]
res.deconv.quantiseq.ordered <- res.deconv.quantiseq[, common_patients]

# remove the following patients identified as outliers in previous PCA analysis
patient.outlier <- c("TCGA.4L.AA1F", "TCGA.J4.A6G1", "TCGA.J9.A8CK")
res.deconv.epic.ordered <- res.deconv.epic.ordered[, !(colnames(res.deconv.epic.ordered) %in% patient.outlier)]
res.deconv.quantiseq.ordered <- res.deconv.quantiseq.ordered[, !(colnames(res.deconv.quantiseq.ordered) %in% patient.outlier)]

# Define shared cell types
shared_cell_types <- c("B cell", "T cell CD4+", "T cell CD8+", "NK cell")

# Check if all shared cell types exist in both datasets
available_epic <- intersect(shared_cell_types, rownames(res.deconv.epic.ordered))
available_quantiseq <- intersect(shared_cell_types, rownames(res.deconv.quantiseq.ordered))
cell_types_to_analyze <- intersect(available_epic, available_quantiseq)

print(paste("Cell types available for analysis:", paste(cell_types_to_analyze, collapse = ", ")))

# ==============================================================================
# Clustering and Survival Analysis (Pipeline from clustering_TCGA.R)
# ==============================================================================

# Load additional libraries
library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(dplyr)

# Load color palette
source("color_palette.R")

# Function to perform clustering and survival analysis
run_clustering_survival <- function(res.deconv, method_name, df.clinical) {
  # Scale: z-score transformation
  res.deconv.t <- t(res.deconv)
  res.mean <- colMeans(res.deconv.t)
  res.sd <- apply(res.deconv.t, 2, sd)
  res.deconv.scale <- t(scale(res.deconv.t, center = res.mean, scale = res.sd))
  
  # Output directory
  out_dir <- paste0("../output/comparison_", method_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Clustering
  # Using the same parameters as in clustering_TCGA.R
  results <- ConsensusClusterPlus(res.deconv.scale, maxK=6, reps=1000, pItem=0.8, pFeature=1, 
                                  title=out_dir, 
                                  clusterAlg="pam", distance="pearson", seed=123456, plot="pdf")
  
  # Get clusters (k=3)
  cluster <- unname(results[[3]][["consensusClass"]])
  
  # Prepare data for survival analysis
  res.deconv.sample <- as.data.frame(t(res.deconv.scale))
  res.deconv.sample$cluster <- as.factor(cluster)
  
  # Merge with clinical data
  # Ensure df.clinical has bcr_months
  if(!"bcr_months" %in% colnames(df.clinical)) {
     df.clinical$bcr_months <- df.clinical$bcr_days/30
  }
  
  df.clinical.cluster <- merge(df.clinical, res.deconv.sample, by = "row.names")
  df.clinical.cluster <- column_to_rownames(df.clinical.cluster, colnames(df.clinical.cluster)[1])
  
  df.clinical.cluster$bcr_months <- as.numeric(df.clinical.cluster$bcr_months)
  df.clinical.cluster$biochemical_recurrence_num <- 1*(df.clinical.cluster$biochemical_recurrence=="YES")
  
  # Survival Analysis
  km.fit <- surv_fit(Surv(bcr_months, biochemical_recurrence_num) ~ cluster, data=df.clinical.cluster)
  
  # Plot
  p.surv <- ggsurvplot(km.fit, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE,
                       title = paste("Survival Analysis -", method_name),
                       # palette = color.three.clusters, # Use if clusters are 3 and mapped
                       xlab = "Time in months", ylab = "BCR free probability")
                       
  pdf(file.path(out_dir, paste0("surv_curv_3_groups_", method_name, ".pdf")))
  print(p.surv)
  dev.off()
  
  return(list(results=results, plot=p.surv, data=df.clinical.cluster))
}

# Run for EPIC
print("Running analysis for EPIC...")
res.epic <- run_clustering_survival(res.deconv.epic.ordered, "epic", df.clinical.bcr)

# Run for quanTIseq
print("Running analysis for quanTIseq...")
res.quantiseq <- run_clustering_survival(res.deconv.quantiseq.ordered, "quantiseq", df.clinical.bcr)

# ==============================================================================
# Venn Diagram Analysis
# ==============================================================================
library(VennDiagram)

# Function to label clusters based on survival (Good, Intermediate, Bad)
label_clusters <- function(df) {
  # Calculate survival probability at median follow-up
  median_time <- median(df$bcr_months, na.rm = TRUE)
  
  fit <- survfit(Surv(bcr_months, biochemical_recurrence_num) ~ cluster, data = df)
  
  # Get survival probability at median time
  surv_summary <- summary(fit, times = median_time)
  
  # Create a dataframe of survival probabilities
  surv_df <- data.frame(
    cluster = gsub("cluster=", "", surv_summary$strata),
    surv_prob = surv_summary$surv
  )
  
  # Order clusters by survival probability (High -> Good, Low -> Bad)
  surv_df <- surv_df[order(surv_df$surv_prob, decreasing = TRUE), ]
  
  # Assign labels
  labels <- c("Good", "Intermediate", "Bad")
  
  # Create mapping
  cluster_mapping <- setNames(labels, surv_df$cluster)
  
  # Apply mapping to dataframe
  df$cluster_label <- cluster_mapping[as.character(df$cluster)]
  
  return(df)
}

# Label clusters for both methods
print("Labeling clusters based on survival...")
df.epic <- label_clusters(res.epic$data)
df.quantiseq <- label_clusters(res.quantiseq$data)

# Create output directory for Venn diagrams
venn_dir <- "../output/Venn_comparison_EPIC_quanTIseq"
dir.create(venn_dir, showWarnings = FALSE, recursive = TRUE)

# Function to draw Venn diagram for a specific label
draw_venn <- function(label, df1, df2, name1="EPIC", name2="quanTIseq") {
  patients1 <- rownames(df1[df1$cluster_label == label, ])
  patients2 <- rownames(df2[df2$cluster_label == label, ])
  
  venn.plot <- venn.diagram(
    x = list(patients1, patients2),
    category.names = c(name1, name2),
    filename = NULL,
    output = TRUE,
    fill = c("#E64B35", "#4BBAD3"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.5,
    main = paste("Overlap of", label, "Prognosis Cluster"),
    disable.logging = TRUE
  )
  
  pdf(file.path(venn_dir, paste0("Venn_", label, ".pdf")))
  grid.draw(venn.plot)
  dev.off()
}

# Draw Venn diagrams for each group
print("Generating Venn diagrams...")
draw_venn("Good", df.epic, df.quantiseq)
draw_venn("Intermediate", df.epic, df.quantiseq)
draw_venn("Bad", df.epic, df.quantiseq)

print("Venn diagrams generated.")

# ==============================================================================
# Additional Analysis: Robustness Check
# ==============================================================================
# 1. Correlation of shared cell types
# 2. Boxplots of quanTIseq estimates grouped by EPIC clusters
# 3. Survival analysis of "Consensus" patients

out_dir_robust <- "../output/Robustness_Check_EPIC_quanTIseq"
dir.create(out_dir_robust, showWarnings = FALSE, recursive = TRUE)

# --- 1. Correlation of shared cell types ---
print("Calculating correlation of shared cell types...")
pdf(file.path(out_dir_robust, "Correlation_Shared_Cell_Types.pdf"))
par(mfrow=c(2,2))
for (cell in cell_types_to_analyze) {
  val_epic <- as.numeric(res.deconv.epic.ordered[cell, ])
  val_quantiseq <- as.numeric(res.deconv.quantiseq.ordered[cell, ])
  
  cor_val <- cor(val_epic, val_quantiseq, method = "pearson")
  
  plot(val_epic, val_quantiseq, 
       main = paste0(cell, "\nPearson r = ", round(cor_val, 2)),
       xlab = "EPIC", ylab = "quanTIseq", pch = 19, col = rgb(0,0,0,0.2))
  abline(lm(val_quantiseq ~ val_epic), col = "red")
}
dev.off()

# --- 2. Boxplots of quanTIseq and EPIC estimates grouped by EPIC clusters ---
print("Generating boxplots of quanTIseq and EPIC estimates by EPIC clusters...")
# We want to see if EPIC clusters (Good/Bad) show distinct profiles in quanTIseq data
# Merge EPIC cluster info into quanTIseq data
# df.epic has 'cluster_label' (Good, Intermediate, Bad)
epic_clusters <- df.epic[, c("cluster_label"), drop=FALSE]
colnames(epic_clusters) <- "EPIC_Cluster"

# --- Boxplots for quanTIseq ---
# Transpose quanTIseq data for plotting
quantiseq_data_t <- as.data.frame(t(res.deconv.quantiseq.ordered))
quantiseq_data_t <- merge(quantiseq_data_t, epic_clusters, by = "row.names")

# Fix column names to be valid R names (remove spaces, +, etc.) to avoid ggpubr errors
# Save original names for titles
original_names <- colnames(quantiseq_data_t)
colnames(quantiseq_data_t) <- make.names(colnames(quantiseq_data_t))

# Identify the cell type columns (excluding Row.names and EPIC_Cluster)
cell_cols <- setdiff(colnames(quantiseq_data_t), c("Row.names", "EPIC_Cluster"))

# Plotting
library(ggpubr)
library(gridExtra) # Ensure gridExtra is loaded for grid.arrange

pdf(file.path(out_dir_robust, "Boxplots_quanTIseq_by_EPIC_Clusters.pdf"), width = 10, height = 8)
plot_list <- list()

for (cell in cell_cols) {
  # Find original name for title (optional, but nicer)
  # This is a bit tricky since make.names is not 1-to-1 reversible easily, 
  # but we can just use the valid name for the plot variable, and clean it up for the title.
  title_name <- gsub("\\.", " ", cell) 
  
  p <- ggboxplot(quantiseq_data_t, x = "EPIC_Cluster", y = cell,
                 color = "EPIC_Cluster", palette = c("#E64B35", "#4BBAD3", "#B3B3B3"),
                 add = "jitter", title = paste("quanTIseq:", title_name)) +
       stat_compare_means(label.y.npc = "top", label.x.npc = "center") + # Add p-value
       stat_summary(fun = median, geom = "text", aes(label = round(..y.., 3)), 
                    vjust = -0.5, color = "black", position = position_dodge(0.8)) + # Add median label
       theme(legend.position = "none") +
       ylab("Cell Fraction")
  plot_list[[cell]] <- p
}

# Arrange plots in a grid
# Split into pages if too many cells
if (length(plot_list) > 0) {
  # Print 4 plots per page
  num_plots <- length(plot_list)
  num_pages <- ceiling(num_plots / 4)
  
  for (i in 1:num_pages) {
    start_idx <- (i - 1) * 4 + 1
    end_idx <- min(i * 4, num_plots)
    
    # Extract sublist of plots
    # Note: plot_list is a named list, we need to access by index for slicing
    current_plots <- plot_list[start_idx:end_idx]
    
    do.call(grid.arrange, c(current_plots, ncol = 2, nrow = 2))
  }
}
dev.off()

# --- Boxplots for EPIC ---
# Transpose EPIC data for plotting
epic_data_t <- as.data.frame(t(res.deconv.epic.ordered))
epic_data_t <- merge(epic_data_t, epic_clusters, by = "row.names")

# Fix column names to be valid R names
colnames(epic_data_t) <- make.names(colnames(epic_data_t))

# Identify the cell type columns (excluding Row.names and EPIC_Cluster)
cell_cols_epic <- setdiff(colnames(epic_data_t), c("Row.names", "EPIC_Cluster"))

pdf(file.path(out_dir_robust, "Boxplots_EPIC_by_EPIC_Clusters.pdf"), width = 10, height = 8)
plot_list_epic <- list()

for (cell in cell_cols_epic) {
  title_name <- gsub("\\.", " ", cell) 
  
  p <- ggboxplot(epic_data_t, x = "EPIC_Cluster", y = cell,
                 color = "EPIC_Cluster", palette = c("#E64B35", "#4BBAD3", "#B3B3B3"),
                 add = "jitter", title = paste("EPIC:", title_name)) +
       stat_compare_means(label.y.npc = "top", label.x.npc = "center") + # Add p-value
       stat_summary(fun = median, geom = "text", aes(label = round(..y.., 3)), 
                    vjust = -0.5, color = "black", position = position_dodge(0.8)) + # Add median label
       theme(legend.position = "none") +
       ylab("Cell Fraction")
  plot_list_epic[[cell]] <- p
}

# Arrange plots in a grid
if (length(plot_list_epic) > 0) {
  num_plots <- length(plot_list_epic)
  num_pages <- ceiling(num_plots / 4)
  
  for (i in 1:num_pages) {
    start_idx <- (i - 1) * 4 + 1
    end_idx <- min(i * 4, num_plots)
    current_plots <- plot_list_epic[start_idx:end_idx]
    do.call(grid.arrange, c(current_plots, ncol = 2, nrow = 2))
  }
}
dev.off()

# --- 3. Survival analysis of "Consensus" patients ---
print("Performing survival analysis on consensus patients...")
# Identify patients who are consistently "Good" or consistently "Bad"
patients_consensus_good <- rownames(df.epic[df.epic$cluster_label == "Good" & df.quantiseq$cluster_label == "Good", ])
patients_consensus_intermediate <- rownames(df.epic[df.epic$cluster_label == "Intermediate" & df.quantiseq$cluster_label == "Intermediate", ])
patients_consensus_bad <- rownames(df.epic[df.epic$cluster_label == "Bad" & df.quantiseq$cluster_label == "Bad", ])

# Create a new dataframe for consensus analysis
df_consensus <- df.clinical.bcr[c(patients_consensus_good, patients_consensus_intermediate, patients_consensus_bad), ]
df_consensus$Consensus_Cluster <- ifelse(rownames(df_consensus) %in% patients_consensus_good, "Good",
                                        ifelse(rownames(df_consensus) %in% patients_consensus_intermediate, "Intermediate", "Bad"))

# Ensure bcr_months exists
if(!"bcr_months" %in% colnames(df_consensus)) {
  df_consensus$bcr_months <- df_consensus$bcr_days/30
}
df_consensus$bcr_months <- as.numeric(df_consensus$bcr_months)
df_consensus$biochemical_recurrence_num <- 1*(df_consensus$biochemical_recurrence=="YES")

# Set factor levels to ensure correct ordering in plot and table
df_consensus$Consensus_Cluster <- factor(df_consensus$Consensus_Cluster, levels = c("Good", "Intermediate", "Bad"))

# Survival Plot
km.fit.consensus <- surv_fit(Surv(bcr_months, biochemical_recurrence_num) ~ Consensus_Cluster, data=df_consensus)
p.surv.consensus <- ggsurvplot(km.fit.consensus, pval = TRUE, risk.table = TRUE, ncensor.plot = FALSE,
                     title = "Survival of Consensus Patients (EPIC & quanTIseq)",
                     palette = color.three.clusters,
                     xlab = "Time in months", ylab = "BCR free probability")

pdf(file.path(out_dir_robust, "Survival_Consensus_Patients.pdf"),
    width = 10,      # inches; try 12–16 depending on your plot
    height = 8)
print(p.surv.consensus)
dev.off()

print("Robustness analysis complete.")

print("Analysis complete. Check output directories for results.")

