
rm(list=ls())

library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tibble)

# Create output directory
out_dir <- "../output/dea_gsea"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Load data
message("Loading data...")
load("../data/clustering_tcga_removing_badguys.RData") # Contains df.clinical.bcr.cluster
load("../data/prad_data_tcga.RData") # Contains prad.gene

# Align samples
# df.clinical.bcr.cluster rownames are patient IDs
# prad.gene colnames are patient IDs
common_samples <- intersect(rownames(df.clinical.bcr.cluster), colnames(prad.gene))
message(paste("Number of common samples:", length(common_samples)))

clinical_data <- df.clinical.bcr.cluster[common_samples, ]
expression_data <- prad.gene[, common_samples]

# Check if expression data needs log transformation
# Assuming RNASeq2GeneNorm data (RSEM), it's often not log-transformed.
# We usually do log2(x + 1) for limma if it's not counts.
# Let's check max value to guess.
if (max(expression_data, na.rm = TRUE) > 100) {
  message("Applying log2(x+1) transformation...")
  expression_data <- log2(expression_data + 1)
}

# Filter out non-expressed genes or genes with zero variance
# This prevents the "Zero sample variances detected" warning and improves power
message("Filtering low-expressed/invariant genes...")
# Keep genes expressed (value > 0) in at least 10 samples
keep_expressed <- rowSums(expression_data > 0) >= 10
# Also ensure variance is not zero
vars <- apply(expression_data, 1, var)
keep_var <- vars > 0

keep <- keep_expressed & keep_var
expression_data <- expression_data[keep, ]
message(paste("Retained", nrow(expression_data), "genes for analysis."))

# Prepare design matrix
group <- factor(clinical_data$cluster, levels = c("TCE", "EPCE", "TASCE"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit linear model
message("Fitting linear model...")
fit <- lmFit(expression_data, design)

# Define contrasts
# Comparing against TCE (Good) and between others
contrasts <- makeContrasts(
  EPCE_vs_TCE = EPCE - TCE,
  TASCE_vs_TCE = TASCE - TCE,
  TASCE_vs_EPCE = TASCE - EPCE,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Function to perform GSEA and save results
run_dea_gsea <- function(fit_obj, contrast_name) {
  message(paste("Processing contrast:", contrast_name))
  
  # --- DEA ---
  res_table <- topTable(fit_obj, coef = contrast_name, number = Inf)
  res_table$gene <- rownames(res_table)
  
  # Save DEA results
  write.csv(res_table, file.path(out_dir, paste0("DEA_", contrast_name, ".csv")), row.names = FALSE)
  
  # --- GSEA ---
  # Prepare gene list for GSEA
  # We use t-statistic for ranking
  gene_list <- res_table$t
  names(gene_list) <- res_table$gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Convert gene symbols to Entrez IDs
  ids <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Remove duplicates if any
  ids <- ids[!duplicated(ids$SYMBOL), ]
  
  # Filter gene list to those with Entrez IDs
  gene_list_entrez <- gene_list[ids$SYMBOL]
  names(gene_list_entrez) <- ids$ENTREZID
  
  # Run GSEA (GO)
  message("Running GSEA (GO)...")
  gse_go <- gseGO(geneList = gene_list_entrez,
                  ont = "ALL",
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  pvalueCutoff = 0.05,
                  verbose = FALSE,
                  pAdjustMethod = "BH")
  
  if (!is.null(gse_go) && nrow(gse_go) > 0) {
    save(gse_go, file = file.path(out_dir, paste0("GSEA_GO_", contrast_name, ".RData")))
    write.csv(as.data.frame(gse_go), file.path(out_dir, paste0("GSEA_GO_", contrast_name, ".csv")))
    
    # Dotplot
    # Split by ontology for better visualization
    p1 <- dotplot(gse_go, showCategory=10, split="ONTOLOGY") + 
      facet_grid(ONTOLOGY~., scale="free") + 
      ggtitle(paste("GSEA GO:", contrast_name))
    ggsave(file.path(out_dir, paste0("GSEA_GO_dotplot_", contrast_name, ".pdf")), p1, width = 12, height = 12)
  } else {
    message("No significant GSEA GO results found.")
  }
  
  # Run GSEA (KEGG)
  message("Running GSEA (KEGG)...")
  gse_kegg <- gseKEGG(geneList = gene_list_entrez,
                      organism = 'hsa',
                      pvalueCutoff = 0.05,
                      verbose = FALSE)
  
  if (!is.null(gse_kegg) && nrow(gse_kegg) > 0) {
    save(gse_kegg, file = file.path(out_dir, paste0("GSEA_KEGG_", contrast_name, ".RData")))
    write.csv(as.data.frame(gse_kegg), file.path(out_dir, paste0("GSEA_KEGG_", contrast_name, ".csv")))
    
    # Dotplot
    p2 <- dotplot(gse_kegg, showCategory=20) + ggtitle(paste("GSEA KEGG:", contrast_name))
    ggsave(file.path(out_dir, paste0("GSEA_KEGG_dotplot_", contrast_name, ".pdf")), p2, width = 10, height = 8)
  } else {
    message("No significant GSEA KEGG results found.")
  }
}

# Run for all contrasts
contrast_names <- colnames(contrasts)
for (cn in contrast_names) {
  run_dea_gsea(fit2, cn)
}

message("Analysis complete. Results saved to ", out_dir)
