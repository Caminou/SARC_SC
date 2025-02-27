###### scATOMIC
library(scATOMIC)
library(plyr)
library(dplyr)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(Matrix)
library(Seurat)
library(agrmt)
library(cutoff.scATOMIC)
library(copykat)
library(ggplot2)
library(scATOMIC)

# Set the Python path in R
use_condaenv("r-reticulate", required = TRUE)
magic <- import("magic", convert=FALSE)

setwd("~/camino_sandbox/R_pipelines/CNV/")
output_dir <- "scatomic_PDSC"
dir.create(output_dir, showWarnings = FALSE)


for (ident in (unique(PDSC@meta.data$orig.ident))) {
  message("Reading...", ident)
  setwd("~/camino_sandbox/R_pipelines/CNV/scatomic_PDSC/")
  # Subset the Seurat object by orig.ident
  seurat_subset <- subset(PDSC, subset = orig.ident == ident)
  # Extract the RNA count matrix (as a sparse matrix)
  ident_dir <- file.path("~/camino_sandbox/R_pipelines/CNV/scatomic_PDSC", ident)
  dir.create(ident_dir, showWarnings = FALSE)
  setwd(ident_dir)
  rna_matrix <- (seurat_subset@assays$RNA$counts)
  cell_predictions  = run_scATOMIC(rna_matrix, mc.core=100)
  results_scatomic<- create_summary_matrix(prediction_list = cell_predictions, use_CNVs = F, modify_results = T, mc.cores = 100, raw_counts = rna_matrix, min_prop = 0.5, known_cancer_type = "sarcoma")
  # Create a directory for the current orig.ident
  
  # Save the results to the corresponding directory
  # Adjust the file name and format as needed
  saveRDS(results_scatomic, file = file.path(paste0(ident, "_scatomic_results.rds")))
  message("Finished...", ident)
  rm(rna_matrix, seurat_subset, results_scatomic, cell_predictions)
  gc()
}

