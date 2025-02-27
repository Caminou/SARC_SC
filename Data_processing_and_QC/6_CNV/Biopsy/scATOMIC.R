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
output_dir <- "scATOMIC_outputs"
dir.create(output_dir, showWarnings = FALSE)

tumor_ASPC= readRDS("~/camino_sandbox/R_pipelines/ASPC_tcells_endot.rds")
tumor_ASPC = subset(tumor_ASPC, Sarcoma_type %in% c("lipoma"), invert=T)
DimPlot(tumor_ASPC, reduction="umap.harmony", group.by="cell_type")

for (ident in (unique(tumor_ASPC@meta.data$orig.ident))) {
  message("Reading...", ident)
  setwd("~/camino_sandbox/R_pipelines/CNV/")
  # Subset the Seurat object by orig.ident
  seurat_subset <- subset(tumor_ASPC, subset = orig.ident == ident)
  # Extract the RNA count matrix (as a sparse matrix)
  ident_dir <- file.path(output_dir, ident)
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


### Merge all results

base_path = c("~/camino_sandbox/R_pipelines/CNV/scATOMIC_outputs/")
file_list <- list.dirs(path = base_path, full.names = TRUE)

# Initialize an empty data frame to store results
scatomic_results= data.frame()

# Loop through each directory in file_list
for (i in file_list) {
  # Get all .rds files in the current directory
  rds_files <- list.files(path = i, pattern = "\\.rds$", full.names = TRUE)
  
  # Loop through each .rds file found in the directory
  for (file_path in rds_files) {
    # Read the file
    data <- readRDS(file_path)  # No need for header in readRDS
    
    # Remove unwanted columns
    data$nCount_RNA <- NULL
    data$nFeature_RNA <- NULL
    
    # Combine the data with the existing results
    scatomic_results <- dplyr::bind_rows(scatomic_results, data)
  }
  write.csv(scatomic_results, paste0(base_path, "scatomic_results.csv"), row.names = FALSE)
}

scatomic_results = read.csv("~/camino_sandbox/R_pipelines/CNV/scATOMIC_outputs/scatomic_results.csv")

