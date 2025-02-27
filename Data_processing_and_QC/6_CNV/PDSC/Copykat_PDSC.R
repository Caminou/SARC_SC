library(magrittr)
library(glue)
library(Seurat)
library(plyr)
library(dplyr)
library(data.table)
library(randomForest)
library(caret)
library(parallel)
library(reticulate)
library(Rmagic)
library(Matrix)
library(agrmt)
library(copykat)
library(ggplot2)
library(stringr)
library(future)
library(future.apply)
library(doParallel)


PDSC_f = readRDS("~/camino_sandbox/R_pipelines/SARC_f_Soup_doublet/PDSC_f_embedding.rds")
split_data <- do.call(rbind, strsplit(PDSC_f@meta.data[["orig.ident"]], "_"))
PDSC_f$patient <- split_data[, 1]  # First part is the patient
PDSC_f$passage <- split_data[, 2]  # Second part is the passage
info = read.csv("~/camino_sandbox/R_pipelines/SARC_f_Soup_doublet/predicted_cell_type_for_PDSC.csv", row.names = 1, header=T)
PDSC_f =AddMetaData(PDSC_f, info)
Biopsy = readRDS("~/camino_sandbox/R_pipelines/Merge_data_f_annotations_corrected_removal_doublets.rds")
all_cnv_results = read.csv("~/camino_sandbox/R_pipelines/CNV/CNV_final_results/all_cnv_results.csv")
all_cnv_results$X.1 = NULL
rownames(all_cnv_results)=all_cnv_results$Cell_barcodes

metadata = Biopsy@meta.data
# Filter out columns in all_cnv_results that are also in metadata
all_cnv_results <- all_cnv_results[, !colnames(all_cnv_results) %in% colnames(metadata)]
all_cnv_results$Cell_barcodes = all_cnv_results$cell_names
all_cnv_results$X.1 = NULL
all_cnv_results$cell_names = NULL

metadata= merge(metadata, all_cnv_results, by="Cell_barcodes", all=TRUE)
rownames(metadata)=metadata$Cell_barcodes
metadata = metadata[Cells(Biopsy),]
Biopsy@meta.data=metadata

## Running Endothelial as baseline
# Merge Biopsy and PDSC
Biopsy_endo = subset(Biopsy, orig.ident %in% c("hSC92_t"))
Biopsy_endo = subset(Biopsy_endo, cell_type %in%c("Endothelial"))
Biopsy_endo[["SCT"]]=NULL
DefaultAssay(PDSC_f)="RNA"
PDSC_f[["SCT"]]=NULL
common.features <- intersect(rownames(Biopsy_endo), rownames(PDSC_f))
Merge = merge(Biopsy_endo[common.features,], PDSC_f[common.features,])
Merge = JoinLayers(Merge)
#run copykat
Merge = subset(Merge, patient %in% c("hSC92"))
for (ident in (unique(Merge@meta.data$patient))) {
  message("Reading...", ident)
  setwd("~/camino_sandbox/R_pipelines/CNV/Copykat/Copykat_PDSC")
  output_dir="~/camino_sandbox/R_pipelines/CNV/Copykat/Copykat_PDSC"
  # Subset the Seurat object by orig.ident
  seurat_subset <- subset(Merge, subset = patient == ident)
  normal_cells = as.vector(Cells(subset(seurat_subset, cell_type %in%  c("Endothelial"))))
  # Extract the RNA count matrix (as a sparse matrix)
  rna_matrix <- (seurat_subset@assays$RNA$counts)
  # Create a directory for the current orig.ident
  ident_dir <- file.path(output_dir, ident)
  dir.create(ident_dir, showWarnings = FALSE)
  setwd(ident_dir)
  # Run the copykat function
  copykat_results = copykat(rna_matrix, norm.cell.names = normal_cells,genome="hg20",n.cores=5,id.type="S")
  # Save the results to the corresponding directory
  # Adjust the file name and format as needed
  saveRDS(copykat_results, file = file.path(paste0(ident, "_copykat_results.rds")))
  message("Finished...", ident)
  rm(copykat_results, rna_matrix, seurat_subset,normal_cells)
  gc()
}
 
