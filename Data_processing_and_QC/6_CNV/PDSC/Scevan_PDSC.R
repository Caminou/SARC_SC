library(SCEVAN)
library(Seurat)
library(magrittr)
library(glue)
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
Biopsy_endo = subset(Biopsy, orig.ident %in% c("hSC72_t","hSC81_t", "hSC92_t"))
Biopsy_endo = subset(Biopsy_endo, cell_type %in%c("Endothelial"))
Biopsy_endo[["SCT"]]=NULL
DefaultAssay(PDSC_f)="RNA"
PDSC_f[["SCT"]]=NULL
common.features <- intersect(rownames(Biopsy_endo), rownames(PDSC_f))
Merge = merge(Biopsy_endo[common.features,], PDSC_f[common.features,])
Merge = JoinLayers(Merge)
#run scevan
setwd("~/camino_sandbox/R_pipelines/CNV/SCEVAN_PDSC/")
output_dir <- "SCEVAN_outputs"
dir.create(output_dir, showWarnings = FALSE)
DefaultAssay(PDSC_f)="RNA"

for (ident in unique(Merge@meta.data$patient)) {
  message("Reading...", ident)
  setwd("~/camino_sandbox/R_pipelines/CNV/SCEVAN_PDSC/")
  # Subset the Seurat object by orig.ident
  seurat_subset <- subset(Merge, subset = patient == ident)
  rna_matrix = seurat_subset@assays$RNA$counts
  for (i in unique(seurat_subset$orig.ident)){
    if (grepl("_t", i)) {
      next  # Skips the current iteration and moves to the next one
    }
    else{
      message("Reading...", i)
      dir.create(paste0("~/camino_sandbox/R_pipelines/CNV/SCEVAN_PDSC/",i), showWarnings = FALSE)
      setwd(paste0("~/camino_sandbox/R_pipelines/CNV/SCEVAN_PDSC/",i))
      norm_cells <- Cells(subset(seurat_subset, cell_type %in% c("Endothelial")))
      passage_cells <- Cells(subset(seurat_subset, orig.ident %in% i))
      cell_chunks <- split(passage_cells, ceiling(seq_along(passage_cells) / (2000 - length(norm_cells))))
      # Process each chunk
      for (chunk_idx in seq_along(cell_chunks)) {
        chunk_cells <- cell_chunks[[chunk_idx]]
        all_cells = (c(norm_cells, chunk_cells))
        chunk_name <- paste0(ident, "_chunk_", chunk_idx)
        message("Processing chunk: ", chunk_name, ": ", length(chunk_cells),"/",ncol(seurat_subset))
        chunk_subset = subset(seurat_subset, cells = all_cells)
        # Extract the RNA count matrix (as a sparse matrix)
        rna_matrix <- chunk_subset@assays$RNA$counts
        
        # set new folder
        new_dir <- file.path(getwd(), chunk_name)
        
        # Create the directory (if it doesn't already exist)
        if (!dir.exists(new_dir)) {
          dir.create(new_dir)
        }
        
        # Set the new directory as the working directory
        setwd(new_dir)
        
        # Verify the working directory
        getwd()
        
        # Run the copykat function (or pipelineCNA in your case)
        scevan_results <- pipelineCNA(
          rna_matrix, 
          par_cores = 70, 
          norm_cell = norm_cells, 
          SUBCLONES = FALSE, 
          SCEVANsignatures = FALSE
        )
        saveRDS(scevan_results, file = file.path(paste0(chunk_name, "_scevan_results.rds")))
        message("Finished chunk: ", chunk_name)
        rm(scevan_results, rna_matrix, subset_seurat)
        gc()
      }
    }
  rm(scevan_results, rna_matrix, subset_seurat)
  gc()
  }
}



### write results for each chunk

files =list.files("~/camino_sandbox/R_pipelines/CNV/SCEVAN_PDSC/hSC81_p4/")
files <- files[grep(".rds", files)]
results = data.frame()

for (file in files) {
  message(paste0("Reading file:", file))
  data = readRDS(paste0("~/camino_sandbox/R_pipelines/CNV/SCEVAN_PDSC/hSC81_p4/",file))
  data$Cell_barcodes = rownames(data)
  results = rbind(results, data)
}

write.csv(results, "~/camino_sandbox/R_pipelines/CNV/SCEVAN_PDSC/hSC81_p4/merge_results.csv")

