
#################################  CNV pipeline using SCEVAN ##################################

Tumor_ASPC_T = readRDS("~/camino_sandbox/R_pipelines/ASPC_tcells_endot.rds")
Tumor_ASPC_T = subset(Tumor_ASPC_T, Sarcoma_type %in% c("lipoma"), invert=T)
Tumor_ASPC_T = subset(Tumor_ASPC_T, cell_type %in% c("ASPC","ASPC_preadipocyte","Mesothelial","Muscle","Endothelial","Adipocyte_C1QA"))
## SCEVAN
setwd("~/camino_sandbox/R_pipelines/CNV")
output_dir <- "SCEVAN_outputs"
dir.create(output_dir, showWarnings = FALSE)
library(SCEVAN)
library(Seurat)
DefaultAssay(Tumor_ASPC_T)="RNA"
for (ident in (unique(Tumor_ASPC_T@meta.data$orig.ident))) {
  message("Reading...", ident)
  setwd("~/camino_sandbox/R_pipelines/CNV/")
  # Subset the Seurat object by orig.ident
  seurat_subset <- subset(Tumor_ASPC_T, subset = orig.ident == ident)
  # Extract the RNA count matrix (as a sparse matrix)
  rna_matrix <- (seurat_subset@assays$RNA$counts)
  # Create a directory for the current orig.ident
  ident_dir <- file.path(output_dir, ident)
  dir.create(ident_dir, showWarnings = FALSE)
  setwd(ident_dir)
  norm_cells = Cells(subset(seurat_subset, cell_type %in% c("Endothelial")))
  # Run the copykat function
  scevan_results <- pipelineCNA(rna_matrix, par_cores = 70, norm_cell =norm_cells, SUBCLONES = FALSE, SCEVANsignatures = FALSE)
  # Save the results to the corresponding directory
  # Adjust the file name and format as needed
  saveRDS(scevan_results, file = file.path(paste0(ident, "_scevan_results.rds")))
  message("Finished...", ident)
  rm(scevan_results, rna_matrix, seurat_subset)
  gc()
}
