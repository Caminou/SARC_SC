library(Seurat)
library(dplyr)
library(future)
options(future.globals.maxSize=40000*1024^2)


### Load PDSC data and add doublet metadata ###

PDSC = readRDS("~/camino_sandbox/R_pipelines/SARC_f_Soup_doublet/SARC_PDSC_f_Soup_doublet.rds")

data= read.table("~/camino_sandbox/SC_linux/SARC_PDSC_DoubletFinder.csv")
metadata = PDSC@meta.data
metadata$Cell_barcodes=rownames(metadata)
data$Cell_barcodes=rownames(data)
metadata <- merge(metadata, data, all = TRUE, by="Cell_barcodes")
rownames(metadata)=metadata$Cell_barcodes
metadata = metadata[Cells(PDSC),]
table(metadata$scDblFinder.class, metadata$doublet_finder)

PDSC@meta.data = metadata
table(PDSC$scDblFinder.class, PDSC$doublet_finder)

PDSC_f = subset(PDSC, scDblFinder.class %in% c("doublet") & doublet_finder %in% c("Doublet"), invert=T)
table(PDSC_f$scDblFinder.class, PDSC_f$doublet_finder)


############ Label transfer from Biopsies ###########################

### Load Biopsy data and add all metadata

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

DimPlot(Biopsy, reduction="umap.harmony", raster=F, group.by = "cell_type")

PDSC_f <- subset(PDSC_f, features = rownames(Biopsy))

Biopsy@reductions[["umap.harmony"]]@misc[["model"]][["a"]] =0.9921756  #after rerunning object with return.model=T
Biopsy@reductions[["umap.harmony"]]@misc[["model"]][["b"]] = 1.112253  #after rerunning object with return.model=T
        
        ## Label transfer requires only 1 integrated SCT model, thus running SCT in Joined RNA assay
        ## Up to 450GB RAM required!!!

        Merge[["RNA"]] <- JoinLayers(Merge_data[["RNA"]]) 
        Merge = SCTransform(Merge, assay = "RNA")
        Merge_s <- subset(Merge, cells = Cells(Biopsy))
        metadata = as.data.frame(Biopsy@meta.data)
        metadata <- metadata[rownames(metadata) %in% Cells(Merge_s), ]
        Merge_s@meta.data = metadata

## First, predict the cell_type labels
PDSC_f <- JoinLayers(PDSC_f, assay="RNA")
PDSC_f <- SCTransform(PDSC_f, assay = "RNA")
PDSC_f <- RunPCA(PDSC_f)
DefaultAssay(PDSC_f)="SCT"
DefaultAssay(Merge)="SCT"

## Get anchors for the PCA 
anchors <- FindTransferAnchors(reference = Merge, query = PDSC_f, dims = 1:50, normalization.method="SCT",reference.reduction = "pca")

PDSC_f <- IntegrateEmbeddings(anchorset = anchors,reference = Merge,query = PDSC_f,new.reduction.name = "ref.pca")
PDSC_f <- ProjectUMAP(query = PDSC_f,
       query.reduction = "ref.pca",
       reference = Merge,
       reference.reduction = "pca",
       reduction.model = "umap.harmony")
DimPlot(querPDSC_fy_f, reduction="ref.umap", group.by="predicted.cell_type")

# Cell_type prediction
DefaultAssay(PDSC_f)="SCT"
DefaultAssay(Merge_s)="SCT"
anchors <- FindTransferAnchors(reference = Merge_s, query = PDSC_f, dims = 1:30,reference.reduction = "pca")
# Cell_type prediction
PDSC_f <- TransferData(
  anchorset = anchors,
  reference = Merge_s,
  query = PDSC_f,
  refdata = list(cell_type = "cell_type")
)
# Assuming PDSC_f@meta.data is a data frame
split_data <- do.call(rbind, strsplit(PDSC_f@meta.data[["orig.ident"]], "_"))
PDSC_f$patient <- split_data[, 1]  # First part is the patient
PDSC_f$passage <- split_data[, 2]  # Second part is the passage

DimPlot(Biopsy, group.by="cell_type", reduction = "umap.harmony", raster=F) + DimPlot(PDSC_f, group.by="predicted.cell_type", reduction="ref.umap", raster=F)
DimPlot(PDSC_f, group.by="predicted.cell_type", reduction = "ref.umap", raster=F, split.by="passage")
table(PDSC_f$predicted.id, PDSC_f$passage)
DefaultAssay(PDSC_f)="RNA"
FeaturePlot(PDSC_f, features=c("COL1A1","TAGLN","DCN"),reduction = "ref.umap", raster=F)
saveRDS(PDSC_f, "~/camino_sandbox/R_pipelines/SARC_f_Soup_doublet/PDSC_f_embedding.rds")
predicted_cell_type = as.data.frame(PDSC_f$predicted.cell_type)
write.csv(predicted_cell_type, "~/camino_sandbox/R_pipelines/SARC_f_Soup_doublet/predicted_cell_type_for_PDSC.csv")
