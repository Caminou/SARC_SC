
###### inferCNV ###########
library(infercnv)
library(Seurat)
setwd("~/camino_sandbox/R_pipelines/CNV")
options(bitmapType="Xlib") #set for figure visualization!
tumor_ASPC_T = readRDS("~/camino_sandbox/R_pipelines/CNV/SARC_Tumor_Gruel_ASPC_T.rds")

tumor_ASPC_T = subset(tumor_ASPC_T, Sarcoma_type %in% c("lipoma"), invert=T)
tumor_ASPC_T = subset(tumor_ASPC_T, cell_type %in% c("Endothelial","ASPC","ASPC_preadipocyte","Muscle","Mesothelial","Adipocyte_C1QA"))

## Read gene_position file
gene_post = read.table("~/camino_sandbox/R_pipelines/CNV/hg38_gencode_v27.txt", row.names = 1)
matrix_counts= GetAssayData(tumor_ASPC_T, layer="counts")
cellnames = colnames(matrix_counts)
annotations = as.data.frame(tumor_ASPC_T$cell_type)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=matrix_counts,
                                    annotations_file=annotations,
                                    delim="\t",
                                    gene_order_file=gene_post,
                                    ref_group_names=c("Endothelial"))
output_dir = c("InferCNV_outputs")

dir.create(output_dir, showWarnings = FALSE)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=file.path(output_dir), 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE, leiden_resolution = 0.01,
                             HMM=TRUE,HMM_type = "i6", up_to_step = 15)



# plot_cnv(infercnv_obj,
         out_dir=file.path(output_dir),
         obs_title="Observations (Cells)",
         ref_title="References (Cells)",
         cluster_by_groups=TRUE,
         x.center=1,
         x.range="auto",
         hclust_method='ward.D',
         color_safe_pal=FALSE,
         output_filename="infercnv",
         output_format="png",
         png_res=300,
         dynamic_resize=0
)

plot_subclusters(infercnv_obj,
                 out_dir=file.path(output_dir),
                 output_filename="subclusters_as_annotations")
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=file.path(output_dir), 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)



scores=apply(infercnv_obj@expr.data,2,function(x){ sum(x < 0.95 | x > 1.05)/length(x) })
