
sample_dir <- "~/camino_sandbox/R_pipelines/CNV/Copykat/Copykat_outputs_Endot/"

library(tibble)
library(ggplot2)
library(reshape2)
library(dplyr)

# List of sample IDs to process
sample_ids <- unique(Tumor_ASPC_T$orig.ident.x)# Add all your sample IDs here

# Initialize an empty data frame to store results for all samples
all_chr12_data <- data.frame(
  chrom = numeric(),
  chrompos = numeric(),
  abspos = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each sample
for (sample_id in sample_ids) {
  message("Reading sample: ", sample_id)
  # Define file paths for the count matrix and CNA matrix
  data = read.table(paste0(sample_dir, sample_id,"/_copykat_CNA_results.txt"), header=T)
  chr12_data <- data %>% filter(chrom=="12") # Check if files exist
  rownames(chr12_data) = chr12_data$chrompos
  matrix_data <- as.matrix(chr12_data)
  # Combine with the overall data frame
  all_chr12_data <- all_chr12_data %>%
    full_join(chr12_data, by = c("chrom","chrompos","abspos"))
}


rownames(all_chr12_data) = all_chr12_data$chrompos



ordered_all_chr12_data = all_chr12_data %>% arrange(as.numeric(all_chr12_data$chrompos))


ordered_all_chr12_data$chrom <- NULL
ordered_all_chr12_data$abspos <- NULL
rownames(ordered_all_chr12_data) <- ordered_all_chr12_data$chrompos


# Step 2: Reshape data into long format
df_long <- melt(ordered_all_chr12_data, id.vars=c("chrompos"), variable.name = c("Cell"), value.name = "CN_ratio")

# Step 3: Add metadata information

df_long$Cell <- gsub("\\.", "-", df_long$Cell)

metadata = as.data.frame(Tumor_ASPC_T@meta.data)
df_long_with_metadata <- df_long %>%
  left_join(metadata, by = c("Cell" = "Cell_barcodes"))

# Step 4: Calculate mean CN ratio by position and cell type
df_cn_position_celltype <- df_long_with_metadata %>%
  group_by(chrompos, final_annotation, orig.ident.x) %>%
  summarize(Mean_CNR = mean(CN_ratio, na.rm = TRUE))

# Step 5: Plot mean CN ratio by position for each cell type
ggplot(df_cn_position_celltype, aes(x = chrompos, y = Mean_CNR, color=final_annotation)) +
  geom_line() +  # Line plot to show CN ratio changes
  geom_point(size = 0.1) +  # Points to highlight individual positions
  labs(x = "Chromosome 12 Position", y = "Mean CN Ratio", title = "CN Ratio along Chromosome 12 by Cell Type") +
  theme_minimal() +
  facet_wrap(~ orig.ident.x) # Separate plots for each cell type
