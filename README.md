# Single-cell RNAseq Analysis on Sarcoma

Single-cell RNAseq Analysis on Sarcoma
All samples have been processed from RAW .fastq files using cellranger_8.0.0 and GRCh38-2024-A.

## Experiment Design
Samples are divided into two different datasets:
- <b>SARC dataset (n=11)</b>
  - Tumor biopsies (n=5)
  - PDSC
      - passage 0 (n=3)
      - passage 4 (n=3)
- <b>Gruel dataset (n=31)</b>
  - Tumor biopsies available at [GSE221494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221494).

### 1. Accessing Jupyter
```bash
# SSH into the server with X or Y forwarding
ssh -Y user@threadripper.X.X
# activate micromamba shell
eval "$(micromamba shell hook --shell bash)"
micromamba activate r-reticulate
jupyter notebook
# Forward the local port 
ssh -L 8889:localhost:8889 user@threadripper.X.X
```
### 2. To reproduce in other servers:
To reproduce the environment on a different server, follow these steps:
1. Clone this repository:
   
```bash
git clone https://github.com/xingyiwoo/SARC_scRNAseq.git
```
2. Set up the environments with necessary packages:
```bash
# R environment (R must be version 4.1.1)
conda env create -f environment_r.yml r_env
conda activate r_env
          
# python environment
conda env create -f environment_py.yml
conda activate py_env
```
-------------------------------------------------------------

### Data_processing_and_QC
  - [X] [1_Wrangling and Demultiplexing](https://github.com/Caminou/SARC_SC/Data_processing_and_QC/1_Demultiplexing)
  - [x] [2_QC](https://caminou.github.io/SARC_SC/QC/QC.html)  
  - [X] [3_Ambient RNA removal and Doublet identification](https://github.com/Caminou/SARC_SC/Data_processing_and_QC/3_Ambient_RNA_doublet)
##### Data Integration
  - [X] [Tumor biopsies](https://github.com/Caminou/SARC_SC/Data_processing_and_QC/4_Integration)
    - [x] Harmony
    - [x] SCVI
    - [x] Wnn
  - [X] [Gruel et al](https://github.com/Caminou/SARC_SC/Data_processing_and_QC/4_Integration)
    - [x] Harmony
    - [X] SCVI
    - [x] Wnn
  - [X] [Tumor and Gruel together](https://github.com/Caminou/SARC_SC/Data_processing_and_QC/4_Integration)
    - [x] Harmony
    - [X] SCVI
    - [x] Wnn
##### [Cell-type prediction and clean-up](https://github.com/Caminou/SARC_SC/Data_processing_and_QC/5_Cell_annotation)
  - [X] Azimuth
  - [X] SCimilarity from HCA
##### 6) [CNV for Tumor-cell identification](https://github.com/Caminou/SARC_SC/Data_processing_and_QC/5_CNV)
  - [X] Copykat
  - [X] SCEVAN
  - [X] InferCNV
  - [X] ScATOMIC
##### 7) PDSC mapping and label transfer
  - [X] Seurat label transfer
  - [X] Seurat mapping
##### 7.1) CNV in PDSC cells using control Biopsy as control baseline
  - [X] Copykat
  - [X] SCEVAN
  - [X] InferCNV
  - [X] ScATOMIC
##### 8) Grouping tumor cells by MDM2 expression
  - [X] Gene signatures
  - [X] Cell proportion among patients 
##### 8) DGE :tada:
  - [X] Tumor biopsy vs PDSC
  - [X] ASPC vs tumor cells
  - [X] MDM2 high vs MDM2 low
##### 9) Cell-Cell interaction analysis
  - [ ] CellPhone
##### 10) fGWAS
##### 11) Spatial validation?


