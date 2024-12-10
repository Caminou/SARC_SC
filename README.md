# Single-cell RNAseq Analysis on Sarcoma


Full dataset is divided into SARC data (Tumor (n=5) + PDSC (n=6)) and Raw data from Tumor samples from Gruel et al (n=31) under [GSE221494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221494).

Both datasets have been processed, QC and analyzed in the same way:

#####  1) Data preparation and Quality control
- [x] [Wrangling and Demultiplexing](https://github.com/Caminou/SARC_SC/Load_Seurat/)
- [x] [QC](https://caminou.github.io/SARC_SC/QC/QC.html)
- [X] [Ambient RNA removal](https://github.com/Caminou/SARC_SC/blob/main/Ambient_RNA/)
- [X] [Doublet annotation](https://github.com/Caminou/SARC_SC/blob/main/Doublet_removal)
##### 2) Data Integration
- [X] [Tumor biopsies](https://github.com/Caminou/SARC_SC/tree/main/Integration/Tumor)
  - [x] Harmony
  - [x] SCVI
  - [x] Wnn
- [X] [Gruel et al](https://github.com/Caminou/SARC_SC/tree/main/Integration/Gruel)
  - [x] Harmony
  - [X] SCVI
  - [x] Wnn
- [X] Tumor and Gruel together
  - [x] Harmony
  - [X] SCVI
  - [x] Wnn
##### 4) Cell-type prediction
- [X] [Azimuth](https://github.com/Caminou/SARC_SC/blob/main/Integration/Integration_Tumor.Rmd#L64)
- [X] SCimilarity from HCA
##### 5) Clean cell-type annotation
- [X] Done, removal of cluster doublets and other cells
##### 6) CNV for Tumor-cell identification
- [ ] Copykat
- [ ] SCEVAN
- [ ] InferCNV
- [ ] ScATOMIC
##### 7) PDSC mapping and label transfer
- [ ] Seurat label transfer
- [ ] Seurat mapping
##### 8) Grouping tumor cells by MDM2 expression
- [ ] Gene signatures
- [ ] Cell proportion among patients 
##### 8) DGE :tada:
- [ ] Tumor biopsy vs PDSC
- [ ] ASPC vs tumor cells
- [ ] MDM2 high vs MDM2 low
##### 9) Cell-Cell interaction analysis
- [ ] CellPhone
##### 10) fGWAS
##### 11) Spatial validation?


