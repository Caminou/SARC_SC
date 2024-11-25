# Single-cell RNAseq Analysis on Sarcoma


Full dataset is divided into SARC data (Tumor (n=5) + PDSC (n=6)) and Raw data from Tumor samples from Gruel et al (n=31) under [GSE221494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221494).

Both datasets have been processed, QC and analyzed in the same way:

#####  1) Data preparation and Quality control
- [x] [Wrangling and Demultiplexing](https://github.com/Caminou/SRC_SC/tree/main/QC/Load_wrangling_and_QC.Rmd)
- [x] [QC](https://github.com/Caminou/SRC_SC/tree/main/QC)
- [X] [Ambient RNA removal](https://github.com/Caminou/SARC_SC/blob/main/Ambient_RNA/)
- [X] [Doublet removal](https://github.com/Caminou/SARC_SC/blob/main/Doublet_removal)
##### 2) Data Integration
- [ ] Tumor biopsies
  - [x] [Harmony](https://github.com/Caminou/SARC_SC/blob/main/Integration/Azimuth_and_integration.Rmd#L124)
  - [ ] SCVI
  - [x] Wnn
- [ ] Gruel et al
  - [x] Harmony
  - [ ] SCVI
  - [x] Wnn
##### 4) Cell-type prediction
- [X] Azimuth
- [ ] SCimilarity from HCA
- [ ] Gruel et al vs SARC
##### 5) Clean cell-type annotation
##### 6) PDSC mapping and label transfer
- [ ] Seurat label transfer
- [ ] Seurat mapping
##### 7) CNV for Tumor-cell identification
- [ ] Copykat
- [ ] SCEVAN
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


