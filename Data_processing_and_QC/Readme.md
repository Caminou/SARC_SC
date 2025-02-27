##  Data preparation and Quality control
1) Wrangling and Demultiplexing using [HTOdemux](https://rdrr.io/github/satijalab/seurat/man/HTODemux.html)
2) QC filtering of nCount and nFeature (0.01 and 0.99 quantiles) and MT percentage > 20% by Run (orig.ident)
3) Removal of Ambient RNA using [SoupX](https://github.com/constantAmateur/SoupX) based on Sambomics [SoupX notebook](https://github.com/mousepixels/sanbomics_scripts/tree/main/soupX)
    - Filtered matrices from CellRanger output after quantile and mt.percent QC are normalized and a basic clustering by orig.ident is produced in order to improve RNA ambient capture
    - Raw feature barcode matrices from CR output are combined with the filtered matrices from CR creating a soup.channel object (table_of_droplets, table_of_counts)
    - Counts are adjusted leading to the reduction of initial filtered matrix between 0.01 to 0.2.
  Doublet characterization using [Doublet Finder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) and [scDblFinder](https://github.com/plger/scDblFinder) --> to note, doublets will not be removed for clustering however will be filtered out later on!!
4) Integration testing
  a) Gruel by itself
    - Harmony
    - scVI
    - Mnn
  b) SARC_tumor by itself
    - Harmony
    - scVI
    - Mnn
  c) Gruel and SARC_tumor together
    - Harmony
    - scVI
    - Mnn
5) Cell annotation and subclustering