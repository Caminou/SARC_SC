# RNA Soup removal
In order to distinguishing the actual cell expression profiles from background noise (ambient RNA) [SoupX](https://github.com/constantAmateur/SoupX).

This notebook is based on Sambomics [SoupX notebook] (https://github.com/mousepixels/sanbomics_scripts/tree/main/soupX)

### Steps:

1) Filtered matrices from CellRanger output after quantile and mt.percent QC are normalized and a basic clustering by orig.ident is produced in order to improve RNA ambient capture
2) Raw feature barcode matrices from CR output are combined with the filtered matrices from CR creating a soup.channel object (table_of_droplets, table_of_counts)
3) Counts are adjusted leading to the reduction of initial filtered matrix between 0.01 to 0.2.
