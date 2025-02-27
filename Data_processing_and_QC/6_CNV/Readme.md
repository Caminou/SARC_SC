
# Cancer cell prediction based on CNV inference
In order to determine and extract Cancer cells with "high-confidence" from initial cell (cluster-based) annotation, three distinct CNV variation algorithms have been used, using Endothelial cells as the baseline CNV status
1. **Copykat**: CNV inference based, Endothelial cells used as CNV baseline
2. **SCEVAN**: CNV inference based, Endothelial cells used as CNV baseline
3. **scATOMIC**: Reference transcriptome based, no baseline but parameter of known cancer set to "sarcoma"
