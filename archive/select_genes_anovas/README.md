# Archived: select-gene ANOVAs

These files are not used in the manuscript and are kept for the record only.

`select_genes_anovas.R` ran separate ANOVAs on FPKM values for five genes (Hph, sima, tgo, Vhl, ERR). Testing a handful of genes outside the genome-wide model is not statistically sound, so these results were discarded. Its input file, `hypoxia_select_genes.csv`, was made by hand.

The statistics reported for these genes in the paper (Figure S6 and Table S9) come from the genome-wide DESeq2 analysis performed by Novogene. See `code/figure_S6_HIF_genes.R` and `code/table_S9_hypoxia_genes.R`.
