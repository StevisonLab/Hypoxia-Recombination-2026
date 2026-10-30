setwd("~/GitHub/Hypoxia-Recombination-2026/")

# ============================================================
# Table S9: expression of known hypoxia-responsive genes
#
# All statistics are taken from the genome-wide DESeq2 results
# (Novogene; Benjamini-Hochberg adjusted P) in inputs/*_deg.csv.
# No new tests are performed. Mean FPKM is across all 24 samples.
# ============================================================

library(dplyr)

genes <- read.csv(text = "
Category|Gene|FlyBase ID|Expected response to hypoxia|Reference
HIF pathway|sima|FBgn0266411|mRNA unchanged (protein stabilized)|Lavista-Llanos et al. 2002
HIF pathway|tgo|FBgn0264075|Constitutive|Romero et al. 2007
HIF pathway|Vhl|FBgn0041174|Regulates Sima stability|Romero et al. 2007
HIF pathway|fatiga (Hph)|FBgn0264785|Up (fgaB isoform)|Acevedo et al. 2010
HIF pathway|scylla|FBgn0041094|Up|Reiling & Hafen 2004; Liu et al. 2006; Li et al. 2013
HIF pathway|Tribbles|FBgn0028978|Up (Sima target)|Noguchi et al. 2022; Li et al. 2013
HIF pathway|Ldh (ImpL3)|FBgn0001258|Up|Li et al. 2013
HIF pathway|Pyk|FBgn0267385|Glycolytic enzyme; not hypoxia-induced in wild-type larvae|Li et al. 2013
HIF pathway|Musashi|FBgn0011666|mRNA unchanged (protein down)|Bertolin et al. 2016
FOXO-associated|Thor (4E-BP)|FBgn0261560|Up|Liu et al. 2006; Li et al. 2013
FOXO-associated|Hipk|FBgn0035142|Up|Ding et al. 2022
FOXO-associated|Relish|FBgn0014018|Up|Barretto et al. 2020
Metabolism|Pepck (CG10924)|FBgn0034356|Up|Mossman et al. 2017
Metabolism|hairy|FBgn0001168|Up|Mossman et al. 2017; Zhou et al. 2008
Metabolism|SdhC|FBgn0037873|Up in one female background|Mossman et al. 2017
Metabolism|SdhD|FBgn0039112|Up in one female background|Mossman et al. 2017
Acute screens|MFS14|FBgn0010651|Up|Azad et al. 2009; Mossman et al. 2017
Acute screens|CG31324|FBgn0051324|Up|Mossman et al. 2017
Acute screens|CG1600|FBgn0033188|Up|Azad et al. 2009; Liu et al. 2006
Acute screens|GstD2 (GstD21)|FBgn0010038|Up|Liu et al. 2006
Acute screens|astray|FBgn0023129|Up|Liu et al. 2006
Acute screens|ade5|FBgn0020513|Up|Liu et al. 2006
Acute screens|Mdr49|FBgn0004512|Up|Azad et al. 2009
Stress and autophagy|Hsp70Aa|FBgn0013275|Up|Azad et al. 2009
Stress and autophagy|Hsp70Ab|FBgn0013276|Up|Azad et al. 2009
Stress and autophagy|Hsp23|FBgn0001224|Up|Azad et al. 2009
Stress and autophagy|Atg8a|FBgn0052672|Up|Valko et al. 2022
", sep = "|", strip.white = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Each DEG file reports log2FoldChange = log2(first/second treatment in the
# file name). `sign` re-expresses it for the first-named treatment in the
# column heading (Hypoxia = Chamber 8 kPa; Normoxia = Chamber 21 kPa).
comparisons <- tribble(
  ~column,                       ~file,                          ~sign,
  "Pooled: 8 kPa vs Control",    "HypoxiavsControl",                1,
  "Pooled: 21 kPa vs Control",   "NormoxiavsControl",               1,
  "Pooled: 8 kPa vs 21 kPa",     "NormoxiavsHypoxia",              -1,
  "Day 2: 8 kPa vs Control",     "Control_Day2vsHypoxia_Day2",     -1,
  "Day 2: 21 kPa vs Control",    "Control_Day2vsNormoxia_Day2",    -1,
  "Day 2: 8 kPa vs 21 kPa",      "Hypoxia_Day2vsNormoxia_Day2",     1,
  "Day 6: 8 kPa vs Control",     "Control_Day6vsHypoxia_Day6",     -1,
  "Day 6: 21 kPa vs Control",    "Control_Day6vsNormoxia_Day6",    -1,
  "Day 6: 8 kPa vs 21 kPa",      "Hypoxia_Day6vsNormoxia_Day6",     1
)

read_deg <- function(file) {
  d <- read.csv(file.path("inputs", paste0(file, "_deg.csv")), check.names = FALSE,
                fileEncoding = "UTF-8-BOM")
  d[match(genes$`FlyBase ID`, d$gene_id), ]
}
degs <- lapply(comparisons$file, read_deg)
stopifnot(!anyNA(degs[[1]]$gene_id))

# Mean FPKM across all 24 samples: the six per-day files together contain
# every sample (each one twice), so keep each sample column once
fpkm <- do.call(cbind, lapply(degs[4:9], function(d) d[, grep("_fpkm$", names(d))]))
fpkm <- fpkm[, !duplicated(names(fpkm))]
stopifnot(ncol(fpkm) == 24)
genes$`Mean FPKM` <- round(rowMeans(fpkm), 1)

# Significant results as "+0.16 (padj=0.017)"; otherwise "ns"
for (i in seq_len(nrow(comparisons))) {
  l2fc <- comparisons$sign[i] * degs[[i]]$log2FoldChange
  padj <- degs[[i]]$padj
  genes[[comparisons$column[i]]] <- ifelse(!is.na(padj) & padj < 0.05,
    sprintf("%+.2f (padj=%s)", l2fc, as.character(signif(padj, 2))), "ns")
}

write.csv(genes, "outputs/table_S9_hypoxia_genes.csv", row.names = FALSE)
