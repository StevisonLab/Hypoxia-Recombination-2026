setwd("~/GitHub/Hypoxia-Recombination-2026/")

# ============================================================
# Tables S10 and S11: KEGG enrichment and pathway-level fold changes for every
# pairwise comparison, pooled and by sampling day
#
# Table S10: KEGG enrichment from Novogene's clusterProfiler results
#     (inputs/*KEGGenrich*.csv and inputs/KEGG-day2,6/). No new tests.
# Table S11: for four KEGG pathways (ribosome, oxidative phosphorylation,
#     glycolysis/gluconeogenesis, TCA cycle), the median and interquartile
#     range of DESeq2 log2 fold changes across all pathway genes in the
#     DESeq2 results, and the number of genes with adjusted P < 0.05.
#     These are summaries only; no new tests.
#
# All comparisons are oriented as in the main text: 8 kPa vs Control,
# 21 kPa vs Control and 8 kPa vs 21 kPa, with direction and fold change
# for the first-named treatment. Novogene's file names use the opposite
# order for some comparisons, so those are flipped here.
#
# Pathway membership (inputs/kegg_gene_sets.csv) is from the KEGG REST API
# (rest.kegg.jp/link/dme/<pathway>, downloaded 11/20/2026), with KEGG CG IDs
# converted to FlyBase IDs using FlyBase fbgn_annotation_ID (FB2026_03).
# ============================================================

lab <- c(Hypoxia = "8 kPa", Normoxia = "21 kPa", Control = "Control")

# Novogene comparison name -> manuscript comparison and whether to flip it
orient <- function(cmp) {
  p <- strsplit(sub("_Day[26]", "", cmp), "vs")[[1]]
  p <- sub("_Day[26]", "", p)
  flip <- (p[1] == "Control") || (p[1] == "Normoxia" && p[2] == "Hypoxia")
  first <- if (flip) p[2] else p[1]; second <- if (flip) p[1] else p[2]
  list(label = paste(lab[[first]], "vs", lab[[second]]), flip = flip)
}
day_of <- function(f) if (grepl("Day2", f)) "Day 2" else if (grepl("Day6", f)) "Day 6" else "Pooled"

# ---------------------------------------------------------------- Table S10: KEGG enrichment
files <- c(Sys.glob("inputs/*KEGGenrich_significant.csv"), Sys.glob("inputs/KEGG-day2,6/*.csv"))
s10a <- do.call(rbind, lapply(files, function(f) {
  x <- read.csv(f, check.names = FALSE)
  if (!nrow(x)) return(NULL)
  b <- basename(f); cmp <- sub("\\..*", "", b); o <- orient(cmp)
  dir <- sub(".*\\.(up|down|all)_KEGG.*", "\\1", b)
  if (o$flip && dir != "all") dir <- c(up = "down", down = "up")[[dir]]
  data.frame(Comparison = o$label, Day = day_of(b),
             Direction = c(up = "Higher", down = "Lower", all = "All DEGs")[[dir]],
             KEGG_ID = x$KEGGID, Pathway = x$Description, Genes = x$Count,
             GeneRatio = x$GeneRatio, padj = signif(x$padj, 3))
}))
s10a <- s10a[order(s10a$Day, s10a$Comparison, s10a$Direction, s10a$padj), ]
write.csv(s10a, "outputs/table_S10_KEGG_by_day.csv", row.names = FALSE)

# ---------------------------------------------------------------- Table S11: pathway fold changes
sets <- read.csv("inputs/kegg_gene_sets.csv")
degs <- c("HypoxiavsControl", "NormoxiavsControl", "NormoxiavsHypoxia",
          "Control_Day2vsHypoxia_Day2", "Control_Day2vsNormoxia_Day2", "Hypoxia_Day2vsNormoxia_Day2",
          "Control_Day6vsHypoxia_Day6", "Control_Day6vsNormoxia_Day6", "Hypoxia_Day6vsNormoxia_Day6")
s10b <- do.call(rbind, lapply(degs, function(cmp) {
  d <- read.csv(paste0("inputs/", cmp, "_deg.csv"), check.names = FALSE)
  o <- orient(cmp); fc <- if (o$flip) -d$log2FoldChange else d$log2FoldChange
  do.call(rbind, lapply(split(sets, sets$pathway), function(s) {
    i <- d$gene_id %in% s$flybase_id
    q <- quantile(fc[i], c(0.25, 0.5, 0.75), na.rm = TRUE)
    data.frame(Comparison = o$label, Day = day_of(cmp), Pathway = s$pathway[1], KEGG_ID = s$kegg_id[1],
               Genes_tested = sum(i), Median_log2FC = round(q[[2]], 3),
               IQR = paste0(round(q[[1]], 2), " to ", round(q[[3]], 2)),
               Higher_padj05 = sum(i & d$padj < 0.05 & fc > 0, na.rm = TRUE),
               Lower_padj05 = sum(i & d$padj < 0.05 & fc < 0, na.rm = TRUE))
  }))
}))
s10b$Day <- factor(s10b$Day, c("Pooled", "Day 2", "Day 6"))
s10b <- s10b[order(s10b$Pathway, s10b$Day, s10b$Comparison), ]
write.csv(s10b, "outputs/table_S11_pathway_fold_changes.csv", row.names = FALSE)
print(s10b, row.names = FALSE)
