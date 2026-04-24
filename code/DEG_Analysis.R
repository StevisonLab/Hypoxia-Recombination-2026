setwd("~/GitHub/Hypoxia-Recombination-2026/")

# PART 1: LOAD YOUR DATA
library(UpSetR)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(circlize)
library(RColorBrewer)

Hypoxia_vs_Control <- read.csv("inputs/HypoxiavsControl_deg.csv")
Normoxia_vs_Control <- read.csv("inputs/NormoxiavsControl_deg.csv")
Normoxia_vs_Hypoxia <- read.csv("inputs/NormoxiavsHypoxia_deg.csv")

Normoxia_vs_Control$Comparison <- "Chamber 21% vs Control 21%"
Normoxia_vs_Hypoxia$Comparison <- "Chamber 21% vs Chamber 8%"  
Hypoxia_vs_Control$Comparison <- "Chamber 8% vs Control 21%"
# PART 2: FILTER SIGNIFICANT DEGs

p_threshold <- 0.05

filter_deg <- function(df) {
  df %>%
    mutate(
      log2FC = log2FoldChange,
      neglog10p = -log10(padj)
    ) %>%
    filter(
      !is.na(log2FC),
      !is.na(padj),
      !is.infinite(neglog10p),
      padj < p_threshold
    )
}
Hypoxia_vs_Control_sig  <- filter_deg(Hypoxia_vs_Control)
Normoxia_vs_Control_sig <- filter_deg(Normoxia_vs_Control)
Normoxia_vs_Hypoxia_sig <- filter_deg(Normoxia_vs_Hypoxia)

all_genes <- bind_rows(
  Hypoxia_vs_Control_sig,
  Normoxia_vs_Control_sig,
  Normoxia_vs_Hypoxia_sig
  
)

cat("\n=== DATA LOADING SUMMARY ===\n")
cat("Control vs Normoxia: ", nrow(Normoxia_vs_Control_sig), " significant DEGs\n")
cat("Normoxia vs Hypoxia: ", nrow(Normoxia_vs_Hypoxia_sig), " significant DEGs\n")
cat("Control vs Hypoxia: ", nrow(Hypoxia_vs_Control_sig), " significant DEGs\n")
cat("Total unique genes: ",
    length(unique(all_genes$gene_id)), "\n\n")

# PART 3: UPSET PLOT 
#Figure 1a
# Ensure your list is ready
deg_list <- list( 
  "Chamber 8% vs Control 21%"  = Hypoxia_vs_Control_sig$gene_id,
  "Chamber 21% vs Control 21%" = Normoxia_vs_Control_sig$gene_id,
  "Chamber 21% vs Chamber 8%"  = Normoxia_vs_Hypoxia_sig$gene_id
)

png("outputs/1_upset_plot_DEGs_padj.png", width = 1100, height = 900, res = 150)
upset(
  fromList(deg_list),
  
  # Force all 3 sets to show in this exact order
  sets = c(
    "Chamber 21% vs Chamber 8%",
    "Chamber 21% vs Control 21%",
    "Chamber 8% vs Control 21%"
  ),
  keep.order = TRUE,
  
  # This shows empty intersections (needed to display the empty set)
  # empty.intersections = "on",
  
  main.bar.color = "gray38",
  matrix.color = "red",
  sets.bar.color = c("gray30", "gray40", "gray60"),  # 3 colors, matches 3 sets
  
  text.scale = c(1.8, 1.5, 1.5, 1.3, 1.8, 1.5),
  point.size = 5,
  line.size = 2,
  mb.ratio = c(0.6, 0.4),
  mainbar.y.label = "Gene Intersections",
  sets.x.label = "DEGs per Comparison"
)
dev.off()


# EXTRACT & SAVE INTERSECTION GENES
# Define correct sets
set_8_vs_ctrl   <- Hypoxia_vs_Control_sig$gene_id
set_21_vs_ctrl  <- Normoxia_vs_Control_sig$gene_id
set_21_vs_8     <- Normoxia_vs_Hypoxia_sig$gene_id

# Intersection: Chamber 8% vs Control 21% AND Chamber 21% vs Control 21%
common_genes <- intersect(set_8_vs_ctrl, set_21_vs_ctrl)

cat("Number of genes in intersection:", length(common_genes), "\n")

# Extract full information
common_genes_info <- all_genes %>%
  filter(gene_id %in% common_genes) %>%
  mutate(
    Direction = ifelse(log2FoldChange > 0, "Up", "Down"),
    gene_label = ifelse(is.na(gene_name) | gene_name == "",
                        gene_id, gene_name)
  ) %>%
  arrange(gene_id, Comparison)

# Save to CSV
write.csv(common_genes_info,  "outputs/common_140 genes.csv",  row.names = FALSE)

# PART 8: DETAILED STATISTICS AND COMMON DEGs

cat("\n=== DETAILED ANALYSIS ===\n\n")

# Create DEG lists per comparison
deg_list <- list(
  Normoxia_vs_Control = Normoxia_vs_Control_sig$gene_id,
  Normoxia_vs_Hypoxia = Normoxia_vs_Hypoxia_sig$gene_id,
  Hypoxia_vs_Control  = Hypoxia_vs_Control_sig$gene_id
)

# Intersection analysis
common_all <- Reduce(intersect, deg_list)
common_Control_norm <- intersect(deg_list$Normoxia_vs_Control,
                                 deg_list$Hypoxia_vs_Control)
common_norm_hyp <- intersect(deg_list$Normoxia_vs_Hypoxia,
                             deg_list$Hypoxia_vs_Control)
common_Control_norm_only <- intersect(deg_list$Normoxia_vs_Control,
                                      deg_list$Normoxia_vs_Hypoxia)

cat("DEG OVERLAP ANALYSIS:\n")
cat("Genes in ALL three comparisons:", length(common_all), "\n")
cat("Genes in Control vs Normoxia ONLY:",
    length(setdiff(deg_list$Normoxia_vs_Control,
                   c(deg_list$Normoxia_vs_Hypoxia,
                     deg_list$Hypoxia_vs_Control))), "\n")
cat("Genes in Normoxia vs Hypoxia ONLY:",
    length(setdiff(deg_list$Normoxia_vs_Hypoxia,
                   c(deg_list$Normoxia_vs_Control,
                     deg_list$Hypoxia_vs_Control))), "\n")
cat("Genes in Control vs Hypoxia ONLY:",
    length(setdiff(deg_list$Hypoxia_vs_Control,
                   c(deg_list$Normoxia_vs_Control,
                     deg_list$Normoxia_vs_Hypoxia))), "\n\n")

# Show common DEGs
if (length(common_all) > 0) {
  
  cat("DEGs COMMON TO ALL COMPARISONS:\n")
  cat("─────────────────────────────────────────\n")
  
  common_info <- all_genes %>%
    filter(gene_id %in% common_all) %>%
    mutate(
      Direction = ifelse(log2FoldChange > 0, "Up", "Down"),
      gene_label = ifelse(is.na(gene_name) | gene_name == "",
                          gene_id, gene_name)
    ) %>%
    select(gene_id, gene_label, Comparison,
           log2FoldChange, padj, Direction) %>%
    arrange(gene_id, Comparison)
  
  print(common_info, n = Inf)
  cat("\n")
  
  # Save common DEGs
  write.csv(common_info,
            "common_DEGs_all_comparisons_padj.csv",
            row.names = FALSE)
  
  cat("✓ Common DEGs saved: common_DEGs_all_comparisons.csv\n\n")
}

# Top DEGs by significance in each comparison
cat("TOP 5 MOST SIGNIFICANT DEGs PER COMPARISON:\n")
cat("─────────────────────────────────────────\n\n")

for (comp in c("Normoxia_vs_Control",
               "Normoxia_vs_Hypoxia",
               "Hypoxia_vs_Control")) {
  
  cat(comp, ":\n")
  
  top_genes <- all_genes %>%
    filter(Comparison == comp) %>%
    arrange(padj) %>%
    head(5) %>%
    mutate(
      Direction = ifelse(log2FoldChange > 0, "Up", "Down"),
      gene_label = ifelse(is.na(gene_name) | gene_name == "",
                          gene_id, gene_name)
    ) %>%
    select(gene_id, gene_label,
           log2FoldChange, padj, Direction)
  
  print(top_genes)
  cat("\n")
}

# PART 9: EXPORT SUMMARY TABLE

summary_table <- all_genes %>%
  mutate(
    Direction = ifelse(log2FoldChange > 0, "Up", "Down"),
    gene_label = ifelse(is.na(gene_name) | gene_name == "",
                        gene_id, gene_name)
  ) %>%
  select(
    gene_id,
    gene_label,
    Comparison,
    log2FoldChange,
    padj,
    Direction
  ) %>%
  arrange(Comparison, padj)

write.csv(summary_table,  "outputs/summary_all_significant_DEGs_padj.csv", row.names = FALSE)

cat("✓ Summary table saved: summary_all_significant_DEGs.csv\n")

# FINAL SUMMARY

cat("Total unique DEGs analyzed:",
    length(unique(all_genes$gene_id)), "\n")
cat("Total DEG-comparison pairs:",
    nrow(all_genes), "\n")

print("✓ DEG overlap and summary analysis completed successfully!")
