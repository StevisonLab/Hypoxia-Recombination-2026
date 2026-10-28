# ==========================
# UP/DOWN REGULATION SUMMARY
# ==========================
library(dplyr)

# Thresholds — must match volcano plots and UpSet plot
padj_cutoff      <- 0.05
log2fc_threshold <- 1

# --------------------------
# Import CSV files
# --------------------------
# Day 2
Hypoxia_vs_Normoxia_Day2_deg <- read.csv("inputs/Hypoxia_Day2vsNormoxia_Day2_deg.csv", row.names = 1)
Control_vs_Normoxia_Day2_deg <- read.csv("inputs/Control_Day2vsNormoxia_Day2_deg.csv", row.names = 1)
Control_vs_Hypoxia_Day2_deg  <- read.csv("inputs/Control_Day2vsHypoxia_Day2_deg.csv",  row.names = 1)

# Day 6
Hypoxia_vs_Normoxia_Day6_deg <- read.csv("inputs/Hypoxia_Day6vsNormoxia_Day6_deg.csv", row.names = 1)
Control_vs_Normoxia_Day6_deg <- read.csv("inputs/Control_Day6vsNormoxia_Day6_deg.csv", row.names = 1)
Control_vs_Hypoxia_Day6_deg  <- read.csv("inputs/Control_Day6vsHypoxia_Day6_deg.csv",  row.names = 1)

# --------------------------
# Function to count up/down genes
# --------------------------
count_deg_direction <- function(df, comparison_name,
                                padj_cutoff = 0.05,
                                log2fc_threshold = 1) {   # <-- NEW parameter
  df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      direction = case_when(
        padj < padj_cutoff & log2FoldChange >  log2fc_threshold ~ "Upregulated",
        padj < padj_cutoff & log2FoldChange < -log2fc_threshold ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    ) %>%
    count(direction) %>%
    mutate(comparison = comparison_name)
}

# --------------------------
# APPLY TO ALL COMPARISONS
# --------------------------
results <- bind_rows(
  
  # Day 2
  count_deg_direction(Control_vs_Hypoxia_Day2_deg,  "Day2: Control vs Hypoxia",  padj_cutoff, log2fc_threshold),
  count_deg_direction(Control_vs_Normoxia_Day2_deg, "Day2: Control vs Normoxia", padj_cutoff, log2fc_threshold),
  count_deg_direction(Hypoxia_vs_Normoxia_Day2_deg, "Day2: Hypoxia vs Normoxia", padj_cutoff, log2fc_threshold),
  
  # Day 6
  count_deg_direction(Control_vs_Hypoxia_Day6_deg,  "Day6: Control vs Hypoxia",  padj_cutoff, log2fc_threshold),
  count_deg_direction(Control_vs_Normoxia_Day6_deg, "Day6: Control vs Normoxia", padj_cutoff, log2fc_threshold),
  count_deg_direction(Hypoxia_vs_Normoxia_Day6_deg, "Day6: Hypoxia vs Normoxia", padj_cutoff, log2fc_threshold)
)

# --------------------------
# CLEAN TABLE FORMAT
# --------------------------
summary_table <- results %>%
  tidyr::pivot_wider(
    names_from  = direction,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    Total_Significant = Upregulated + Downregulated
  ) %>%
  select(comparison, Upregulated, Downregulated, Total_Significant)

# Print with threshold reminder
cat("Thresholds applied: |log2FC| >", log2fc_threshold, "& padj <", padj_cutoff, "\n\n")
print(summary_table)

# --------------------------
# SAVE TABLE
# --------------------------
write.csv(summary_table, "outputs/deg_up_down_day_summary.csv", row.names = FALSE)
cat("✓ Saved: outputs/deg_up_down_day_summary.csv\n")
