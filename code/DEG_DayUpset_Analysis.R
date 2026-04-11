# ==========================
# UPSET PLOTS PER DAY
# ==========================

library(UpSetR)
library(dplyr)

# Day 2
Hypoxia_vs_Normoxia_Day2_deg <- read.csv("inputs/Hypoxia_Day2vsNormoxia_Day2_deg.csv", row.names = 1)
Control_vs_Normoxia_Day2_deg <- read.csv("inputs/Control_Day2vsNormoxia_Day2_deg.csv", row.names = 1)
Control_vs_Hypoxia_Day2_deg  <- read.csv("inputs/Control_Day2vsHypoxia_Day2_deg.csv",  row.names = 1)

# Day 6
Hypoxia_vs_Normoxia_Day6_deg <- read.csv("inputs/Hypoxia_Day6vsNormoxia_Day6_deg.csv", row.names = 1)
Control_vs_Normoxia_Day6_deg <- read.csv("inputs/Control_Day6vsNormoxia_Day6_deg.csv", row.names = 1)
Control_vs_Hypoxia_Day6_deg  <- read.csv("inputs/Control_Day6vsHypoxia_Day6_deg.csv",  row.names = 1)

p_threshold <- 0.05

# --------------------------
# Function to filter DEGs
# --------------------------
filter_deg <- function(df) {
  df %>%
    filter(!is.na(log2FoldChange),
           !is.na(padj),
           padj < p_threshold)
}

# ==========================
# DAY 2
# ==========================

CH_D2 <- filter_deg(Control_vs_Hypoxia_Day2_deg)
CN_D2 <- filter_deg(Control_vs_Normoxia_Day2_deg)
HN_D2 <- filter_deg(Hypoxia_vs_Normoxia_Day2_deg)

deg_list_D2 <- list(
  "Chamber 8% vs Control 21%"  = rownames(CH_D2),
  "Chamber 21% vs Control 21%" = rownames(CN_D2),
  "Chamber 21% vs Chamber 8%"  = rownames(HN_D2)
)

png("outputs/upset_Day2.png", width = 1100, height = 900, res = 150)

upset(
  fromList(deg_list_D2),
  
  sets = c(
    "Chamber 21% vs Chamber 8%",
    "Chamber 21% vs Control 21%",
    "Chamber 8% vs Control 21%"
  ),
  keep.order = TRUE,
  
  main.bar.color = "gray38",
  matrix.color = "red",
  sets.bar.color = c("gray30", "gray40", "gray60"),
  
  text.scale = c(1.8, 1.5, 1.5, 1.3, 1.8, 1.5),
  point.size = 5,
  line.size = 2,
  mb.ratio = c(0.6, 0.4),
  mainbar.y.label = "Gene Intersections",
  sets.x.label = "DEGs per Comparison"
)

dev.off()

# ==========================
# DAY 6
# ==========================

CH_D6 <- filter_deg(Control_vs_Hypoxia_Day6_deg)
CN_D6 <- filter_deg(Control_vs_Normoxia_Day6_deg)
HN_D6 <- filter_deg(Hypoxia_vs_Normoxia_Day6_deg)

deg_list_D6 <- list(
  "Chamber 8% vs Control 21%"  = rownames(CH_D6),
  "Chamber 21% vs Control 21%" = rownames(CN_D6),
  "Chamber 21% vs Chamber 8%"  = rownames(HN_D6)
)

png("outputs/upset_Day6.png", width = 1100, height = 900, res = 150)

upset(
  fromList(deg_list_D6),
  
  sets = c(
    "Chamber 21% vs Chamber 8%",
    "Chamber 21% vs Control 21%",
    "Chamber 8% vs Control 21%"
  ),
  keep.order = TRUE,
  
  main.bar.color = "gray38",
  matrix.color = "blue",   # optional: different color for day 6
  
  sets.bar.color = c("gray30", "gray40", "gray60"),
  
  text.scale = c(1.8, 1.5, 1.5, 1.3, 1.8, 1.5),
  point.size = 5,
  line.size = 2,
  mb.ratio = c(0.6, 0.4),
  mainbar.y.label = "Gene Intersections",
  sets.x.label = "DEGs per Comparison"
)

dev.off()
