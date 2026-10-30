setwd("~/GitHub/Hypoxia-Recombination-2026/")

# ============================================================
# Figure S6: ovarian expression of core HIF pathway genes
#
# Per-sample FPKM for the four genes comes from a Novomagic gene-subset
# export (24 April 2026). The statistic shown above each panel is the
# pooled Chamber 8 kPa vs Control log2 fold change and Benjamini-Hochberg
# adjusted P from the genome-wide DESeq2 analysis (Novogene). No new
# tests are performed; FPKM is plotted for visualization only.
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

genes <- tribble(
  ~ID,           ~Symbol, ~Role,
  "FBgn0266411", "sima",  "HIF-alpha",
  "FBgn0264075", "tgo",   "HIF-beta / ARNT",
  "FBgn0264785", "fga",   "prolyl hydroxylase, PHD",
  "FBgn0041174", "Vhl",   "von Hippel-Lindau"
)

# Labels in the PO2 notation used in Figures 1-3 and S3-S5 (no "approx.").
# Wrap the labels in list(): as.expression() on a single call would split it.
po2 <- function(kpa) bquote("Chamber" ~ italic(P)[O[2]] ~ .(paste(kpa, "kPa")))
trt_levels <- c("Control", "Chamber21", "Chamber8")
trt_labels <- as.expression(list("Control", po2(21), po2(8)))
trt_colors <- c(Control = "#bf812d", Chamber21 = "#01665e", Chamber8 = "#c51b7d")

novomagic <- read.delim("inputs/novomagic_KeyGeneID_20260424_121126043_get_annot.txt",
                        check.names = FALSE, na.strings = c("--", "NA", "")) %>%
  filter(ID %in% genes$ID)
stopifnot(nrow(novomagic) == 4)

# Pooled Chamber 8 kPa vs Control result from the genome-wide DESeq2 analysis
# (Novogene column "HypoxiavsControl" = log2(Hypoxia / Control))
stats <- novomagic %>%
  transmute(ID,
            log2FC = as.numeric(HypoxiavsControl_log2FoldChange),
            padj   = as.numeric(HypoxiavsControl_padj))

fpkm <- novomagic %>%
  select(ID, ends_with("_fpkm")) %>%
  pivot_longer(-ID, names_to = "sample", values_to = "FPKM") %>%
  extract(sample, c("trt", "Day", "Rep"), "^(Con|Nor|Hyp)_399_(\\d)_(\\d)_fpkm$") %>%
  mutate(Treatment = factor(recode(trt, Con = "Control", Nor = "Chamber21", Hyp = "Chamber8"),
                            levels = trt_levels),
         Day = factor(paste("Day", Day), levels = c("Day 2", "Day 6")),
         FPKM = as.numeric(FPKM))
stopifnot(all(table(fpkm$ID) == 24))

fpkm_summary <- fpkm %>%
  group_by(ID, Treatment, Day) %>%
  summarise(Mean = mean(FPKM), SD = sd(FPKM), .groups = "drop")

fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))

make_panel <- function(id) {
  g  <- filter(genes, ID == id)
  st <- filter(stats, ID == id)
  note <- sprintf("8 kPa vs Control (pooled): log2FC = %+.2f, padj = %s%s",
                  st$log2FC, fmt_p(st$padj), ifelse(st$padj < 0.05, " *", ""))
  pd <- position_dodge(width = 0.55)
  ggplot(filter(fpkm_summary, ID == id), aes(Day, Mean, colour = Treatment)) +
    geom_point(data = filter(fpkm, ID == id), aes(y = FPKM),
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.55, seed = 1),
               size = 1.4, alpha = 0.45, show.legend = FALSE) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.18, linewidth = 0.6,
                  position = pd) +
    geom_point(size = 2.6, position = pd) +
    scale_colour_manual(values = trt_colors, breaks = trt_levels, labels = trt_labels,
                        name = NULL) +
    labs(title = bquote(italic(.(g$Symbol)) ~ "(" * .(g$Role) * ")"),
         subtitle = note, x = NULL, y = "FPKM") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "bold", size = 10.5),
          plot.subtitle = element_text(size = 7.5, colour = "grey25"),
          legend.position = "none")
}

# Show the legend from one panel only; patchwork can't merge legends whose
# labels are plotmath expressions
panels <- lapply(genes$ID, make_panel)
panels[[4]] <- panels[[4]] + theme(legend.position = "bottom")

figure_S6 <- wrap_plots(panels, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  theme = theme(legend.position = "bottom")) &
  theme(plot.tag = element_text(face = "bold"))

ggsave("outputs/figure_S6_HIF_genes.pdf", figure_S6, width = 7, height = 6.4, device = cairo_pdf)
ggsave("outputs/figure_S6_HIF_genes.png", figure_S6, width = 7, height = 6.4, dpi = 300)
