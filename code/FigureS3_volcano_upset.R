# ============================================================
# Figure S3: per-day differential expression
#   (a-c) Day 2 volcano plots; (d-f) Day 6 volcano plots
#   (g)   overlap of DEGs among the three comparisons, Day 2 and
#         Day 6 side by side (complements the pooled UpSet, Fig. 2a)
#
# Builds the whole figure in R with panel letters, shared axes and
# one legend. Threshold matches the Methods: DESeq2 BH-adjusted
# P < 0.05, no fold-change cutoff. Run from the repository root.
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)   # >= 1.2 for free()

setwd("~/GitHub/Hypoxia-Recombination-2026/")
p_threshold <- 0.05

# ------------------------------------------------------------
# Treatments and labels (PO2 notation as in Figures 1-2, no "approx.")
# Normoxia = Chamber PO2 21 kPa; Hypoxia = Chamber PO2 8 kPa;
# Control = outside-chamber control
# ------------------------------------------------------------
po2 <- function(kpa) bquote("Chamber" ~ italic(P)[O[2]] ~ .(paste(kpa, "kPa")))

treatment_colors <- c("Chamber 8 kPa"  = "#c51b7d",
                      "Chamber 21 kPa" = "#01665e",
                      "Control"        = "#bf812d",
                      "Not significant" = "gray75")
# Wrap each label in list(): as.expression() on a single call would split it
treatment_labels <- as.expression(list(bquote("Higher in" ~ .(po2(8))),
                                       bquote("Higher in" ~ .(po2(21))),
                                       "Higher in Control",
                                       "Not significant"))

# Each DEG file reports log2FoldChange = log2(first/second treatment
# in the file name). `flip` puts the treatment named in `right` on
# the positive side.
volcano_sets <- tribble(
  ~panel, ~Day,    ~file,                                  ~right,           ~left,            ~flip,
  "a",    "Day 2", "Control_Day2vsHypoxia_Day2_deg.csv",   "Chamber 8 kPa",  "Control",        TRUE,
  "b",    "Day 2", "Control_Day2vsNormoxia_Day2_deg.csv",  "Chamber 21 kPa", "Control",        TRUE,
  "c",    "Day 2", "Hypoxia_Day2vsNormoxia_Day2_deg.csv",  "Chamber 21 kPa", "Chamber 8 kPa",  TRUE,
  "d",    "Day 6", "Control_Day6vsHypoxia_Day6_deg.csv",   "Chamber 8 kPa",  "Control",        TRUE,
  "e",    "Day 6", "Control_Day6vsNormoxia_Day6_deg.csv",  "Chamber 21 kPa", "Control",        TRUE,
  "f",    "Day 6", "Hypoxia_Day6vsNormoxia_Day6_deg.csv",  "Chamber 21 kPa", "Chamber 8 kPa",  TRUE
)

# Column titles (one per comparison, shown on the Day 2 row)
comparison_titles <- list(
  a = bquote(atop(.(po2(8)),  "vs. Control")),
  b = bquote(atop(.(po2(21)), "vs. Control")),
  c = bquote(atop(.(po2(21)), "vs." ~ .(po2(8))))
)

read_deg <- function(file) {
  read.csv(file.path("inputs", file), check.names = FALSE,
           fileEncoding = "UTF-8-BOM")
}

theme_fig <- theme_classic(base_size = 10) +
  theme(axis.text  = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5),
        legend.title = element_text(face = "bold"),
        plot.tag   = element_text(size = 14, face = "bold"))

# ============================================================
# (a-f) Volcano plots
# ============================================================
volcano_data <- volcano_sets %>%
  rowwise() %>%
  reframe(read_deg(file) %>%
            select(gene_name, log2FoldChange, padj) %>%
            mutate(panel = panel, Day = Day,
                   log2FC = if (flip) -log2FoldChange else log2FoldChange,
                   right = right, left = left)) %>%
  filter(!is.na(padj), !is.na(log2FC)) %>%
  mutate(neglog10p = -log10(padj),
         Status = case_when(padj < p_threshold & log2FC > 0 ~ right,
                            padj < p_threshold & log2FC < 0 ~ left,
                            TRUE ~ "Not significant"),
         Status = factor(Status, levels = names(treatment_colors)))

# Shared, symmetric limits across all six panels
x_max <- ceiling(max(abs(volcano_data$log2FC)))
y_max <- max(volcano_data$neglog10p[is.finite(volcano_data$neglog10p)]) * 1.05

# Label the 5 most significant genes on each side of each panel
volcano_labels <- volcano_data %>%
  filter(Status != "Not significant") %>%
  group_by(panel, Status) %>%
  slice_min(padj, n = 5, with_ties = FALSE) %>%
  ungroup()

make_volcano <- function(p) {
  d  <- filter(volcano_data, panel == p)
  lb <- filter(volcano_labels, panel == p)
  ggplot(d, aes(x = log2FC, y = neglog10p, color = Status)) +
    geom_point(data = filter(d, Status == "Not significant"),
               size = 0.7, alpha = 0.5, show.legend = TRUE) +
    geom_point(data = filter(d, Status != "Not significant"),
               size = 1.1, alpha = 0.8, show.legend = TRUE) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed",
               color = "gray40", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray40", linewidth = 0.3) +
    geom_text_repel(data = lb, aes(label = gene_name), size = 2.2,
                    fontface = "italic", color = "black",
                    segment.color = "gray50", segment.size = 0.25,
                    min.segment.length = 0, max.overlaps = Inf, seed = 1,
                    ylim = c(3, NA), box.padding = 0.35) +
    scale_color_manual(values = treatment_colors, breaks = names(treatment_colors),
                       limits = names(treatment_colors), labels = treatment_labels,
                       drop = FALSE) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE,
                                override.aes = list(size = 2.5, alpha = 1))) +
    coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, y_max)) +
    labs(x = expression(log[2]~"fold change"),
         y = expression(-log[10]~"adjusted"~italic(P)),
         color = NULL) +
    theme_fig
}

volcanoes <- lapply(volcano_sets$panel, make_volcano)
names(volcanoes) <- volcano_sets$panel

# Comparison titles on the Day 2 row; day shown in the left y-axis title;
# x-axis titles only on the Day 6 row
for (p in c("a", "b", "c")) {
  volcanoes[[p]] <- volcanoes[[p]] + ggtitle(comparison_titles[[p]]) +
    theme(axis.title.x = element_blank())
}
day_ylab <- function(day) bquote(atop(bold(.(day)), -log[10] ~ "adjusted" ~ italic(P)))
volcanoes$a <- volcanoes$a + labs(y = day_ylab("Day 2"))
volcanoes$d <- volcanoes$d + labs(y = day_ylab("Day 6"))
for (p in c("b", "c", "e", "f")) {
  volcanoes[[p]] <- volcanoes[[p]] + theme(axis.title.y = element_blank())
}

# ============================================================
# (g) UpSet plot, Day 2 and Day 6 side by side
# ============================================================
# Intersections are computed within each day, so each bar is the
# number of genes DE in exactly that combination of comparisons on
# that day. Set sizes match table S8; bars on each day sum to the
# union in table S8 (1,735 on Day 2; 1,222 on Day 6).
upset_sets <- tribble(
  ~Comparison,                         ~Day,    ~file,
  "Chamber 8 kPa vs Control",          "Day 2", "Control_Day2vsHypoxia_Day2_deg.csv",
  "Chamber 21 kPa vs Control",         "Day 2", "Control_Day2vsNormoxia_Day2_deg.csv",
  "Chamber 8 kPa vs Chamber 21 kPa",   "Day 2", "Hypoxia_Day2vsNormoxia_Day2_deg.csv",
  "Chamber 8 kPa vs Control",          "Day 6", "Control_Day6vsHypoxia_Day6_deg.csv",
  "Chamber 21 kPa vs Control",         "Day 6", "Control_Day6vsNormoxia_Day6_deg.csv",
  "Chamber 8 kPa vs Chamber 21 kPa",   "Day 6", "Hypoxia_Day6vsNormoxia_Day6_deg.csv"
)
set_levels <- c("Chamber 21 kPa vs Control",        # bottom row
                "Chamber 8 kPa vs Control",
                "Chamber 8 kPa vs Chamber 21 kPa")  # top row
set_labels <- as.expression(list(               # same order as set_levels
  bquote(.(po2(21)) ~ "vs. Control"),
  bquote(.(po2(8))  ~ "vs. Control"),
  bquote(.(po2(8))  ~ "vs." ~ .(po2(21)))
))

long <- upset_sets %>%
  rowwise() %>%
  reframe(gene_id = read_deg(file) %>%
            filter(!is.na(padj), padj < p_threshold) %>%
            pull(gene_id) %>% unique(),
          Comparison = Comparison, Day = Day)

bars <- long %>%
  group_by(Day, gene_id) %>%
  summarise(Intersection = paste(sort(unique(Comparison)), collapse = " & "),
            .groups = "drop") %>%
  count(Day, Intersection) %>%
  complete(Day, Intersection, fill = list(n = 0))
intersection_order <- bars %>% group_by(Intersection) %>%
  summarise(total = sum(n)) %>% arrange(desc(total)) %>% pull(Intersection)
bars <- mutate(bars, Intersection = factor(Intersection, levels = intersection_order))

# Checks against table S8
print(count(long, Day, Comparison))
print(bars %>% group_by(Day) %>% summarise(union = sum(n)))

day_fills <- c("Day 2" = "#9ecae1", "Day 6" = "#3182bd")
bar_dodge <- position_dodge(width = 0.8)

p_bars <- ggplot(bars, aes(x = Intersection, y = n, fill = Day)) +
  geom_col(width = 0.75, position = bar_dodge, color = "black", linewidth = 0.2) +
  geom_text(aes(label = n), position = bar_dodge, vjust = -0.4, size = 2.3) +
  scale_fill_manual(values = day_fills) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Intersection size", fill = NULL) +
  theme_fig +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.92, 0.8))

dots <- expand_grid(Intersection = factor(intersection_order, levels = intersection_order),
                    Comparison = factor(set_levels, levels = set_levels)) %>%
  mutate(Member = mapply(function(i, c) c %in% strsplit(i, " & ")[[1]],
                         as.character(Intersection), as.character(Comparison)))

p_dots <- ggplot(dots, aes(x = Intersection, y = Comparison)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = c(0.5, 2.5), ymax = c(1.5, 3.5),
           fill = "gray95") +
  geom_line(data = filter(dots, Member), aes(group = Intersection), linewidth = 0.7) +
  geom_point(aes(color = Member), size = 2.5) +
  scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "gray80"), guide = "none") +
  scale_y_discrete(labels = set_labels) +
  labs(x = NULL, y = NULL) +
  theme_fig +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank())

sets <- count(long, Day, Comparison) %>%
  mutate(Comparison = factor(Comparison, levels = set_levels))
p_sets <- ggplot(sets, aes(y = Comparison, x = n, fill = Day)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.75),
           color = "black", linewidth = 0.2, show.legend = FALSE) +
  scale_x_reverse(expand = expansion(mult = c(0.05, 0))) +
  scale_fill_manual(values = day_fills) +
  labs(x = "Set size", y = NULL) +
  theme_fig +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# free() keeps the long set labels from pushing the bar axis title away
tag_g <- ggplot() + theme_void() + labs(tag = "g") +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0, 1))
upset_panel <- tag_g + free(p_bars, type = "label", side = "l") +
  p_sets + p_dots +
  plot_layout(widths = c(1, 3), heights = c(3, 1))

# ============================================================
# Assemble
# ============================================================
for (p in c("a", "b", "c", "d", "e")) {
  volcanoes[[p]] <- volcanoes[[p]] + theme(legend.position = "none")
}
volcanoes$f <- volcanoes$f + theme(legend.position = "bottom")
# Set the legend position on the layout only. Applying it with `&` would
# switch the legend back on in every panel, and patchwork can't merge
# legends whose labels are plotmath expressions.
volcano_grid <- wrap_plots(volcanoes, ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a",
                  theme = theme(legend.position = "bottom"))

figS3 <- wrap_elements(full = volcano_grid) / wrap_elements(full = upset_panel) +
  plot_layout(heights = c(2, 1.15))

ggsave("outputs/figure_S3_volcano_upset.pdf", figS3, width = 8.5, height = 10,
       device = cairo_pdf)
ggsave("outputs/figure_S3_volcano_upset.png", figS3, width = 8.5, height = 10,
       dpi = 300)
