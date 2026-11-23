# code

R scripts for the hypoxia experiment. Unless noted, scripts set the working directory to the repository root (`setwd("~/GitHub/Hypoxia-Recombination-2026/")`); edit that line to match where you cloned the repository. Inputs are read from `inputs/` and results are written to `outputs/`.

## Figure and table map

| Paper item | Script | Output(s) |
|---|---|---|
| Figure 1, Figure S1 | — (diagrams) | — |
| Figure 2a (pooled UpSet) | `DEG_Analysis.R` | `1_upset_plot_DEGs_padj.png` |
| Figure 2b (KEGG) | `KEGG_Analysis.R` | `dotplot_COMBINED_KEGG_allcomparisons.png` |
| Figure 2 volcano panels | `DEG_Volcano_Analysis.R` | `volcano_*_padj.png` |
| Table 1 | `full_dataset_analysis.R` | `full_dataset_output.csv` |
| Figure 3; Tables 2–4; Tables S3–S7 | `8b_hypoxia_per_vial_fits.R` → `8c_hypoxia_build_anova_inputs.R` → `9_hypoxia_figure_3_script.R` (after script 2) | `figure_three_*` tables, `figure_three_c_boxplot.pdf/.png` |
| Figure S2 (a, standard curves; b, viability vs. egg density) | `standard_curve_script.R` (after script 9) | `figure_S2_standard_curves.pdf/.png` |
| Figure S3 | `FigureS3_volcano_upset.R` | `figure_S3_volcano_upset.pdf/.png` |
| Figures S4, S5 | `Go_Analysis.R` | `dotplot_UP_regulated_GO.png`, `dotplot_DOWN_regulated_GO.png` |
| Figure S6 | `figure_S6_HIF_genes.R` | `figure_S6_HIF_genes.pdf/.png` |
| Table S8 | `DEG_Analysis.R`, `DEG_Daysummary_Analysis.R` | `summary_all_significant_DEGs_padj.csv`, `deg_up_down_day_summary.csv` |
| Table S9 | `table_S9_hypoxia_genes.R` | `table_S9_hypoxia_genes.csv` |
| Tables S10, S11 | `table_S10_S11_KEGG_by_day.R` | `table_S10_KEGG_by_day.csv`, `table_S11_pathway_fold_changes.csv` |

Figure panels from several scripts (Figures 2 and 3) were assembled outside R.

## Fecundity, survival and recombination

To reproduce the recombination results, run these in order:

1. `2_hypoxia_dataset_compilation_script.R` builds the per-vial class counts (`outputs/13_multi_locus_individual_vials_dataset.csv`).
2. `8b_hypoxia_per_vial_fits.R` fits the counting model to each vial. This is the slow step (about 3 minutes per vial), so it can be split across workers:

    ```
    cd inputs
    Rscript ../code/8b_hypoxia_per_vial_fits.R 1 4   # worker 1 of 4; start workers 2-4 the same way
    ```

    Vials that already have a fit are skipped, so an interrupted run can be restarted.
3. `8c_hypoxia_build_anova_inputs.R` combines the fits and writes the two `inputs/21_output_for_anova*.csv` files, matching each fit to its vial by vial ID.
4. `9_hypoxia_figure_3_script.R`, in the same R session as script 2.

The combined per-vial fits are included (`outputs/21_per_vial_H3_mle_output.csv`), so steps 1–2 can be skipped to rerun only the statistics.

| Script | What it does |
|---|---|
| `1_hypoxia_driver_script.R` | Intended to run scripts 2–8 in order. Not yet updated for this repository. |
| `2_hypoxia_dataset_compilation_script.R` | Compiles the raw egg and adult counts. Applies the digital egg-count correction, using the Type II regression coefficients from `standard_curve_script.R`. Builds the marker-free and multiply-marked viability datasets and the single-locus and per-vial multi-locus datasets (`outputs/D_`–`F_`, `4_`–`13b_`). The cross- and brood-pooling section at the end is unfinished. |
| `3_hypoxia_single_locus_custom_functions.R` | Functions for the single-locus marker-recovery tests. |
| `4_hypoxia_multi_locus_custom_functions.R` | Likelihood functions for the extended *Cx(Co)<sup>m</sup>* counting model (Koury et al. 2026). |
| `5_hypoxia_organismal_analysis.R` | Marker-free viability regression. |
| `6_hypoxia_single_locus_analysis.R` | Single-locus goodness-of-fit (G) tests by data collector (Tables S1–S2). |
| `7_hypoxia_multi_locus_analysis.R` | Counting-model fits (H0–H2) to the full experiment. |
| `8_hypoxia_pooling_analysis.R` | Counting-model (H3) fits to treatment-, cross- and brood-pooled data, with ANOVAs of the pooled estimates (not reported in the paper). Requires `dfoptim`. Each fit takes about 2 minutes per experimental unit. |
| `full_dataset_analysis.R` | Counting-model (H3) fit to the treatment-pooled data (`inputs/full_dataset.csv`): map length and interference for Table 1. |
| `8b_hypoxia_per_vial_fits.R` | Fits the H3 counting model to each of the 247 vials. Each fit is saved under its vial ID in `outputs/21_per_vial_fits/` (not included in the repository; the combined file is). Run from `inputs/`. Requires `dfoptim`. |
| `8c_hypoxia_build_anova_inputs.R` | Combines the per-vial fits into `outputs/21_per_vial_H3_mle_output.csv`, joins them to the brood, treatment and data-collector labels by vial ID, and writes `inputs/21_output_for_anova.csv` and `inputs/21_output_for_anova_single_interval.csv`. Run from `inputs/`. |
| `9_hypoxia_figure_3_script.R` | Fecundity ANOVA, survival ANCOVA, and per-vial ANOVAs of total map length, the five intervals and crossover interference (all Type III), with estimated marginal means and contrasts. Draws the Figure 3 panels. Run from `inputs/` in the same R session after script 2, which creates the datasets it uses. Requires `car`, `emmeans`, `ggpattern`. |
| `9_hypoxia_tables_and_figures_script.R` | Earlier version of the Figure 3 script (`figure_two_*` outputs). |
| `Melika_RecombinationFrequency.R` | Observed recombination fractions per interval and vial (`outputs/14_recombination_frequency_per_vial.csv`). |
| `standard_curve_script.R` | Figure S2. Type II (major axis) regression of each observer's digital egg counts on microscope counts. Checks that the coefficients match those used in script 2. Requires `lmodel2`. |

## RNA-seq

Differential expression (DESeq2) and GO/KEGG enrichment (clusterProfiler) were performed by Novogene. The scripts below plot and summarize Novogene's results and do not re-run those analyses.

| Script | What it does |
|---|---|
| `DEG_Analysis.R` | Pooled-day DEG counts and overlaps. Draws the Figure 2a UpSet plot. |
| `DEG_Volcano_Analysis.R` | Pooled-day volcano plots (Figure 2). |
| `KEGG_Analysis.R` | KEGG enrichment dot plot (Figure 2b). |
| `Go_Analysis.R` | GO enrichment dot plots for up- and down-regulated genes (Figures S4, S5). |
| `DEG_Daysummary_Analysis.R` | Up/down DEG counts per comparison and day (Table S8). |
| `FigureS3_volcano_upset.R` | Figure S3: per-day volcano plots (a–f) and per-day UpSet plot (g). |
| `DEG_Volcano_perday_Analysis.R`, `DEG_DayUpset_Analysis.R` | Earlier per-day volcano and UpSet plots, superseded by `FigureS3_volcano_upset.R`. |
| `figure_S6_HIF_genes.R` | Figure S6: FPKM of the core HIF pathway genes (Novomagic export), annotated with the genome-wide DESeq2 results. |
| `table_S9_hypoxia_genes.R` | Table S9: DESeq2 results for known hypoxia-responsive genes and for the homologous recombination genes discussed in the text, pooled and by day, plus mean FPKM across all 24 samples. |
| `table_S10_S11_KEGG_by_day.R` | Table S10: KEGG enrichment for every comparison, pooled and by day. Table S11: median log2 fold changes across all genes in the ribosome, oxidative phosphorylation, glycolysis and TCA cycle pathways (membership in `inputs/kegg_gene_sets.csv`). |

## R packages

`dplyr`, `tidyr`, `ggplot2`, `patchwork`, `ggrepel`, `cowplot`, `UpSetR`, `ComplexHeatmap`, `circlize`, `RColorBrewer`, `ggnewscale`, `gtable`, `readr`, `car`, `emmeans`, `ggpattern`, `lmodel2`, `dfoptim`.
