# outputs

Results written by the scripts in `code/` (see `code/README.md` for which script makes what).

## Figures

| File | Paper item |
|---|---|
| `1_upset_plot_DEGs_padj.png` | Figure 2a |
| `dotplot_COMBINED_KEGG_allcomparisons.png` | Figure 2b |
| `volcano_*_padj.png`, `volcano_combined_*.png` | Figure 2 volcano panels |
| `figure_S2_standard_curves.pdf/.png` | Figure S2 |
| `figure_S3_volcano_upset.pdf/.png` | Figure S3 |
| `dotplot_UP_regulated_GO.png`, `dotplot_DOWN_regulated_GO.png` | Figures S4, S5 |
| `figure_S6_HIF_genes.pdf/.png` | Figure S6 |
| `upset_Day2.png`, `upset_Day6.png`, `upset_plot_combined_D2_D6.png` | Earlier per-day UpSet plots (superseded by Figure S3g) |

## Tables

| File(s) | Paper item |
|---|---|
| `figure_three_a_fecundity_*` | Fecundity: Table 2, Figure 3a, Tables S3–S4 |
| `figure_three_b_*_survival_*` | Egg-to-adult survival: Table 3, Figure 3b, Tables S5–S6 |

Tables S3–S6 are the within-design rows of the `*_contrasts_table.csv` files. Experiment NRR_2 is the marker-free cross and NRR_1 the multiply-marked cross.
| `figure_three_c_length_*`, `figure_three_c_x1`–`x5_*` | X-chromosome map length and the five intervals: Table 4, Figure 3c, Table S7 |
| `summary_all_significant_DEGs_padj.csv`, `deg_up_down_day_summary.csv` | Table S8 (pooled and per-day DEG counts) |
| `table_S9_hypoxia_genes.csv` | Table S9 |
| `common_140 genes.csv`, `common_genes_8pct_vs_21pct_vs_ctrl.csv` | Genes shared between comparisons (from `DEG_Analysis.R`) |

> The map-length, interval and interference results (`figure_three_c_*`) use the re-fitted per-vial estimates and Type III tests.

## Intermediate datasets (script 2)

| File(s) | Contents |
|---|---|
| `D_compiled_raw_marker_free_dataset.csv`, `E_derived_…`, `F_marker_free_viability_dataset.csv` | Marker-free crosses: compiled counts, corrected egg counts, viability |
| `4_compiled_raw_multiply_marked_dataset.csv`, `5_derived_…`, `6_multiply_marked_viability_dataset.csv` | Multiply-marked crosses: the same steps |
| `7_`–`12_*_single_locus_dataset.csv` | Single-locus datasets for each marker (*sc*, *cv*, *v*, *f*, *car*, *y*<sup>+</sup>), by data collector and combined (ALL) |
| `13_multi_locus_individual_vials_dataset.csv`, `13b_cleaned_…`, `13c_…` | Per-vial phenotypic class counts for the counting model |
| `14_`–`16_*` | Cross-, brood- and experiment-pooled class counts |
| `14_recombination_frequency_per_vial.csv` | Observed recombination fractions (`Melika_RecombinationFrequency.R`) |

## Counting-model fits (scripts 7 and 8)

| File(s) | Contents |
|---|---|
| `21_per_vial_H3_mle_output.csv` | Per-vial H3 fits for all 247 vials, keyed by vial ID |
| `21_multi_locus_individual_vials_H3_mle_output.csv`, `full_dataset_output.csv` | Treatment-pooled H3 fits |
| `22_multi_locus_cross_pooled_H3_mle_output.csv`, `23_multi_locus_brood_pooled_H3_mle_output.csv` | Cross- and brood-pooled H3 fits |
| `table_S1_marker_free_viability_regression.csv`, `table_S10_`–`table_S12_*_anova.csv` | Viability regression and pooled-data ANOVAs (file numbering predates the current supplement; not reported in the paper) |

## Other files

- `figure_two_*`: outputs of the earlier Figure 3 script (`9_hypoxia_tables_and_figures_script.R`).
- `coefficient_of_conincidence_dataset.csv`: coefficient of coincidence results.
- `*.xlsx`: spreadsheet working files (post hoc tests, survival ANOVA, full-dataset exploration, coefficient of coincidence, combined per-vial fits).
