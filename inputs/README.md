# inputs

Raw data and analysis input files for the hypoxia experiment. Data-collector initials in file names: CH and EH (digital egg counts), and MS, NRR and SK (adult phenotype scoring). Treatments: Control = outside-chamber *P*<sub>O<sub>2</sub></sub> ≈ 21 kPa; Normoxia = Chamber *P*<sub>O<sub>2</sub></sub> ≈ 21 kPa; Hypoxia = Chamber *P*<sub>O<sub>2</sub></sub> ≈ 8 kPa.

## Raw counts

| File | Contents |
|---|---|
| `A_CH_hypoxia_raw_marker_free_egg_counts.csv`, `B_EH_hypoxia_raw_marker_free_egg_counts.csv` | Digital egg counts for marker-free crosses, by observer. |
| `C_NRR_hypoxia_raw_marker_free_adult_counts.csv` | Adult counts for marker-free crosses. |
| `1_CH_hypoxia_raw_multiply_marked_egg_counts.csv`, `2_EH_hypoxia_raw_multiply_marked_egg_counts.csv` | Digital egg counts for multiply-marked crosses, by observer. |
| `3_MS_…`, `3_NRR_…`, `3_SK_hypoxia_raw_multiply_phenotypic_class_counts.csv` | Adult phenotypic class counts for multiply-marked crosses, by data collector. |
| `3_ALL_hypoxia_raw_multiply_marked_phenotypic_class_counts.csv` | The three collectors' phenotypic class counts combined. |

## Calibration

| File | Contents |
|---|---|
| `S1_hypoxia_standard_curves.csv` | Calibration dataset (Schneider et al. 2025): microscope egg counts (SK) and digital counts (CH, EH) for the same 30 vials. Used for Figure S2 and the egg-count correction in script 2. |

## RNA-seq (from Novogene)

| File | Contents |
|---|---|
| `HypoxiavsControl_deg.csv`, `NormoxiavsControl_deg.csv`, `NormoxiavsHypoxia_deg.csv` | Genome-wide DESeq2 results, days pooled. `log2FoldChange` is log<sub>2</sub>(first / second treatment in the file name). Includes per-sample read counts and FPKM. |
| `Control_Day2vsHypoxia_Day2_deg.csv` and the five other `*_Day2*`/`*_Day6*` files | The same, for each sampling day. |
| `*.up_GOenrich_significant.csv`, `*.down_GOenrich_significant.csv`, `*.all_KEGGenrich_significant.csv` | Significant GO and KEGG terms (clusterProfiler, BH-adjusted *P* < 0.05), days pooled. |
| `GO-day2,6/`, `KEGG-day2,6/` | Significant GO and KEGG terms for each sampling day. |
| `novomagic_KeyGeneID_20260424_121126043_get_annot.txt` | Novomagic gene-subset export (24 April 2026): per-sample FPKM and DESeq2 results for 71 selected genes. Unmodified download (tab-separated). Used for Figure S6 and Table S9. |

## Analysis input files

These files were assembled outside the scripts from earlier outputs. They are being replaced by files written in code, starting with the per-vial recombination estimates.

| File | Used by |
|---|---|
| `21_output_for_anova.csv`, `21_output_for_anova_single_interval.csv` | Per-vial counting-model estimates for the map-length and interval analyses (script 9). |
| `full_dataset.csv` | Treatment-pooled phenotypic class counts (scripts 8 and `full_dataset_analysis.R`). |
| `14_multi_locus_pooled_by_cross_dataset.csv`, `15_multi_locus_pooled_by_brood_dataset.csv` | Cross- and brood-pooled class counts (script 8). |
| `pooled_by_cross_anova_input.csv` | Cross-pooled map-length estimates (script 8). |
| `interference_anova_inputs.csv` | Inputs for the crossover interference analysis. |
| `coincidence by genetic length formula.xlsx` | Coefficient of coincidence calculations. |
