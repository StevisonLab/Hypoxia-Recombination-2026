# Archived: simple ANOVAs

These files are not used in the manuscript and are kept for the record only.

`simple anovas for hypoxia.R` ran exploratory linear models:
- X-chromosome map length for per-vial, cross-pooled and brood-pooled estimates;
- fecundity, as total eggs per vial;
- egg-to-adult viability, as untransformed proportions without the egg-density covariate.

The analyses reported in the paper (Tables 2–4, Figure 3, Tables S3–S7) come from `code/9_hypoxia_figure_3_script.R`. Its fecundity model uses eggs per female per day, and its viability model uses angular-transformed viability with egg density as a covariate, as described in the Methods. The pooled map-length analyses are not reported.

Files: the script, the hand-made inputs only it used (`fecundity_anova_input.csv`, `viability_anova_input.csv`, `pooled_by_brood_anova_input.csv`), and the outputs only it wrote. `inputs/pooled_by_cross_anova_input.csv` and `outputs/table_S10`–`S12` remain in place because `code/8_hypoxia_pooling_analysis.R` also uses those file names.
