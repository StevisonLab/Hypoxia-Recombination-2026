# Builds the script 9 input files from the per-vial fits, joining on vial ID. Run from "inputs".
# Vial labels (brood, treatment, data collector) come from the existing input file and match the raw data.
if (dir.exists("../outputs/21_per_vial_fits")) {
  fits <- do.call(rbind, lapply(list.files("../outputs/21_per_vial_fits", full.names = TRUE), read.csv, check.names = FALSE))
  write.csv(fits, "../outputs/21_per_vial_H3_mle_output.csv", row.names = FALSE)
}
fits <- read.csv("../outputs/21_per_vial_H3_mle_output.csv", check.names = FALSE)
for_anova <- read.csv("21_output_for_anova.csv", fileEncoding = "UTF-8-BOM")
stopifnot(nrow(fits) == 247, !anyDuplicated(fits$vial_unique_id), setequal(fits$vial_unique_id, for_anova$vial_unique_id))
estimates <- c("MLE_X_length", paste0("MLE_x", 1:5), "MLE_p")
for_anova[, estimates] <- fits[match(for_anova$vial_unique_id, fits$vial_unique_id), estimates]
single_interval <- for_anova
single_interval[, paste0("MLE_x", 1:5)] <- single_interval[, paste0("MLE_x", 1:5)] * 100 # interval lengths in cM
write.csv(for_anova, "21_output_for_anova.csv", row.names = FALSE)
write.csv(single_interval, "21_output_for_anova_single_interval.csv", row.names = FALSE)
