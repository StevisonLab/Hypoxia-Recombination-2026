# ============================================================
# Recombination Frequency Calculation (Pooled & Total RF)
# ============================================================

library(dplyr)
library(tidyr)

# Read file
multi_locus_data <- read.csv("~/GitHub/Hypoxia-Recombination-2026/outputs/13_multi_locus_individual_vials_dataset.csv",
                             header = TRUE)

# ============================================================
# 1. Define Crossover Key
# ============================================================

labels_f <- c(
  "NCO_0_f","NCO_0_f",
  "SCO_1_f","SCO_1_f","SCO_2_f","SCO_2_f","SCO_3_f","SCO_3_f","SCO_4_f","SCO_4_f","SCO_5_f","SCO_5_f",
  "DCO_1,2_f","DCO_1,2_f","DCO_1,3_f","DCO_1,3_f","DCO_1,4_f","DCO_1,4_f","DCO_1,5_f","DCO_1,5_f",
  "DCO_2,3_f","DCO_2,3_f","DCO_2,4_f","DCO_2,4_f","DCO_2,5_f","DCO_2,5_f","DCO_3,4_f","DCO_3,4_f",
  "DCO_3,5_f","DCO_3,5_f","DCO_4,5_f","DCO_4,5_f",
  "TCO_1,2,3_f","TCO_1,2,3_f","TCO_1,2,4_f","TCO_1,2,4_f","TCO_1,2,5_f","TCO_1,2,5_f",
  "TCO_1,3,4_f","TCO_1,3,4_f","TCO_1,3,5_f","TCO_1,3,5_f","TCO_1,4,5_f","TCO_1,4,5_f",
  "TCO_2,3,4_f","TCO_2,3,4_f","TCO_2,3,5_f","TCO_2,3,5_f","TCO_2,4,5_f","TCO_2,4,5_f",
  "TCO_3,4,5_f","TCO_3,4,5_f",
  "QCO_1,2,3,4_f","QCO_1,2,3,4_f","QCO_1,2,3,5_f","QCO_1,2,3,5_f",
  "QCO_1,2,4,5_f","QCO_1,2,4,5_f","QCO_1,3,4,5_f","QCO_1,3,4,5_f",
  "QCO_2,3,4,5_f","QCO_2,3,4,5_f",
  "FCO_1,2,3,4,5_f","FCO_1,2,3,4,5_f"
)
labels_m <- gsub("_f", "_m", labels_f)

crossover_key <- data.frame(
  class_num = 1:128,
  crossover_type = c(labels_f, labels_m),
  stringsAsFactors = FALSE
)

# ============================================================
# 2. Interval Mapping & Recombinant Flag
# ============================================================

interval_membership <- data.frame(class_num = 1:128, I1=0, I2=0, I3=0, I4=0, I5=0, is_recombinant=0)

for (i in 1:nrow(crossover_key)) {
  # Extract the numbers
  nums <- as.integer(regmatches(crossover_key$crossover_type[i], 
                                gregexpr("[0-9]+", crossover_key$crossover_type[i]))[[1]])
  
  # If any number > 0, it's a recombinant fly
  if (any(nums > 0)) {
    interval_membership$is_recombinant[i] <- 1
    
    # Flip the mapping for intervals as requested: 5->1 (y-cv), 1->5 (car-y+)
    flipped_nums <- 6 - nums
    for (val in flipped_nums) {
      interval_membership[i, paste0("I", val)] <- 1
    }
  }
}

# ============================================================
# 3. Processing and Grouping
# ============================================================

class_cols <- paste0("class_", 1:128)

counts_long <- multi_locus_data %>%
  select(ID, all_of(class_cols)) %>%
  pivot_longer(cols = all_of(class_cols), names_to = "class_col", values_to = "count") %>%
  mutate(class_num = as.integer(gsub("class_", "", class_col))) %>%
  left_join(crossover_key, by = "class_num") %>%
  left_join(interval_membership, by = "class_num")

# Per-vial Calculation
recomb_by_interval_vial <- counts_long %>%
  group_by(ID) %>%
  summarise(
    total_flies = sum(count, na.rm = TRUE),
    total_recombinants = sum(count[is_recombinant == 1], na.rm = TRUE),
    across(I1:I5, ~sum(count[.x == 1], na.rm = TRUE), .names = "CO_{.col}"),
    .groups = "drop"
  ) %>%
  mutate(
    Total_RF   = total_recombinants / total_flies,
    RF_y_cv    = CO_I1 / total_flies,
    RF_cv_v    = CO_I2 / total_flies,
    RF_v_f     = CO_I3 / total_flies,
    RF_f_car   = CO_I4 / total_flies,
    RF_car_yP  = CO_I5 / total_flies
  )


# ============================================================
# 4. Final Output
# ============================================================

recomb_by_interval_vial <- recomb_by_interval_vial %>% extract(ID,
                                                               into = c("experiment", "vialID", "Brood"),
                                                               regex = "^(.*)_([0-9]+)_([A-Za-z])$"
)

write.csv(recomb_by_interval_vial, "outputs/14_recombination_frequency_per_vial.csv", row.names = FALSE)
 

