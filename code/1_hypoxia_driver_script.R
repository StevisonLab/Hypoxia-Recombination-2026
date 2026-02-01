setwd("Raw Data/") # Please use Raw Data as working directory

# Compile, format, transform all raw data into appropriate input files
source(file = "../Code/2_missing_data_dataset_compilation_script.R")

# Load custom functions for single and multi-locus analyses
source(file = "../Code/3_missing_data_single_locus_custom_functions.R")
source(file = "../Code/4_missing_data_multi_locus_custom_functions.R")

# Load package dfoptim
library(dfoptim)

# Perform organismal viability analyses
source(file = "../Code/5_missing_data_organismal_analysis.R")

# Perform single locus viability analyses
source(file = "../Code/6_missing_data_single_locus_analysis.R")

# Perform multi-locus viability analyses
# This script will take a LONG time!
source(file = "../Code/7_missing_data_multi_locus_analysis.R")

# Perform alternate pooling viability analyses
# This script will also take a LONG time!
source(file = "../Code/8_missing_data_pooling_analysis.R")
