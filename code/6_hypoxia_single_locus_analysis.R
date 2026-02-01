#setwd("Raw Data/") # Please use Raw Data as working directory

# Perform single-locus goodness-of-fit G-tests for scute calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- scute_single_locus_dataset
scute_G_test_table <- single_locus_analysis()
write.csv(scute_G_test_table,"../Summaries and Outputs/table_S3_scute_single_locus_G_tests.csv",row.names = TRUE)

# Perform single-locus goodness-of-fit G-tests for crossveinless calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- crossveinless_single_locus_dataset
crossveinless_G_test_table <- single_locus_analysis()
write.csv(crossveinless_G_test_table,"../Summaries and Outputs/table_S4_crossveinless_single_locus_G_tests.csv",row.names = TRUE)

# Perform single-locus goodness-of-fit G-tests for vermilion calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- vermilion_single_locus_dataset
vermilion_G_test_table <- single_locus_analysis()
write.csv(vermilion_G_test_table,"../Summaries and Outputs/table_S5_vermilion_single_locus_G_tests.csv",row.names = TRUE)

# Perform single-locus goodness-of-fit G-tests for forked calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- forked_single_locus_dataset
forked_G_test_table <- single_locus_analysis()
write.csv(forked_G_test_table,"../Summaries and Outputs/table_S6_forked_single_locus_G_tests.csv",row.names = TRUE)

# Perform single-locus goodness-of-fit G-tests for carnation calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- carnation_single_locus_dataset
carnation_G_test_table <- single_locus_analysis()
write.csv(carnation_G_test_table,"../Summaries and Outputs/table_S7_carnation_single_locus_G_tests.csv",row.names = TRUE)

# Perform single-locus goodness-of-fit G-tests for yellow plus calling custom function from "single_locus_analysis_custom_function.R"
single_locus_input <- yellow_plus_single_locus_dataset
yellow_plus_G_test_table <- single_locus_analysis()
write.csv(yellow_plus_G_test_table,"../Summaries and Outputs/table_S8_yellow_plus_single_locus_G_tests.csv",row.names = TRUE)
