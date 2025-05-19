# merge data together

'%ni%' <- Negate('%in%')
library(dplyr)

load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_manipulated_data.RData")

load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/pan_compiled_data.RData")
pan_compiled_data <- pan_compiled_data[!duplicated(pan_compiled_data$s),]

veg_table_merge <- read.delim("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/VegIIDs.txt")

id_bridge_merge <- read.delim("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/ukb48818bridge31063.txt", header = FALSE, sep = " ")
colnames(id_bridge_merge) <- c("FID", "s")

load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/townsend_data.RData")

id_pan <- left_join(pan_compiled_data, id_bridge_merge, by = "s")
id_man_pan <- left_join(veg_manipulated_data, id_pan, by = "FID")
id_man_pan_veg <- left_join(id_man_pan, veg_table_merge, by = "FID")
veg_merged_data <- left_join(id_man_pan_veg, townsend_data, by = "FID")

save(veg_merged_data, file = "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_merged_data.RData")


