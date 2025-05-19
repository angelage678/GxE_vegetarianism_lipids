# cleaning/excluding data

'%ni%' <- Negate('%in%')
library(dplyr)

load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_bd_pheno.RData")
veg_excluded_data <- data.frame(veg_bd_pheno)

#withdrawn
withdrawn <- read.csv("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/w48818_2023-04-25.csv", header = FALSE)
veg_excluded_data <- veg_excluded_data[veg_excluded_data$FID %ni% withdrawn$V1, ]

#sex don't match
veg_excluded_data <- veg_excluded_data[veg_excluded_data$Sex == veg_excluded_data$Genetic_sex, ]

#sex chrom aneuploidy, is.na bc only two types of values are 1 and NA and don't want 1
veg_excluded_data <- veg_excluded_data[is.na(veg_excluded_data$Sex_chrom_aneuploidy), ]

#outlier, same reasoning as above
veg_excluded_data <- veg_excluded_data[is.na(veg_excluded_data$Outlier_het_miss_genotype), ]

#high degree kinship
veg_excluded_data <- veg_excluded_data[veg_excluded_data$Genetic_kinship != 10, ]

#missing pan ancestry
veg_excluded_data <- veg_excluded_data[!is.na(veg_excluded_data$pop), ]

#remove vegetarian data missing
veg_excluded_data_self <- veg_excluded_data[!(is.na(veg_excluded_data$SelfID)), ]
veg_excluded_data_strict <- veg_excluded_data[!(is.na(veg_excluded_data$Strict)), ]

save(veg_excluded_data, file = "veg_excluded_data.RData")
save(veg_excluded_data_strict, file = "veg_excluded_data_strict.RData")
save(veg_excluded_data_self, file = "veg_excluded_data_self.RData")
