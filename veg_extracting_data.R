# Extracts data from UKBB

'%ni%' <- Negate('%in%')
library(dplyr)

#read in the data
ukbb_data <- read.table("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/ukb672926.tab", header = TRUE, sep = "\t", fill = TRUE)
#save the data as an RData file
save(ukbb_data, file = "ukbb_data_2.RData")
#read in the RData file
load("ukbb_data_2.RData")
#extra data
veg_extracted_data <- ukbb_data %>% select(f.eid,
                                          f.30690.0.0, f.30780.0.0, f.30760.0.0, f.30870.0.0,
                                          f.31.0.0, f.21003.0.0, f.22000.0.0, f.54.0.0, f.21001.0.0, f.20117.0.0, f.20116.0.0, f.22032.0.0,
                                          f.20003.0.0, f.20003.0.1, f.20003.0.2, f.20003.0.3, f.20003.0.4, f.20003.0.5, f.20003.0.6, f.20003.0.7, f.20003.0.8, f.20003.0.9, f.20003.0.10, 
                                          f.20003.0.11, f.20003.0.12, f.20003.0.13, f.20003.0.14, f.20003.0.15, f.20003.0.16, f.20003.0.17, f.20003.0.18, f.20003.0.19, f.20003.0.20, 
                                          f.20003.0.21, f.20003.0.22, f.20003.0.23, f.20003.0.24, f.20003.0.25, f.20003.0.26, f.20003.0.27, f.20003.0.28, f.20003.0.29, f.20003.0.30, 
                                          f.20003.0.31, f.20003.0.32, f.20003.0.33, f.20003.0.34, f.20003.0.35, f.20003.0.36, f.20003.0.37, f.20003.0.38, f.20003.0.39, f.20003.0.40, 
                                          f.20003.0.41, f.20003.0.42, f.20003.0.43, f.20003.0.44, f.20003.0.45, f.20003.0.46, f.20003.0.47, 
                                          f.22001.0.0, f.22019.0.0, f.22027.0.0, f.22021.0.0)
#modify column names
colnames(veg_extracted_data) <- c(#id 
  "FID", 
  #lipids
  "Tot_Chol", "LDL", "HDL", "TAGs", 
  #covariates 
  "Sex", "Age", "Array", "Assessment_centres", "BMI", "Alcohol_status", "Smoking_status", "Physical_activity", 
  "Statin0", "Statin1", "Statin2", "Statin3", "Statin4", "Statin5", "Statin6", "Statin7", "Statin8", "Statin9", "Statin10", 
  "Statin11", "Statin12", "Statin13", "Statin14", "Statin15", "Statin16", "Statin17", "Statin18", "Statin19", "Statin20", 
  "Statin21", "Statin22", "Statin23", "Statin24", "Statin25", "Statin26", "Statin27", "Statin28", "Statin29", "Statin30", 
  "Statin31", "Statin32", "Statin33", "Statin34", "Statin35", "Statin36", "Statin37", "Statin38", "Statin39", "Statin40", 
  "Statin41", "Statin42", "Statin43", "Statin44", "Statin45", "Statin46", "Statin47", 
  #exclusion criteria
  "Genetic_sex", "Sex_chrom_aneuploidy", "Outlier_het_miss_genotype", "Genetic_kinship")
#save extracted data
save(veg_extracted_data, file = "veg_extracted_data.RData")