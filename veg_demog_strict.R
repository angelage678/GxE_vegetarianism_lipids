# extracting demographic information for strict vegetarian

'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library("lubridate")

Pathway=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files")
Pathway_out=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files")

load(paste("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/","veg_excluded_data_strict_PRS.RData",sep = ""))

bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS)

#factorize covariates
bd_pheno_strict$Smoking_status_cov <- as.factor(bd_pheno_strict$smoking_standardized)
bd_pheno_strict$Smoking_status_cov<- factor(bd_pheno_strict$Smoking_status_cov, levels = c("Never","Previous","Current"))
bd_pheno_strict$Alcohol_status_cov <- as.factor(bd_pheno_strict$alcohol_standardized)
bd_pheno_strict$Alcohol_status_cov<- factor(bd_pheno_strict$Alcohol_status_cov, levels = c("Never","Previous","Current"))
bd_pheno_strict$Physical_activity_cov <- as.factor(bd_pheno_strict$physical_activity_standardized)
bd_pheno_strict$Physical_activity_cov<- factor(bd_pheno_strict$Physical_activity_cov, levels = c("Low","Moderate","High"))

#STRICT
#Split into ancestry then split into ancestry w/ veg/not
EUR_strict <- bd_pheno_strict[bd_pheno_strict$pop == "EUR", ]
EUR_strict_veg <- EUR_strict[EUR_strict$Strict == 1, ]
EUR_strict_not <- EUR_strict[EUR_strict$Strict == 0, ]
EAS_strict <- bd_pheno_strict[bd_pheno_strict$pop == "EAS", ]
EAS_strict_veg <- EAS_strict[EAS_strict$Strict == 1, ]
EAS_strict_not <- EAS_strict[EAS_strict$Strict == 0, ]
CSA_strict <- bd_pheno_strict[bd_pheno_strict$pop == "CSA", ]
CSA_strict_veg <- CSA_strict[CSA_strict$Strict == 1, ]
CSA_strict_not <- CSA_strict[CSA_strict$Strict == 0, ]
AFR_strict <- bd_pheno_strict[bd_pheno_strict$pop == "AFR", ]
AFR_strict_veg <- AFR_strict[AFR_strict$Strict == 1, ]
AFR_strict_not <- AFR_strict[AFR_strict$Strict == 0, ]

#Age, Sex, BMI, TC, LDL, HDL, Trig, Smoking (Never, Previous, Current), Alcohol (Never, Previous, Current), Physical activity (Low, Moderate, High), Statin use (yes)
row_label <- c("",
               "",
               "Age, years (SD)", 
               "Sex, female (%)", 
               "Body mass index, kg/m^2 (SD)", 
               "Total cholesterol, mmol/L (SD)", 
               "LDL cholesterol, mmol/L (SD)", 
               "HDL cholesterol, mmol/L (SD)", 
               "Triglycerides, mmol/L (SD)", 
               "Smoking status (%)",
               "    Never",
               "    Previous",
               "    Current",
               "Alcohol status (%)",
               "    Never",
               "    Previous",
               "    Current",
               "Physical activity (%)",
               "    Low",
               "    Moderate",
               "    High",
               "Statin use, yes (%)")
EUR_veg_demog = c("European",
                  paste("Vegetarian"," (n=",(nrow(EUR_strict_veg)),")", sep=""),
               paste(format(round(mean(EUR_strict_veg$Age),digit=0))," (",(round(sd(EUR_strict_veg$Age),digit=1)), ")",sep=""),
               paste((nrow(EUR_strict_veg[EUR_strict_veg$Sex == 0, ]))," (",(round(nrow(EUR_strict_veg[EUR_strict_veg$Sex == 0, ])/(nrow(EUR_strict_veg))*100)), ")",sep=""),
               paste(format(round(mean(EUR_strict_veg[!is.na(EUR_strict_veg$BMI), ]$BMI),digit=0))," (",(round(sd(EUR_strict_veg[!is.na(EUR_strict_veg$BMI), ]$BMI),digit=1)), ")",sep=""),
               paste(format(round(mean(EUR_strict_veg[!is.na(EUR_strict_veg$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(EUR_strict_veg[!is.na(EUR_strict_veg$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
               paste(format(round(mean(EUR_strict_veg[!is.na(EUR_strict_veg$LDL), ]$LDL),digit=2))," (",(round(sd(EUR_strict_veg[!is.na(EUR_strict_veg$LDL), ]$LDL),digit=2)), ")",sep=""),
               paste(format(round(mean(EUR_strict_veg[!is.na(EUR_strict_veg$HDL), ]$HDL),digit=2))," (",(round(sd(EUR_strict_veg[!is.na(EUR_strict_veg$HDL), ]$HDL),digit=2)), ")",sep=""),
               paste(format(round(mean(EUR_strict_veg[!is.na(EUR_strict_veg$TAGs), ]$TAGs),digit=2))," (",(round(sd(EUR_strict_veg[!is.na(EUR_strict_veg$TAGs), ]$TAGs),digit=2)), ")",sep=""),
               "",
               paste((table(EUR_strict_veg$smoking_standardized)[2])," (",(round(table(EUR_strict_veg$smoking_standardized)[2]/(nrow(EUR_strict_veg))*100)), ")",sep=""),
               paste((table(EUR_strict_veg$smoking_standardized)[3])," (",(round(table(EUR_strict_veg$smoking_standardized)[3]/(nrow(EUR_strict_veg))*100)), ")",sep=""),
               paste((table(EUR_strict_veg$smoking_standardized)[1])," (",(round(table(EUR_strict_veg$smoking_standardized)[1]/(nrow(EUR_strict_veg))*100)), ")",sep=""),
               "",
               paste((table(EUR_strict_veg$alcohol_standardized)[2])," (",(round(table(EUR_strict_veg$alcohol_standardized)[2]/(nrow(EUR_strict_veg))*100)), ")",sep=""),
               paste((table(EUR_strict_veg$alcohol_standardized)[3])," (",(round(table(EUR_strict_veg$alcohol_standardized)[3]/(nrow(EUR_strict_veg))*100)), ")",sep=""),
               paste((table(EUR_strict_veg$alcohol_standardized)[1])," (",(round(table(EUR_strict_veg$alcohol_standardized)[1]/(nrow(EUR_strict_veg))*100)), ")",sep=""),
               "",
               paste((table(EUR_strict_veg$physical_activity_standardized)[2])," (",(round(table(EUR_strict_veg$physical_activity_standardized)[2]/(nrow(EUR_strict_veg))*100)), ")",sep=""),
               paste((table(EUR_strict_veg$physical_activity_standardized)[3])," (",(round(table(EUR_strict_veg$physical_activity_standardized)[3]/(nrow(EUR_strict_veg))*100)), ")",sep=""),
               paste((table(EUR_strict_veg$physical_activity_standardized)[1])," (",(round(table(EUR_strict_veg$physical_activity_standardized)[1]/(nrow(EUR_strict_veg))*100)), ")",sep=""),
               (paste((nrow(EUR_strict_veg[EUR_strict_veg$statin_sum_standardized == 1, ]))," (",(round(nrow(EUR_strict_veg[EUR_strict_veg$statin_sum_standardized == 1, ])/(nrow(EUR_strict_veg))*100)), ")",sep="")))

EUR_not_demog = c("European",
                  paste("Not Vegetarian"," (n=",(nrow(EUR_strict_not)),")", sep=""),
                  paste(format(round(mean(EUR_strict_not$Age),digit=0))," (",(round(sd(EUR_strict_not$Age),digit=1)), ")",sep=""),
                  paste((nrow(EUR_strict_not[EUR_strict_not$Sex == 0, ]))," (",(round(nrow(EUR_strict_not[EUR_strict_not$Sex == 0, ])/(nrow(EUR_strict_not))*100)), ")",sep=""),
                  paste(format(round(mean(EUR_strict_not[!is.na(EUR_strict_not$BMI), ]$BMI),digit=0))," (",(round(sd(EUR_strict_not[!is.na(EUR_strict_not$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(EUR_strict_not[!is.na(EUR_strict_not$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(EUR_strict_not[!is.na(EUR_strict_not$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_strict_not[!is.na(EUR_strict_not$LDL), ]$LDL),digit=2))," (",(round(sd(EUR_strict_not[!is.na(EUR_strict_not$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_strict_not[!is.na(EUR_strict_not$HDL), ]$HDL),digit=2))," (",(round(sd(EUR_strict_not[!is.na(EUR_strict_not$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_strict_not[!is.na(EUR_strict_not$TAGs), ]$TAGs),digit=2))," (",(round(sd(EUR_strict_not[!is.na(EUR_strict_not$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(EUR_strict_not$smoking_standardized)[2])," (",(round(table(EUR_strict_not$smoking_standardized)[2]/(nrow(EUR_strict_not))*100)), ")",sep=""),
                  paste((table(EUR_strict_not$smoking_standardized)[3])," (",(round(table(EUR_strict_not$smoking_standardized)[3]/(nrow(EUR_strict_not))*100)), ")",sep=""),
                  paste((table(EUR_strict_not$smoking_standardized)[1])," (",(round(table(EUR_strict_not$smoking_standardized)[1]/(nrow(EUR_strict_not))*100)), ")",sep=""),
                  "",
                  paste((table(EUR_strict_not$alcohol_standardized)[2])," (",(round(table(EUR_strict_not$alcohol_standardized)[2]/(nrow(EUR_strict_not))*100)), ")",sep=""),
                  paste((table(EUR_strict_not$alcohol_standardized)[3])," (",(round(table(EUR_strict_not$alcohol_standardized)[3]/(nrow(EUR_strict_not))*100)), ")",sep=""),
                  paste((table(EUR_strict_not$alcohol_standardized)[1])," (",(round(table(EUR_strict_not$alcohol_standardized)[1]/(nrow(EUR_strict_not))*100)), ")",sep=""),
                  "",
                  paste((table(EUR_strict_not$physical_activity_standardized)[2])," (",(round(table(EUR_strict_not$physical_activity_standardized)[2]/(nrow(EUR_strict_not))*100)), ")",sep=""),
                  paste((table(EUR_strict_not$physical_activity_standardized)[3])," (",(round(table(EUR_strict_not$physical_activity_standardized)[3]/(nrow(EUR_strict_not))*100)), ")",sep=""),
                  paste((table(EUR_strict_not$physical_activity_standardized)[1])," (",(round(table(EUR_strict_not$physical_activity_standardized)[1]/(nrow(EUR_strict_not))*100)), ")",sep=""),
                  (paste((nrow(EUR_strict_not[EUR_strict_not$statin_sum_standardized == 1, ]))," (",(round(nrow(EUR_strict_not[EUR_strict_not$statin_sum_standardized == 1, ])/(nrow(EUR_strict_not))*100)), ")",sep="")))

EAS_veg_demog = c("East Asian",
                  paste("Vegetarian"," (n=",(nrow(EAS_strict_veg)),")", sep=""),
                  paste(format(round(mean(EAS_strict_veg$Age),digit=0))," (",(round(sd(EAS_strict_veg$Age),digit=1)), ")",sep=""),
                  paste((nrow(EAS_strict_veg[EAS_strict_veg$Sex == 0, ]))," (",(round(nrow(EAS_strict_veg[EAS_strict_veg$Sex == 0, ])/(nrow(EAS_strict_veg))*100)), ")",sep=""),
                  paste(format(round(mean(EAS_strict_veg[!is.na(EAS_strict_veg$BMI), ]$BMI),digit=0))," (",(round(sd(EAS_strict_veg[!is.na(EAS_strict_veg$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(EAS_strict_veg[!is.na(EAS_strict_veg$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(EAS_strict_veg[!is.na(EAS_strict_veg$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_strict_veg[!is.na(EAS_strict_veg$LDL), ]$LDL),digit=2))," (",(round(sd(EAS_strict_veg[!is.na(EAS_strict_veg$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_strict_veg[!is.na(EAS_strict_veg$HDL), ]$HDL),digit=2))," (",(round(sd(EAS_strict_veg[!is.na(EAS_strict_veg$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_strict_veg[!is.na(EAS_strict_veg$TAGs), ]$TAGs),digit=2))," (",(round(sd(EAS_strict_veg[!is.na(EAS_strict_veg$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(EAS_strict_veg$smoking_standardized)[2])," (",(round(table(EAS_strict_veg$smoking_standardized)[2]/(nrow(EAS_strict_veg))*100)), ")",sep=""),
                  paste((table(EAS_strict_veg$smoking_standardized)[3])," (",(round(table(EAS_strict_veg$smoking_standardized)[3]/(nrow(EAS_strict_veg))*100)), ")",sep=""),
                  paste((table(EAS_strict_veg$smoking_standardized)[1])," (",(round(table(EAS_strict_veg$smoking_standardized)[1]/(nrow(EAS_strict_veg))*100)), ")",sep=""),
                  "",
                  paste((table(EAS_strict_veg$alcohol_standardized)[2])," (",(round(table(EAS_strict_veg$alcohol_standardized)[2]/(nrow(EAS_strict_veg))*100)), ")",sep=""),
                  paste((table(EAS_strict_veg$alcohol_standardized)[3])," (",(round(table(EAS_strict_veg$alcohol_standardized)[3]/(nrow(EAS_strict_veg))*100)), ")",sep=""),
                  paste((table(EAS_strict_veg$alcohol_standardized)[1])," (",(round(table(EAS_strict_veg$alcohol_standardized)[1]/(nrow(EAS_strict_veg))*100)), ")",sep=""),
                  "",
                  paste((table(EAS_strict_veg$physical_activity_standardized)[2])," (",(round(table(EAS_strict_veg$physical_activity_standardized)[2]/(nrow(EAS_strict_veg))*100)), ")",sep=""),
                  paste((table(EAS_strict_veg$physical_activity_standardized)[3])," (",(round(table(EAS_strict_veg$physical_activity_standardized)[3]/(nrow(EAS_strict_veg))*100)), ")",sep=""),
                  paste((table(EAS_strict_veg$physical_activity_standardized)[1])," (",(round(table(EAS_strict_veg$physical_activity_standardized)[1]/(nrow(EAS_strict_veg))*100)), ")",sep=""),
                  (paste((nrow(EAS_strict_veg[EAS_strict_veg$statin_sum_standardized == 1, ]))," (",(round(nrow(EAS_strict_veg[EAS_strict_veg$statin_sum_standardized == 1, ])/(nrow(EAS_strict_veg))*100)), ")",sep="")))

EAS_not_demog = c("East Asian",
                  paste("Not Vegetarian"," (n=",(nrow(EAS_strict_not)),")", sep=""),
                  paste(format(round(mean(EAS_strict_not$Age),digit=0))," (",(round(sd(EAS_strict_not$Age),digit=1)), ")",sep=""),
                  paste((nrow(EAS_strict_not[EAS_strict_not$Sex == 0, ]))," (",(round(nrow(EAS_strict_not[EAS_strict_not$Sex == 0, ])/(nrow(EAS_strict_not))*100)), ")",sep=""),
                  paste(format(round(mean(EAS_strict_not[!is.na(EAS_strict_not$BMI), ]$BMI),digit=0))," (",(round(sd(EAS_strict_not[!is.na(EAS_strict_not$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(EAS_strict_not[!is.na(EAS_strict_not$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(EAS_strict_not[!is.na(EAS_strict_not$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_strict_not[!is.na(EAS_strict_not$LDL), ]$LDL),digit=2))," (",(round(sd(EAS_strict_not[!is.na(EAS_strict_not$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_strict_not[!is.na(EAS_strict_not$HDL), ]$HDL),digit=2))," (",(round(sd(EAS_strict_not[!is.na(EAS_strict_not$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_strict_not[!is.na(EAS_strict_not$TAGs), ]$TAGs),digit=2))," (",(round(sd(EAS_strict_not[!is.na(EAS_strict_not$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(EAS_strict_not$smoking_standardized)[2])," (",(round(table(EAS_strict_not$smoking_standardized)[2]/(nrow(EAS_strict_not))*100)), ")",sep=""),
                  paste((table(EAS_strict_not$smoking_standardized)[3])," (",(round(table(EAS_strict_not$smoking_standardized)[3]/(nrow(EAS_strict_not))*100)), ")",sep=""),
                  paste((table(EAS_strict_not$smoking_standardized)[1])," (",(round(table(EAS_strict_not$smoking_standardized)[1]/(nrow(EAS_strict_not))*100)), ")",sep=""),
                  "",
                  paste((table(EAS_strict_not$alcohol_standardized)[2])," (",(round(table(EAS_strict_not$alcohol_standardized)[2]/(nrow(EAS_strict_not))*100)), ")",sep=""),
                  paste((table(EAS_strict_not$alcohol_standardized)[3])," (",(round(table(EAS_strict_not$alcohol_standardized)[3]/(nrow(EAS_strict_not))*100)), ")",sep=""),
                  paste((table(EAS_strict_not$alcohol_standardized)[1])," (",(round(table(EAS_strict_not$alcohol_standardized)[1]/(nrow(EAS_strict_not))*100)), ")",sep=""),
                  "",
                  paste((table(EAS_strict_not$physical_activity_standardized)[2])," (",(round(table(EAS_strict_not$physical_activity_standardized)[2]/(nrow(EAS_strict_not))*100)), ")",sep=""),
                  paste((table(EAS_strict_not$physical_activity_standardized)[3])," (",(round(table(EAS_strict_not$physical_activity_standardized)[3]/(nrow(EAS_strict_not))*100)), ")",sep=""),
                  paste((table(EAS_strict_not$physical_activity_standardized)[1])," (",(round(table(EAS_strict_not$physical_activity_standardized)[1]/(nrow(EAS_strict_not))*100)), ")",sep=""),
                  (paste((nrow(EAS_strict_not[EAS_strict_not$statin_sum_standardized == 1, ]))," (",(round(nrow(EAS_strict_not[EAS_strict_not$statin_sum_standardized == 1, ])/(nrow(EAS_strict_not))*100)), ")",sep="")))

CSA_veg_demog = c("Central/South Asian",
                  paste("Vegetarian"," (n=",(nrow(CSA_strict_veg)),")", sep=""),
                  paste(format(round(mean(CSA_strict_veg$Age),digit=0))," (",(round(sd(CSA_strict_veg$Age),digit=1)), ")",sep=""),
                  paste((nrow(CSA_strict_veg[CSA_strict_veg$Sex == 0, ]))," (",(round(nrow(CSA_strict_veg[CSA_strict_veg$Sex == 0, ])/(nrow(CSA_strict_veg))*100)), ")",sep=""),
                  paste(format(round(mean(CSA_strict_veg[!is.na(CSA_strict_veg$BMI), ]$BMI),digit=0))," (",(round(sd(CSA_strict_veg[!is.na(CSA_strict_veg$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(CSA_strict_veg[!is.na(CSA_strict_veg$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(CSA_strict_veg[!is.na(CSA_strict_veg$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_strict_veg[!is.na(CSA_strict_veg$LDL), ]$LDL),digit=2))," (",(round(sd(CSA_strict_veg[!is.na(CSA_strict_veg$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_strict_veg[!is.na(CSA_strict_veg$HDL), ]$HDL),digit=2))," (",(round(sd(CSA_strict_veg[!is.na(CSA_strict_veg$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_strict_veg[!is.na(CSA_strict_veg$TAGs), ]$TAGs),digit=2))," (",(round(sd(CSA_strict_veg[!is.na(CSA_strict_veg$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(CSA_strict_veg$smoking_standardized)[2])," (",(round(table(CSA_strict_veg$smoking_standardized)[2]/(nrow(CSA_strict_veg))*100)), ")",sep=""),
                  paste((table(CSA_strict_veg$smoking_standardized)[3])," (",(round(table(CSA_strict_veg$smoking_standardized)[3]/(nrow(CSA_strict_veg))*100)), ")",sep=""),
                  paste((table(CSA_strict_veg$smoking_standardized)[1])," (",(round(table(CSA_strict_veg$smoking_standardized)[1]/(nrow(CSA_strict_veg))*100)), ")",sep=""),
                  "",
                  paste((table(CSA_strict_veg$alcohol_standardized)[2])," (",(round(table(CSA_strict_veg$alcohol_standardized)[2]/(nrow(CSA_strict_veg))*100)), ")",sep=""),
                  paste((table(CSA_strict_veg$alcohol_standardized)[3])," (",(round(table(CSA_strict_veg$alcohol_standardized)[3]/(nrow(CSA_strict_veg))*100)), ")",sep=""),
                  paste((table(CSA_strict_veg$alcohol_standardized)[1])," (",(round(table(CSA_strict_veg$alcohol_standardized)[1]/(nrow(CSA_strict_veg))*100)), ")",sep=""),
                  "",
                  paste((table(CSA_strict_veg$physical_activity_standardized)[2])," (",(round(table(CSA_strict_veg$physical_activity_standardized)[2]/(nrow(CSA_strict_veg))*100)), ")",sep=""),
                  paste((table(CSA_strict_veg$physical_activity_standardized)[3])," (",(round(table(CSA_strict_veg$physical_activity_standardized)[3]/(nrow(CSA_strict_veg))*100)), ")",sep=""),
                  paste((table(CSA_strict_veg$physical_activity_standardized)[1])," (",(round(table(CSA_strict_veg$physical_activity_standardized)[1]/(nrow(CSA_strict_veg))*100)), ")",sep=""),
                  (paste((nrow(CSA_strict_veg[CSA_strict_veg$statin_sum_standardized == 1, ]))," (",(round(nrow(CSA_strict_veg[CSA_strict_veg$statin_sum_standardized == 1, ])/(nrow(CSA_strict_veg))*100)), ")",sep="")))

CSA_not_demog = c("Central/South Asian",
                  paste("Not Vegetarian"," (n=",(nrow(CSA_strict_not)),")", sep=""),
                  paste(format(round(mean(CSA_strict_not$Age),digit=0))," (",(round(sd(CSA_strict_not$Age),digit=1)), ")",sep=""),
                  paste((nrow(CSA_strict_not[CSA_strict_not$Sex == 0, ]))," (",(round(nrow(CSA_strict_not[CSA_strict_not$Sex == 0, ])/(nrow(CSA_strict_not))*100)), ")",sep=""),
                  paste(format(round(mean(CSA_strict_not[!is.na(CSA_strict_not$BMI), ]$BMI),digit=0))," (",(round(sd(CSA_strict_not[!is.na(CSA_strict_not$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(CSA_strict_not[!is.na(CSA_strict_not$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(CSA_strict_not[!is.na(CSA_strict_not$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_strict_not[!is.na(CSA_strict_not$LDL), ]$LDL),digit=2))," (",(round(sd(CSA_strict_not[!is.na(CSA_strict_not$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_strict_not[!is.na(CSA_strict_not$HDL), ]$HDL),digit=2))," (",(round(sd(CSA_strict_not[!is.na(CSA_strict_not$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_strict_not[!is.na(CSA_strict_not$TAGs), ]$TAGs),digit=2))," (",(round(sd(CSA_strict_not[!is.na(CSA_strict_not$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(CSA_strict_not$smoking_standardized)[2])," (",(round(table(CSA_strict_not$smoking_standardized)[2]/(nrow(CSA_strict_not))*100)), ")",sep=""),
                  paste((table(CSA_strict_not$smoking_standardized)[3])," (",(round(table(CSA_strict_not$smoking_standardized)[3]/(nrow(CSA_strict_not))*100)), ")",sep=""),
                  paste((table(CSA_strict_not$smoking_standardized)[1])," (",(round(table(CSA_strict_not$smoking_standardized)[1]/(nrow(CSA_strict_not))*100)), ")",sep=""),
                  "",
                  paste((table(CSA_strict_not$alcohol_standardized)[2])," (",(round(table(CSA_strict_not$alcohol_standardized)[2]/(nrow(CSA_strict_not))*100)), ")",sep=""),
                  paste((table(CSA_strict_not$alcohol_standardized)[3])," (",(round(table(CSA_strict_not$alcohol_standardized)[3]/(nrow(CSA_strict_not))*100)), ")",sep=""),
                  paste((table(CSA_strict_not$alcohol_standardized)[1])," (",(round(table(CSA_strict_not$alcohol_standardized)[1]/(nrow(CSA_strict_not))*100)), ")",sep=""),
                  "",
                  paste((table(CSA_strict_not$physical_activity_standardized)[2])," (",(round(table(CSA_strict_not$physical_activity_standardized)[2]/(nrow(CSA_strict_not))*100)), ")",sep=""),
                  paste((table(CSA_strict_not$physical_activity_standardized)[3])," (",(round(table(CSA_strict_not$physical_activity_standardized)[3]/(nrow(CSA_strict_not))*100)), ")",sep=""),
                  paste((table(CSA_strict_not$physical_activity_standardized)[1])," (",(round(table(CSA_strict_not$physical_activity_standardized)[1]/(nrow(CSA_strict_not))*100)), ")",sep=""),
                  (paste((nrow(CSA_strict_not[CSA_strict_not$statin_sum_standardized == 1, ]))," (",(round(nrow(CSA_strict_not[CSA_strict_not$statin_sum_standardized == 1, ])/(nrow(CSA_strict_not))*100)), ")",sep="")))

AFR_veg_demog = c("African",
                  paste("Vegetarian"," (n=",(nrow(AFR_strict_veg)),")", sep=""),
                  paste(format(round(mean(AFR_strict_veg$Age),digit=0))," (",(round(sd(AFR_strict_veg$Age),digit=1)), ")",sep=""),
                  paste((nrow(AFR_strict_veg[AFR_strict_veg$Sex == 0, ]))," (",(round(nrow(AFR_strict_veg[AFR_strict_veg$Sex == 0, ])/(nrow(AFR_strict_veg))*100)), ")",sep=""),
                  paste(format(round(mean(AFR_strict_veg[!is.na(AFR_strict_veg$BMI), ]$BMI),digit=0))," (",(round(sd(AFR_strict_veg[!is.na(AFR_strict_veg$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(AFR_strict_veg[!is.na(AFR_strict_veg$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(AFR_strict_veg[!is.na(AFR_strict_veg$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_strict_veg[!is.na(AFR_strict_veg$LDL), ]$LDL),digit=2))," (",(round(sd(AFR_strict_veg[!is.na(AFR_strict_veg$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_strict_veg[!is.na(AFR_strict_veg$HDL), ]$HDL),digit=2))," (",(round(sd(AFR_strict_veg[!is.na(AFR_strict_veg$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_strict_veg[!is.na(AFR_strict_veg$TAGs), ]$TAGs),digit=2))," (",(round(sd(AFR_strict_veg[!is.na(AFR_strict_veg$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(AFR_strict_veg$smoking_standardized)[2])," (",(round(table(AFR_strict_veg$smoking_standardized)[2]/(nrow(AFR_strict_veg))*100)), ")",sep=""),
                  paste((table(AFR_strict_veg$smoking_standardized)[3])," (",(round(table(AFR_strict_veg$smoking_standardized)[3]/(nrow(AFR_strict_veg))*100)), ")",sep=""),
                  paste((table(AFR_strict_veg$smoking_standardized)[1])," (",(round(table(AFR_strict_veg$smoking_standardized)[1]/(nrow(AFR_strict_veg))*100)), ")",sep=""),
                  "",
                  paste((table(AFR_strict_veg$alcohol_standardized)[2])," (",(round(table(AFR_strict_veg$alcohol_standardized)[2]/(nrow(AFR_strict_veg))*100)), ")",sep=""),
                  paste((table(AFR_strict_veg$alcohol_standardized)[3])," (",(round(table(AFR_strict_veg$alcohol_standardized)[3]/(nrow(AFR_strict_veg))*100)), ")",sep=""),
                  paste((table(AFR_strict_veg$alcohol_standardized)[1])," (",(round(table(AFR_strict_veg$alcohol_standardized)[1]/(nrow(AFR_strict_veg))*100)), ")",sep=""),
                  "",
                  paste((table(AFR_strict_veg$physical_activity_standardized)[2])," (",(round(table(AFR_strict_veg$physical_activity_standardized)[2]/(nrow(AFR_strict_veg))*100)), ")",sep=""),
                  paste((table(AFR_strict_veg$physical_activity_standardized)[3])," (",(round(table(AFR_strict_veg$physical_activity_standardized)[3]/(nrow(AFR_strict_veg))*100)), ")",sep=""),
                  paste((table(AFR_strict_veg$physical_activity_standardized)[1])," (",(round(table(AFR_strict_veg$physical_activity_standardized)[1]/(nrow(AFR_strict_veg))*100)), ")",sep=""),
                  (paste((nrow(AFR_strict_veg[AFR_strict_veg$statin_sum_standardized == 1, ]))," (",(round(nrow(AFR_strict_veg[AFR_strict_veg$statin_sum_standardized == 1, ])/(nrow(AFR_strict_veg))*100)), ")",sep="")))

AFR_not_demog = c("African",
                  paste("Not Vegetarian"," (n=",(nrow(AFR_strict_not)),")", sep=""),
                  paste(format(round(mean(AFR_strict_not$Age),digit=0))," (",(round(sd(AFR_strict_not$Age),digit=1)), ")",sep=""),
                  paste((nrow(AFR_strict_not[AFR_strict_not$Sex == 0, ]))," (",(round(nrow(AFR_strict_not[AFR_strict_not$Sex == 0, ])/(nrow(AFR_strict_not))*100)), ")",sep=""),
                  paste(format(round(mean(AFR_strict_not[!is.na(AFR_strict_not$BMI), ]$BMI),digit=0))," (",(round(sd(AFR_strict_not[!is.na(AFR_strict_not$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(AFR_strict_not[!is.na(AFR_strict_not$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(AFR_strict_not[!is.na(AFR_strict_not$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_strict_not[!is.na(AFR_strict_not$LDL), ]$LDL),digit=2))," (",(round(sd(AFR_strict_not[!is.na(AFR_strict_not$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_strict_not[!is.na(AFR_strict_not$HDL), ]$HDL),digit=2))," (",(round(sd(AFR_strict_not[!is.na(AFR_strict_not$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_strict_not[!is.na(AFR_strict_not$TAGs), ]$TAGs),digit=2))," (",(round(sd(AFR_strict_not[!is.na(AFR_strict_not$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(AFR_strict_not$smoking_standardized)[2])," (",(round(table(AFR_strict_not$smoking_standardized)[2]/(nrow(AFR_strict_not))*100)), ")",sep=""),
                  paste((table(AFR_strict_not$smoking_standardized)[3])," (",(round(table(AFR_strict_not$smoking_standardized)[3]/(nrow(AFR_strict_not))*100)), ")",sep=""),
                  paste((table(AFR_strict_not$smoking_standardized)[1])," (",(round(table(AFR_strict_not$smoking_standardized)[1]/(nrow(AFR_strict_not))*100)), ")",sep=""),
                  "",
                  paste((table(AFR_strict_not$alcohol_standardized)[2])," (",(round(table(AFR_strict_not$alcohol_standardized)[2]/(nrow(AFR_strict_not))*100)), ")",sep=""),
                  paste((table(AFR_strict_not$alcohol_standardized)[3])," (",(round(table(AFR_strict_not$alcohol_standardized)[3]/(nrow(AFR_strict_not))*100)), ")",sep=""),
                  paste((table(AFR_strict_not$alcohol_standardized)[1])," (",(round(table(AFR_strict_not$alcohol_standardized)[1]/(nrow(AFR_strict_not))*100)), ")",sep=""),
                  "",
                  paste((table(AFR_strict_not$physical_activity_standardized)[2])," (",(round(table(AFR_strict_not$physical_activity_standardized)[2]/(nrow(AFR_strict_not))*100)), ")",sep=""),
                  paste((table(AFR_strict_not$physical_activity_standardized)[3])," (",(round(table(AFR_strict_not$physical_activity_standardized)[3]/(nrow(AFR_strict_not))*100)), ")",sep=""),
                  paste((table(AFR_strict_not$physical_activity_standardized)[1])," (",(round(table(AFR_strict_not$physical_activity_standardized)[1]/(nrow(AFR_strict_not))*100)), ")",sep=""),
                  (paste((nrow(AFR_strict_not[AFR_strict_not$statin_sum_standardized == 1, ]))," (",(round(nrow(AFR_strict_not[AFR_strict_not$statin_sum_standardized == 1, ])/(nrow(AFR_strict_not))*100)), ")",sep="")))


veg_demographics <- cbind(row_label, EUR_veg_demog, EUR_not_demog, EAS_veg_demog, EAS_not_demog, CSA_veg_demog, CSA_not_demog, AFR_veg_demog, AFR_not_demog)


save(veg_demographics, file = "veg_demographics.RData")

write.table(veg_demographics, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_demographics.txt", col.names = FALSE, append = TRUE,
            row.names = F, quote = FALSE, na = "",sep='\t')


































