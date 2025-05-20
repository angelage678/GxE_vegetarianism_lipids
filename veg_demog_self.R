# extracting demographic information for self-identified vegetarian

'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library("lubridate")

Pathway=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files")
Pathway_out=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files")

load(paste("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/","veg_excluded_data_self.RData",sep = ""))

bd_pheno_self <- data.frame(veg_excluded_data_self)

#factorize covariates
bd_pheno_self$Smoking_status_cov <- as.factor(bd_pheno_self$smoking_standardized)
bd_pheno_self$Smoking_status_cov<- factor(bd_pheno_self$Smoking_status_cov, levels = c("Never","Previous","Current"))
bd_pheno_self$Alcohol_status_cov <- as.factor(bd_pheno_self$alcohol_standardized)
bd_pheno_self$Alcohol_status_cov<- factor(bd_pheno_self$Alcohol_status_cov, levels = c("Never","Previous","Current"))
bd_pheno_self$Physical_activity_cov <- as.factor(bd_pheno_self$physical_activity_standardized)
bd_pheno_self$Physical_activity_cov<- factor(bd_pheno_self$Physical_activity_cov, levels = c("Low","Moderate","High"))

#self
#Split into ancestry then split into ancestry w/ veg/not
EUR_self <- bd_pheno_self[bd_pheno_self$pop == "EUR", ]
EUR_self_veg <- EUR_self[EUR_self$SelfID == 1, ]
EUR_self_not <- EUR_self[EUR_self$SelfID == 0, ]
EAS_self <- bd_pheno_self[bd_pheno_self$pop == "EAS", ]
EAS_self_veg <- EAS_self[EAS_self$SelfID == 1, ]
EAS_self_not <- EAS_self[EAS_self$SelfID == 0, ]
CSA_self <- bd_pheno_self[bd_pheno_self$pop == "CSA", ]
CSA_self_veg <- CSA_self[CSA_self$SelfID == 1, ]
CSA_self_not <- CSA_self[CSA_self$SelfID == 0, ]
AFR_self <- bd_pheno_self[bd_pheno_self$pop == "AFR", ]
AFR_self_veg <- AFR_self[AFR_self$SelfID == 1, ]
AFR_self_not <- AFR_self[AFR_self$SelfID == 0, ]

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
                  paste("Vegetarian"," (n=",(nrow(EUR_self_veg)),")", sep=""),
                  paste(format(round(mean(EUR_self_veg$Age),digit=0))," (",(round(sd(EUR_self_veg$Age),digit=1)), ")",sep=""),
                  paste((nrow(EUR_self_veg[EUR_self_veg$Sex == 0, ]))," (",(round(nrow(EUR_self_veg[EUR_self_veg$Sex == 0, ])/(nrow(EUR_self_veg))*100)), ")",sep=""),
                  paste(format(round(mean(EUR_self_veg[!is.na(EUR_self_veg$BMI), ]$BMI),digit=0))," (",(round(sd(EUR_self_veg[!is.na(EUR_self_veg$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(EUR_self_veg[!is.na(EUR_self_veg$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(EUR_self_veg[!is.na(EUR_self_veg$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_self_veg[!is.na(EUR_self_veg$LDL), ]$LDL),digit=2))," (",(round(sd(EUR_self_veg[!is.na(EUR_self_veg$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_self_veg[!is.na(EUR_self_veg$HDL), ]$HDL),digit=2))," (",(round(sd(EUR_self_veg[!is.na(EUR_self_veg$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_self_veg[!is.na(EUR_self_veg$TAGs), ]$TAGs),digit=2))," (",(round(sd(EUR_self_veg[!is.na(EUR_self_veg$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(EUR_self_veg$smoking_standardized)[2])," (",(round(table(EUR_self_veg$smoking_standardized)[2]/(nrow(EUR_self_veg))*100)), ")",sep=""),
                  paste((table(EUR_self_veg$smoking_standardized)[3])," (",(round(table(EUR_self_veg$smoking_standardized)[3]/(nrow(EUR_self_veg))*100)), ")",sep=""),
                  paste((table(EUR_self_veg$smoking_standardized)[1])," (",(round(table(EUR_self_veg$smoking_standardized)[1]/(nrow(EUR_self_veg))*100)), ")",sep=""),
                  "",
                  paste((table(EUR_self_veg$alcohol_standardized)[2])," (",(round(table(EUR_self_veg$alcohol_standardized)[2]/(nrow(EUR_self_veg))*100)), ")",sep=""),
                  paste((table(EUR_self_veg$alcohol_standardized)[3])," (",(round(table(EUR_self_veg$alcohol_standardized)[3]/(nrow(EUR_self_veg))*100)), ")",sep=""),
                  paste((table(EUR_self_veg$alcohol_standardized)[1])," (",(round(table(EUR_self_veg$alcohol_standardized)[1]/(nrow(EUR_self_veg))*100)), ")",sep=""),
                  "",
                  paste((table(EUR_self_veg$physical_activity_standardized)[2])," (",(round(table(EUR_self_veg$physical_activity_standardized)[2]/(nrow(EUR_self_veg))*100)), ")",sep=""),
                  paste((table(EUR_self_veg$physical_activity_standardized)[3])," (",(round(table(EUR_self_veg$physical_activity_standardized)[3]/(nrow(EUR_self_veg))*100)), ")",sep=""),
                  paste((table(EUR_self_veg$physical_activity_standardized)[1])," (",(round(table(EUR_self_veg$physical_activity_standardized)[1]/(nrow(EUR_self_veg))*100)), ")",sep=""),
                  (paste((nrow(EUR_self_veg[EUR_self_veg$statin_sum_standardized == 1, ]))," (",(round(nrow(EUR_self_veg[EUR_self_veg$statin_sum_standardized == 1, ])/(nrow(EUR_self_veg))*100)), ")",sep="")))

EUR_not_demog = c("European",
                  paste("Not Vegetarian"," (n=",(nrow(EUR_self_not)),")", sep=""),
                  paste(format(round(mean(EUR_self_not$Age),digit=0))," (",(round(sd(EUR_self_not$Age),digit=1)), ")",sep=""),
                  paste((nrow(EUR_self_not[EUR_self_not$Sex == 0, ]))," (",(round(nrow(EUR_self_not[EUR_self_not$Sex == 0, ])/(nrow(EUR_self_not))*100)), ")",sep=""),
                  paste(format(round(mean(EUR_self_not[!is.na(EUR_self_not$BMI), ]$BMI),digit=0))," (",(round(sd(EUR_self_not[!is.na(EUR_self_not$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(EUR_self_not[!is.na(EUR_self_not$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(EUR_self_not[!is.na(EUR_self_not$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_self_not[!is.na(EUR_self_not$LDL), ]$LDL),digit=2))," (",(round(sd(EUR_self_not[!is.na(EUR_self_not$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_self_not[!is.na(EUR_self_not$HDL), ]$HDL),digit=2))," (",(round(sd(EUR_self_not[!is.na(EUR_self_not$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_self_not[!is.na(EUR_self_not$TAGs), ]$TAGs),digit=2))," (",(round(sd(EUR_self_not[!is.na(EUR_self_not$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(EUR_self_not$smoking_standardized)[2])," (",(round(table(EUR_self_not$smoking_standardized)[2]/(nrow(EUR_self_not))*100)), ")",sep=""),
                  paste((table(EUR_self_not$smoking_standardized)[3])," (",(round(table(EUR_self_not$smoking_standardized)[3]/(nrow(EUR_self_not))*100)), ")",sep=""),
                  paste((table(EUR_self_not$smoking_standardized)[1])," (",(round(table(EUR_self_not$smoking_standardized)[1]/(nrow(EUR_self_not))*100)), ")",sep=""),
                  "",
                  paste((table(EUR_self_not$alcohol_standardized)[2])," (",(round(table(EUR_self_not$alcohol_standardized)[2]/(nrow(EUR_self_not))*100)), ")",sep=""),
                  paste((table(EUR_self_not$alcohol_standardized)[3])," (",(round(table(EUR_self_not$alcohol_standardized)[3]/(nrow(EUR_self_not))*100)), ")",sep=""),
                  paste((table(EUR_self_not$alcohol_standardized)[1])," (",(round(table(EUR_self_not$alcohol_standardized)[1]/(nrow(EUR_self_not))*100)), ")",sep=""),
                  "",
                  paste((table(EUR_self_not$physical_activity_standardized)[2])," (",(round(table(EUR_self_not$physical_activity_standardized)[2]/(nrow(EUR_self_not))*100)), ")",sep=""),
                  paste((table(EUR_self_not$physical_activity_standardized)[3])," (",(round(table(EUR_self_not$physical_activity_standardized)[3]/(nrow(EUR_self_not))*100)), ")",sep=""),
                  paste((table(EUR_self_not$physical_activity_standardized)[1])," (",(round(table(EUR_self_not$physical_activity_standardized)[1]/(nrow(EUR_self_not))*100)), ")",sep=""),
                  (paste((nrow(EUR_self_not[EUR_self_not$statin_sum_standardized == 1, ]))," (",(round(nrow(EUR_self_not[EUR_self_not$statin_sum_standardized == 1, ])/(nrow(EUR_self_not))*100)), ")",sep="")))

EAS_veg_demog = c("East Asian",
                  paste("Vegetarian"," (n=",(nrow(EAS_self_veg)),")", sep=""),
                  paste(format(round(mean(EAS_self_veg$Age),digit=0))," (",(round(sd(EAS_self_veg$Age),digit=1)), ")",sep=""),
                  paste((nrow(EAS_self_veg[EAS_self_veg$Sex == 0, ]))," (",(round(nrow(EAS_self_veg[EAS_self_veg$Sex == 0, ])/(nrow(EAS_self_veg))*100)), ")",sep=""),
                  paste(format(round(mean(EAS_self_veg[!is.na(EAS_self_veg$BMI), ]$BMI),digit=0))," (",(round(sd(EAS_self_veg[!is.na(EAS_self_veg$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(EAS_self_veg[!is.na(EAS_self_veg$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(EAS_self_veg[!is.na(EAS_self_veg$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_self_veg[!is.na(EAS_self_veg$LDL), ]$LDL),digit=2))," (",(round(sd(EAS_self_veg[!is.na(EAS_self_veg$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_self_veg[!is.na(EAS_self_veg$HDL), ]$HDL),digit=2))," (",(round(sd(EAS_self_veg[!is.na(EAS_self_veg$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_self_veg[!is.na(EAS_self_veg$TAGs), ]$TAGs),digit=2))," (",(round(sd(EAS_self_veg[!is.na(EAS_self_veg$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(EAS_self_veg$smoking_standardized)[2])," (",(round(table(EAS_self_veg$smoking_standardized)[2]/(nrow(EAS_self_veg))*100)), ")",sep=""),
                  paste((table(EAS_self_veg$smoking_standardized)[3])," (",(round(table(EAS_self_veg$smoking_standardized)[3]/(nrow(EAS_self_veg))*100)), ")",sep=""),
                  paste((table(EAS_self_veg$smoking_standardized)[1])," (",(round(table(EAS_self_veg$smoking_standardized)[1]/(nrow(EAS_self_veg))*100)), ")",sep=""),
                  "",
                  paste((table(EAS_self_veg$alcohol_standardized)[2])," (",(round(table(EAS_self_veg$alcohol_standardized)[2]/(nrow(EAS_self_veg))*100)), ")",sep=""),
                  paste((table(EAS_self_veg$alcohol_standardized)[3])," (",(round(table(EAS_self_veg$alcohol_standardized)[3]/(nrow(EAS_self_veg))*100)), ")",sep=""),
                  paste((table(EAS_self_veg$alcohol_standardized)[1])," (",(round(table(EAS_self_veg$alcohol_standardized)[1]/(nrow(EAS_self_veg))*100)), ")",sep=""),
                  "",
                  paste((table(EAS_self_veg$physical_activity_standardized)[2])," (",(round(table(EAS_self_veg$physical_activity_standardized)[2]/(nrow(EAS_self_veg))*100)), ")",sep=""),
                  paste((table(EAS_self_veg$physical_activity_standardized)[3])," (",(round(table(EAS_self_veg$physical_activity_standardized)[3]/(nrow(EAS_self_veg))*100)), ")",sep=""),
                  paste((table(EAS_self_veg$physical_activity_standardized)[1])," (",(round(table(EAS_self_veg$physical_activity_standardized)[1]/(nrow(EAS_self_veg))*100)), ")",sep=""),
                  (paste((nrow(EAS_self_veg[EAS_self_veg$statin_sum_standardized == 1, ]))," (",(round(nrow(EAS_self_veg[EAS_self_veg$statin_sum_standardized == 1, ])/(nrow(EAS_self_veg))*100)), ")",sep="")))

EAS_not_demog = c("East Asian",
                  paste("Not Vegetarian"," (n=",(nrow(EAS_self_not)),")", sep=""),
                  paste(format(round(mean(EAS_self_not$Age),digit=0))," (",(round(sd(EAS_self_not$Age),digit=1)), ")",sep=""),
                  paste((nrow(EAS_self_not[EAS_self_not$Sex == 0, ]))," (",(round(nrow(EAS_self_not[EAS_self_not$Sex == 0, ])/(nrow(EAS_self_not))*100)), ")",sep=""),
                  paste(format(round(mean(EAS_self_not[!is.na(EAS_self_not$BMI), ]$BMI),digit=0))," (",(round(sd(EAS_self_not[!is.na(EAS_self_not$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(EAS_self_not[!is.na(EAS_self_not$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(EAS_self_not[!is.na(EAS_self_not$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_self_not[!is.na(EAS_self_not$LDL), ]$LDL),digit=2))," (",(round(sd(EAS_self_not[!is.na(EAS_self_not$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_self_not[!is.na(EAS_self_not$HDL), ]$HDL),digit=2))," (",(round(sd(EAS_self_not[!is.na(EAS_self_not$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_self_not[!is.na(EAS_self_not$TAGs), ]$TAGs),digit=2))," (",(round(sd(EAS_self_not[!is.na(EAS_self_not$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(EAS_self_not$smoking_standardized)[2])," (",(round(table(EAS_self_not$smoking_standardized)[2]/(nrow(EAS_self_not))*100)), ")",sep=""),
                  paste((table(EAS_self_not$smoking_standardized)[3])," (",(round(table(EAS_self_not$smoking_standardized)[3]/(nrow(EAS_self_not))*100)), ")",sep=""),
                  paste((table(EAS_self_not$smoking_standardized)[1])," (",(round(table(EAS_self_not$smoking_standardized)[1]/(nrow(EAS_self_not))*100)), ")",sep=""),
                  "",
                  paste((table(EAS_self_not$alcohol_standardized)[2])," (",(round(table(EAS_self_not$alcohol_standardized)[2]/(nrow(EAS_self_not))*100)), ")",sep=""),
                  paste((table(EAS_self_not$alcohol_standardized)[3])," (",(round(table(EAS_self_not$alcohol_standardized)[3]/(nrow(EAS_self_not))*100)), ")",sep=""),
                  paste((table(EAS_self_not$alcohol_standardized)[1])," (",(round(table(EAS_self_not$alcohol_standardized)[1]/(nrow(EAS_self_not))*100)), ")",sep=""),
                  "",
                  paste((table(EAS_self_not$physical_activity_standardized)[2])," (",(round(table(EAS_self_not$physical_activity_standardized)[2]/(nrow(EAS_self_not))*100)), ")",sep=""),
                  paste((table(EAS_self_not$physical_activity_standardized)[3])," (",(round(table(EAS_self_not$physical_activity_standardized)[3]/(nrow(EAS_self_not))*100)), ")",sep=""),
                  paste((table(EAS_self_not$physical_activity_standardized)[1])," (",(round(table(EAS_self_not$physical_activity_standardized)[1]/(nrow(EAS_self_not))*100)), ")",sep=""),
                  (paste((nrow(EAS_self_not[EAS_self_not$statin_sum_standardized == 1, ]))," (",(round(nrow(EAS_self_not[EAS_self_not$statin_sum_standardized == 1, ])/(nrow(EAS_self_not))*100)), ")",sep="")))

CSA_veg_demog = c("Central/South Asian",
                  paste("Vegetarian"," (n=",(nrow(CSA_self_veg)),")", sep=""),
                  paste(format(round(mean(CSA_self_veg$Age),digit=0))," (",(round(sd(CSA_self_veg$Age),digit=1)), ")",sep=""),
                  paste((nrow(CSA_self_veg[CSA_self_veg$Sex == 0, ]))," (",(round(nrow(CSA_self_veg[CSA_self_veg$Sex == 0, ])/(nrow(CSA_self_veg))*100)), ")",sep=""),
                  paste(format(round(mean(CSA_self_veg[!is.na(CSA_self_veg$BMI), ]$BMI),digit=0))," (",(round(sd(CSA_self_veg[!is.na(CSA_self_veg$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(CSA_self_veg[!is.na(CSA_self_veg$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(CSA_self_veg[!is.na(CSA_self_veg$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_self_veg[!is.na(CSA_self_veg$LDL), ]$LDL),digit=2))," (",(round(sd(CSA_self_veg[!is.na(CSA_self_veg$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_self_veg[!is.na(CSA_self_veg$HDL), ]$HDL),digit=2))," (",(round(sd(CSA_self_veg[!is.na(CSA_self_veg$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_self_veg[!is.na(CSA_self_veg$TAGs), ]$TAGs),digit=2))," (",(round(sd(CSA_self_veg[!is.na(CSA_self_veg$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(CSA_self_veg$smoking_standardized)[2])," (",(round(table(CSA_self_veg$smoking_standardized)[2]/(nrow(CSA_self_veg))*100)), ")",sep=""),
                  paste((table(CSA_self_veg$smoking_standardized)[3])," (",(round(table(CSA_self_veg$smoking_standardized)[3]/(nrow(CSA_self_veg))*100)), ")",sep=""),
                  paste((table(CSA_self_veg$smoking_standardized)[1])," (",(round(table(CSA_self_veg$smoking_standardized)[1]/(nrow(CSA_self_veg))*100)), ")",sep=""),
                  "",
                  paste((table(CSA_self_veg$alcohol_standardized)[2])," (",(round(table(CSA_self_veg$alcohol_standardized)[2]/(nrow(CSA_self_veg))*100)), ")",sep=""),
                  paste((table(CSA_self_veg$alcohol_standardized)[3])," (",(round(table(CSA_self_veg$alcohol_standardized)[3]/(nrow(CSA_self_veg))*100)), ")",sep=""),
                  paste((table(CSA_self_veg$alcohol_standardized)[1])," (",(round(table(CSA_self_veg$alcohol_standardized)[1]/(nrow(CSA_self_veg))*100)), ")",sep=""),
                  "",
                  paste((table(CSA_self_veg$physical_activity_standardized)[2])," (",(round(table(CSA_self_veg$physical_activity_standardized)[2]/(nrow(CSA_self_veg))*100)), ")",sep=""),
                  paste((table(CSA_self_veg$physical_activity_standardized)[3])," (",(round(table(CSA_self_veg$physical_activity_standardized)[3]/(nrow(CSA_self_veg))*100)), ")",sep=""),
                  paste((table(CSA_self_veg$physical_activity_standardized)[1])," (",(round(table(CSA_self_veg$physical_activity_standardized)[1]/(nrow(CSA_self_veg))*100)), ")",sep=""),
                  (paste((nrow(CSA_self_veg[CSA_self_veg$statin_sum_standardized == 1, ]))," (",(round(nrow(CSA_self_veg[CSA_self_veg$statin_sum_standardized == 1, ])/(nrow(CSA_self_veg))*100)), ")",sep="")))

CSA_not_demog = c("Central/South Asian",
                  paste("Not Vegetarian"," (n=",(nrow(CSA_self_not)),")", sep=""),
                  paste(format(round(mean(CSA_self_not$Age),digit=0))," (",(round(sd(CSA_self_not$Age),digit=1)), ")",sep=""),
                  paste((nrow(CSA_self_not[CSA_self_not$Sex == 0, ]))," (",(round(nrow(CSA_self_not[CSA_self_not$Sex == 0, ])/(nrow(CSA_self_not))*100)), ")",sep=""),
                  paste(format(round(mean(CSA_self_not[!is.na(CSA_self_not$BMI), ]$BMI),digit=0))," (",(round(sd(CSA_self_not[!is.na(CSA_self_not$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(CSA_self_not[!is.na(CSA_self_not$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(CSA_self_not[!is.na(CSA_self_not$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_self_not[!is.na(CSA_self_not$LDL), ]$LDL),digit=2))," (",(round(sd(CSA_self_not[!is.na(CSA_self_not$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_self_not[!is.na(CSA_self_not$HDL), ]$HDL),digit=2))," (",(round(sd(CSA_self_not[!is.na(CSA_self_not$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_self_not[!is.na(CSA_self_not$TAGs), ]$TAGs),digit=2))," (",(round(sd(CSA_self_not[!is.na(CSA_self_not$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(CSA_self_not$smoking_standardized)[2])," (",(round(table(CSA_self_not$smoking_standardized)[2]/(nrow(CSA_self_not))*100)), ")",sep=""),
                  paste((table(CSA_self_not$smoking_standardized)[3])," (",(round(table(CSA_self_not$smoking_standardized)[3]/(nrow(CSA_self_not))*100)), ")",sep=""),
                  paste((table(CSA_self_not$smoking_standardized)[1])," (",(round(table(CSA_self_not$smoking_standardized)[1]/(nrow(CSA_self_not))*100)), ")",sep=""),
                  "",
                  paste((table(CSA_self_not$alcohol_standardized)[2])," (",(round(table(CSA_self_not$alcohol_standardized)[2]/(nrow(CSA_self_not))*100)), ")",sep=""),
                  paste((table(CSA_self_not$alcohol_standardized)[3])," (",(round(table(CSA_self_not$alcohol_standardized)[3]/(nrow(CSA_self_not))*100)), ")",sep=""),
                  paste((table(CSA_self_not$alcohol_standardized)[1])," (",(round(table(CSA_self_not$alcohol_standardized)[1]/(nrow(CSA_self_not))*100)), ")",sep=""),
                  "",
                  paste((table(CSA_self_not$physical_activity_standardized)[2])," (",(round(table(CSA_self_not$physical_activity_standardized)[2]/(nrow(CSA_self_not))*100)), ")",sep=""),
                  paste((table(CSA_self_not$physical_activity_standardized)[3])," (",(round(table(CSA_self_not$physical_activity_standardized)[3]/(nrow(CSA_self_not))*100)), ")",sep=""),
                  paste((table(CSA_self_not$physical_activity_standardized)[1])," (",(round(table(CSA_self_not$physical_activity_standardized)[1]/(nrow(CSA_self_not))*100)), ")",sep=""),
                  (paste((nrow(CSA_self_not[CSA_self_not$statin_sum_standardized == 1, ]))," (",(round(nrow(CSA_self_not[CSA_self_not$statin_sum_standardized == 1, ])/(nrow(CSA_self_not))*100)), ")",sep="")))

AFR_veg_demog = c("African",
                  paste("Vegetarian"," (n=",(nrow(AFR_self_veg)),")", sep=""),
                  paste(format(round(mean(AFR_self_veg$Age),digit=0))," (",(round(sd(AFR_self_veg$Age),digit=1)), ")",sep=""),
                  paste((nrow(AFR_self_veg[AFR_self_veg$Sex == 0, ]))," (",(round(nrow(AFR_self_veg[AFR_self_veg$Sex == 0, ])/(nrow(AFR_self_veg))*100)), ")",sep=""),
                  paste(format(round(mean(AFR_self_veg[!is.na(AFR_self_veg$BMI), ]$BMI),digit=0))," (",(round(sd(AFR_self_veg[!is.na(AFR_self_veg$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(AFR_self_veg[!is.na(AFR_self_veg$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(AFR_self_veg[!is.na(AFR_self_veg$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_self_veg[!is.na(AFR_self_veg$LDL), ]$LDL),digit=2))," (",(round(sd(AFR_self_veg[!is.na(AFR_self_veg$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_self_veg[!is.na(AFR_self_veg$HDL), ]$HDL),digit=2))," (",(round(sd(AFR_self_veg[!is.na(AFR_self_veg$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_self_veg[!is.na(AFR_self_veg$TAGs), ]$TAGs),digit=2))," (",(round(sd(AFR_self_veg[!is.na(AFR_self_veg$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(AFR_self_veg$smoking_standardized)[2])," (",(round(table(AFR_self_veg$smoking_standardized)[2]/(nrow(AFR_self_veg))*100)), ")",sep=""),
                  paste((table(AFR_self_veg$smoking_standardized)[3])," (",(round(table(AFR_self_veg$smoking_standardized)[3]/(nrow(AFR_self_veg))*100)), ")",sep=""),
                  paste((table(AFR_self_veg$smoking_standardized)[1])," (",(round(table(AFR_self_veg$smoking_standardized)[1]/(nrow(AFR_self_veg))*100)), ")",sep=""),
                  "",
                  paste((table(AFR_self_veg$alcohol_standardized)[2])," (",(round(table(AFR_self_veg$alcohol_standardized)[2]/(nrow(AFR_self_veg))*100)), ")",sep=""),
                  paste((table(AFR_self_veg$alcohol_standardized)[3])," (",(round(table(AFR_self_veg$alcohol_standardized)[3]/(nrow(AFR_self_veg))*100)), ")",sep=""),
                  paste((table(AFR_self_veg$alcohol_standardized)[1])," (",(round(table(AFR_self_veg$alcohol_standardized)[1]/(nrow(AFR_self_veg))*100)), ")",sep=""),
                  "",
                  paste((table(AFR_self_veg$physical_activity_standardized)[2])," (",(round(table(AFR_self_veg$physical_activity_standardized)[2]/(nrow(AFR_self_veg))*100)), ")",sep=""),
                  paste((table(AFR_self_veg$physical_activity_standardized)[3])," (",(round(table(AFR_self_veg$physical_activity_standardized)[3]/(nrow(AFR_self_veg))*100)), ")",sep=""),
                  paste((table(AFR_self_veg$physical_activity_standardized)[1])," (",(round(table(AFR_self_veg$physical_activity_standardized)[1]/(nrow(AFR_self_veg))*100)), ")",sep=""),
                  (paste((nrow(AFR_self_veg[AFR_self_veg$statin_sum_standardized == 1, ]))," (",(round(nrow(AFR_self_veg[AFR_self_veg$statin_sum_standardized == 1, ])/(nrow(AFR_self_veg))*100)), ")",sep="")))

AFR_not_demog = c("African",
                  paste("Not Vegetarian"," (n=",(nrow(AFR_self_not)),")", sep=""),
                  paste(format(round(mean(AFR_self_not$Age),digit=0))," (",(round(sd(AFR_self_not$Age),digit=1)), ")",sep=""),
                  paste((nrow(AFR_self_not[AFR_self_not$Sex == 0, ]))," (",(round(nrow(AFR_self_not[AFR_self_not$Sex == 0, ])/(nrow(AFR_self_not))*100)), ")",sep=""),
                  paste(format(round(mean(AFR_self_not[!is.na(AFR_self_not$BMI), ]$BMI),digit=0))," (",(round(sd(AFR_self_not[!is.na(AFR_self_not$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(AFR_self_not[!is.na(AFR_self_not$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(AFR_self_not[!is.na(AFR_self_not$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_self_not[!is.na(AFR_self_not$LDL), ]$LDL),digit=2))," (",(round(sd(AFR_self_not[!is.na(AFR_self_not$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_self_not[!is.na(AFR_self_not$HDL), ]$HDL),digit=2))," (",(round(sd(AFR_self_not[!is.na(AFR_self_not$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_self_not[!is.na(AFR_self_not$TAGs), ]$TAGs),digit=2))," (",(round(sd(AFR_self_not[!is.na(AFR_self_not$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(AFR_self_not$smoking_standardized)[2])," (",(round(table(AFR_self_not$smoking_standardized)[2]/(nrow(AFR_self_not))*100)), ")",sep=""),
                  paste((table(AFR_self_not$smoking_standardized)[3])," (",(round(table(AFR_self_not$smoking_standardized)[3]/(nrow(AFR_self_not))*100)), ")",sep=""),
                  paste((table(AFR_self_not$smoking_standardized)[1])," (",(round(table(AFR_self_not$smoking_standardized)[1]/(nrow(AFR_self_not))*100)), ")",sep=""),
                  "",
                  paste((table(AFR_self_not$alcohol_standardized)[2])," (",(round(table(AFR_self_not$alcohol_standardized)[2]/(nrow(AFR_self_not))*100)), ")",sep=""),
                  paste((table(AFR_self_not$alcohol_standardized)[3])," (",(round(table(AFR_self_not$alcohol_standardized)[3]/(nrow(AFR_self_not))*100)), ")",sep=""),
                  paste((table(AFR_self_not$alcohol_standardized)[1])," (",(round(table(AFR_self_not$alcohol_standardized)[1]/(nrow(AFR_self_not))*100)), ")",sep=""),
                  "",
                  paste((table(AFR_self_not$physical_activity_standardized)[2])," (",(round(table(AFR_self_not$physical_activity_standardized)[2]/(nrow(AFR_self_not))*100)), ")",sep=""),
                  paste((table(AFR_self_not$physical_activity_standardized)[3])," (",(round(table(AFR_self_not$physical_activity_standardized)[3]/(nrow(AFR_self_not))*100)), ")",sep=""),
                  paste((table(AFR_self_not$physical_activity_standardized)[1])," (",(round(table(AFR_self_not$physical_activity_standardized)[1]/(nrow(AFR_self_not))*100)), ")",sep=""),
                  (paste((nrow(AFR_self_not[AFR_self_not$statin_sum_standardized == 1, ]))," (",(round(nrow(AFR_self_not[AFR_self_not$statin_sum_standardized == 1, ])/(nrow(AFR_self_not))*100)), ")",sep="")))



veg_demographics <- cbind(row_label, EUR_veg_demog, EUR_not_demog, EAS_veg_demog, EAS_not_demog, CSA_veg_demog, CSA_not_demog, AFR_veg_demog, AFR_not_demog)


save(veg_demographics, file = "veg_demographics.RData")

write.table(veg_demographics, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_demographics_self.txt", col.names = FALSE, append = TRUE,
            row.names = F, quote = FALSE, na = "",sep='\t')


































