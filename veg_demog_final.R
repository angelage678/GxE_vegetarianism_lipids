# extracting demographic information for main table

'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library("lubridate")

Pathway=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files")
Pathway_out=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files")

#using self file because it includes both self and strict participants
load(paste("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/","veg_excluded_data_self.RData",sep = ""))

bd_pheno <- data.frame(veg_excluded_data_self)

#factorize covariates
bd_pheno$Smoking_status_cov <- as.factor(bd_pheno$smoking_standardized)
bd_pheno$Smoking_status_cov<- factor(bd_pheno$Smoking_status_cov, levels = c("Never","Previous","Current"))
bd_pheno$Alcohol_status_cov <- as.factor(bd_pheno$alcohol_standardized)
bd_pheno$Alcohol_status_cov<- factor(bd_pheno$Alcohol_status_cov, levels = c("Never","Previous","Current"))
bd_pheno$Physical_activity_cov <- as.factor(bd_pheno$physical_activity_standardized)
bd_pheno$Physical_activity_cov<- factor(bd_pheno$Physical_activity_cov, levels = c("Low","Moderate","High"))

#Split into ancestry 
EUR_data <- bd_pheno[bd_pheno$pop == "EUR", ]
EAS_data <- bd_pheno[bd_pheno$pop == "EAS", ]
CSA_data <- bd_pheno[bd_pheno$pop == "CSA", ]
AFR_data <- bd_pheno[bd_pheno$pop == "AFR", ]

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
               "Vegetarian, yes (%)", #this is self
               "    Strict (%)", 
               "    Missing (%)", 
               "Statin use, yes (%)")
EUR_demog = c("European",
                  paste("Vegetarian"," (n=",(nrow(EUR_data)),")", sep=""),
                  paste(format(round(mean(EUR_data$Age),digit=0))," (",(round(sd(EUR_data$Age),digit=1)), ")",sep=""),
                  paste((nrow(EUR_data[EUR_data$Sex == 0, ]))," (",(round(nrow(EUR_data[EUR_data$Sex == 0, ])/(nrow(EUR_data))*100)), ")",sep=""),
                  paste(format(round(mean(EUR_data[!is.na(EUR_data$BMI), ]$BMI),digit=0))," (",(round(sd(EUR_data[!is.na(EUR_data$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(EUR_data[!is.na(EUR_data$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(EUR_data[!is.na(EUR_data$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_data[!is.na(EUR_data$LDL), ]$LDL),digit=2))," (",(round(sd(EUR_data[!is.na(EUR_data$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_data[!is.na(EUR_data$HDL), ]$HDL),digit=2))," (",(round(sd(EUR_data[!is.na(EUR_data$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EUR_data[!is.na(EUR_data$TAGs), ]$TAGs),digit=2))," (",(round(sd(EUR_data[!is.na(EUR_data$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(EUR_data$smoking_standardized)[2])," (",(round(table(EUR_data$smoking_standardized)[2]/(nrow(EUR_data))*100)), ")",sep=""),
                  paste((table(EUR_data$smoking_standardized)[3])," (",(round(table(EUR_data$smoking_standardized)[3]/(nrow(EUR_data))*100)), ")",sep=""),
                  paste((table(EUR_data$smoking_standardized)[1])," (",(round(table(EUR_data$smoking_standardized)[1]/(nrow(EUR_data))*100)), ")",sep=""),
                  "",
                  paste((table(EUR_data$alcohol_standardized)[2])," (",(round(table(EUR_data$alcohol_standardized)[2]/(nrow(EUR_data))*100)), ")",sep=""),
                  paste((table(EUR_data$alcohol_standardized)[3])," (",(round(table(EUR_data$alcohol_standardized)[3]/(nrow(EUR_data))*100)), ")",sep=""),
                  paste((table(EUR_data$alcohol_standardized)[1])," (",(round(table(EUR_data$alcohol_standardized)[1]/(nrow(EUR_data))*100)), ")",sep=""),
                  "",
                  paste((table(EUR_data$physical_activity_standardized)[2])," (",(round(table(EUR_data$physical_activity_standardized)[2]/(nrow(EUR_data))*100)), ")",sep=""),
                  paste((table(EUR_data$physical_activity_standardized)[3])," (",(round(table(EUR_data$physical_activity_standardized)[3]/(nrow(EUR_data))*100)), ")",sep=""),
                  paste((table(EUR_data$physical_activity_standardized)[1])," (",(round(table(EUR_data$physical_activity_standardized)[1]/(nrow(EUR_data))*100)), ")",sep=""),
                  paste((table(EUR_data$SelfID)[2])," (",(round(table(EUR_data$SelfID)[2]/(nrow(EUR_data))*100)), ")",sep=""),
                  paste((table(EUR_data$Strict)[2])," (",(round(table(EUR_data$Strict)[2]/(nrow(EUR_data))*100)), ")",sep=""),
              paste(summary(EUR_data$Strict)[7]," (",(round(summary(EUR_data$Strict)[7]/(nrow(EUR_data))*100)), ")",sep=""),
                  (paste((nrow(EUR_data[EUR_data$statin_sum_standardized == 1, ]))," (",(round(nrow(EUR_data[EUR_data$statin_sum_standardized == 1, ])/(nrow(EUR_data))*100)), ")",sep="")))

EAS_demog = c("East Asian",
                  paste("Vegetarian"," (n=",(nrow(EAS_data)),")", sep=""),
                  paste(format(round(mean(EAS_data$Age),digit=0))," (",(round(sd(EAS_data$Age),digit=1)), ")",sep=""),
                  paste((nrow(EAS_data[EAS_data$Sex == 0, ]))," (",(round(nrow(EAS_data[EAS_data$Sex == 0, ])/(nrow(EAS_data))*100)), ")",sep=""),
                  paste(format(round(mean(EAS_data[!is.na(EAS_data$BMI), ]$BMI),digit=0))," (",(round(sd(EAS_data[!is.na(EAS_data$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(EAS_data[!is.na(EAS_data$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(EAS_data[!is.na(EAS_data$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_data[!is.na(EAS_data$LDL), ]$LDL),digit=2))," (",(round(sd(EAS_data[!is.na(EAS_data$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_data[!is.na(EAS_data$HDL), ]$HDL),digit=2))," (",(round(sd(EAS_data[!is.na(EAS_data$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(EAS_data[!is.na(EAS_data$TAGs), ]$TAGs),digit=2))," (",(round(sd(EAS_data[!is.na(EAS_data$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(EAS_data$smoking_standardized)[2])," (",(round(table(EAS_data$smoking_standardized)[2]/(nrow(EAS_data))*100)), ")",sep=""),
                  paste((table(EAS_data$smoking_standardized)[3])," (",(round(table(EAS_data$smoking_standardized)[3]/(nrow(EAS_data))*100)), ")",sep=""),
                  paste((table(EAS_data$smoking_standardized)[1])," (",(round(table(EAS_data$smoking_standardized)[1]/(nrow(EAS_data))*100)), ")",sep=""),
                  "",
                  paste((table(EAS_data$alcohol_standardized)[2])," (",(round(table(EAS_data$alcohol_standardized)[2]/(nrow(EAS_data))*100)), ")",sep=""),
                  paste((table(EAS_data$alcohol_standardized)[3])," (",(round(table(EAS_data$alcohol_standardized)[3]/(nrow(EAS_data))*100)), ")",sep=""),
                  paste((table(EAS_data$alcohol_standardized)[1])," (",(round(table(EAS_data$alcohol_standardized)[1]/(nrow(EAS_data))*100)), ")",sep=""),
                  "",
                  paste((table(EAS_data$physical_activity_standardized)[2])," (",(round(table(EAS_data$physical_activity_standardized)[2]/(nrow(EAS_data))*100)), ")",sep=""),
                  paste((table(EAS_data$physical_activity_standardized)[3])," (",(round(table(EAS_data$physical_activity_standardized)[3]/(nrow(EAS_data))*100)), ")",sep=""),
                  paste((table(EAS_data$physical_activity_standardized)[1])," (",(round(table(EAS_data$physical_activity_standardized)[1]/(nrow(EAS_data))*100)), ")",sep=""),
              paste((table(EAS_data$SelfID)[2])," (",(round(table(EAS_data$SelfID)[2]/(nrow(EAS_data))*100)), ")",sep=""),
              paste((table(EAS_data$Strict)[2])," (",(round(table(EAS_data$Strict)[2]/(nrow(EAS_data))*100)), ")",sep=""),
              paste(summary(EAS_data$Strict)[7]," (",(round(summary(EAS_data$Strict)[7]/(nrow(EAS_data))*100)), ")",sep=""),
                  (paste((nrow(EAS_data[EAS_data$statin_sum_standardized == 1, ]))," (",(round(nrow(EAS_data[EAS_data$statin_sum_standardized == 1, ])/(nrow(EAS_data))*100)), ")",sep="")))

CSA_demog = c("Central/South Asian",
                  paste("Vegetarian"," (n=",(nrow(CSA_data)),")", sep=""),
                  paste(format(round(mean(CSA_data$Age),digit=0))," (",(round(sd(CSA_data$Age),digit=1)), ")",sep=""),
                  paste((nrow(CSA_data[CSA_data$Sex == 0, ]))," (",(round(nrow(CSA_data[CSA_data$Sex == 0, ])/(nrow(CSA_data))*100)), ")",sep=""),
                  paste(format(round(mean(CSA_data[!is.na(CSA_data$BMI), ]$BMI),digit=0))," (",(round(sd(CSA_data[!is.na(CSA_data$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(CSA_data[!is.na(CSA_data$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(CSA_data[!is.na(CSA_data$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_data[!is.na(CSA_data$LDL), ]$LDL),digit=2))," (",(round(sd(CSA_data[!is.na(CSA_data$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_data[!is.na(CSA_data$HDL), ]$HDL),digit=2))," (",(round(sd(CSA_data[!is.na(CSA_data$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(CSA_data[!is.na(CSA_data$TAGs), ]$TAGs),digit=2))," (",(round(sd(CSA_data[!is.na(CSA_data$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(CSA_data$smoking_standardized)[2])," (",(round(table(CSA_data$smoking_standardized)[2]/(nrow(CSA_data))*100)), ")",sep=""),
                  paste((table(CSA_data$smoking_standardized)[3])," (",(round(table(CSA_data$smoking_standardized)[3]/(nrow(CSA_data))*100)), ")",sep=""),
                  paste((table(CSA_data$smoking_standardized)[1])," (",(round(table(CSA_data$smoking_standardized)[1]/(nrow(CSA_data))*100)), ")",sep=""),
                  "",
                  paste((table(CSA_data$alcohol_standardized)[2])," (",(round(table(CSA_data$alcohol_standardized)[2]/(nrow(CSA_data))*100)), ")",sep=""),
                  paste((table(CSA_data$alcohol_standardized)[3])," (",(round(table(CSA_data$alcohol_standardized)[3]/(nrow(CSA_data))*100)), ")",sep=""),
                  paste((table(CSA_data$alcohol_standardized)[1])," (",(round(table(CSA_data$alcohol_standardized)[1]/(nrow(CSA_data))*100)), ")",sep=""),
                  "",
                  paste((table(CSA_data$physical_activity_standardized)[2])," (",(round(table(CSA_data$physical_activity_standardized)[2]/(nrow(CSA_data))*100)), ")",sep=""),
                  paste((table(CSA_data$physical_activity_standardized)[3])," (",(round(table(CSA_data$physical_activity_standardized)[3]/(nrow(CSA_data))*100)), ")",sep=""),
                  paste((table(CSA_data$physical_activity_standardized)[1])," (",(round(table(CSA_data$physical_activity_standardized)[1]/(nrow(CSA_data))*100)), ")",sep=""),
              paste((table(CSA_data$SelfID)[2])," (",(round(table(CSA_data$SelfID)[2]/(nrow(CSA_data))*100)), ")",sep=""),
              paste((table(CSA_data$Strict)[2])," (",(round(table(CSA_data$Strict)[2]/(nrow(CSA_data))*100)), ")",sep=""),
              paste(summary(CSA_data$Strict)[7]," (",(round(summary(CSA_data$Strict)[7]/(nrow(CSA_data))*100)), ")",sep=""),
              (paste((nrow(CSA_data[CSA_data$statin_sum_standardized == 1, ]))," (",(round(nrow(CSA_data[CSA_data$statin_sum_standardized == 1, ])/(nrow(CSA_data))*100)), ")",sep="")))

AFR_demog = c("African",
                  paste("Vegetarian"," (n=",(nrow(AFR_data)),")", sep=""),
                  paste(format(round(mean(AFR_data$Age),digit=0))," (",(round(sd(AFR_data$Age),digit=1)), ")",sep=""),
                  paste((nrow(AFR_data[AFR_data$Sex == 0, ]))," (",(round(nrow(AFR_data[AFR_data$Sex == 0, ])/(nrow(AFR_data))*100)), ")",sep=""),
                  paste(format(round(mean(AFR_data[!is.na(AFR_data$BMI), ]$BMI),digit=0))," (",(round(sd(AFR_data[!is.na(AFR_data$BMI), ]$BMI),digit=1)), ")",sep=""),
                  paste(format(round(mean(AFR_data[!is.na(AFR_data$Tot_Chol), ]$Tot_Chol),digit=2))," (",(round(sd(AFR_data[!is.na(AFR_data$Tot_Chol), ]$Tot_Chol),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_data[!is.na(AFR_data$LDL), ]$LDL),digit=2))," (",(round(sd(AFR_data[!is.na(AFR_data$LDL), ]$LDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_data[!is.na(AFR_data$HDL), ]$HDL),digit=2))," (",(round(sd(AFR_data[!is.na(AFR_data$HDL), ]$HDL),digit=2)), ")",sep=""),
                  paste(format(round(mean(AFR_data[!is.na(AFR_data$TAGs), ]$TAGs),digit=2))," (",(round(sd(AFR_data[!is.na(AFR_data$TAGs), ]$TAGs),digit=2)), ")",sep=""),
                  "",
                  paste((table(AFR_data$smoking_standardized)[2])," (",(round(table(AFR_data$smoking_standardized)[2]/(nrow(AFR_data))*100)), ")",sep=""),
                  paste((table(AFR_data$smoking_standardized)[3])," (",(round(table(AFR_data$smoking_standardized)[3]/(nrow(AFR_data))*100)), ")",sep=""),
                  paste((table(AFR_data$smoking_standardized)[1])," (",(round(table(AFR_data$smoking_standardized)[1]/(nrow(AFR_data))*100)), ")",sep=""),
                  "",
                  paste((table(AFR_data$alcohol_standardized)[2])," (",(round(table(AFR_data$alcohol_standardized)[2]/(nrow(AFR_data))*100)), ")",sep=""),
                  paste((table(AFR_data$alcohol_standardized)[3])," (",(round(table(AFR_data$alcohol_standardized)[3]/(nrow(AFR_data))*100)), ")",sep=""),
                  paste((table(AFR_data$alcohol_standardized)[1])," (",(round(table(AFR_data$alcohol_standardized)[1]/(nrow(AFR_data))*100)), ")",sep=""),
                  "",
                  paste((table(AFR_data$physical_activity_standardized)[2])," (",(round(table(AFR_data$physical_activity_standardized)[2]/(nrow(AFR_data))*100)), ")",sep=""),
                  paste((table(AFR_data$physical_activity_standardized)[3])," (",(round(table(AFR_data$physical_activity_standardized)[3]/(nrow(AFR_data))*100)), ")",sep=""),
                  paste((table(AFR_data$physical_activity_standardized)[1])," (",(round(table(AFR_data$physical_activity_standardized)[1]/(nrow(AFR_data))*100)), ")",sep=""),
              paste((table(AFR_data$SelfID)[2])," (",(round(table(AFR_data$SelfID)[2]/(nrow(AFR_data))*100)), ")",sep=""),
              paste((table(AFR_data$Strict)[2])," (",(round(table(AFR_data$Strict)[2]/(nrow(AFR_data))*100)), ")",sep=""),    
              paste(summary(AFR_data$Strict)[7]," (",(round(summary(AFR_data$Strict)[7]/(nrow(AFR_data))*100)), ")",sep=""),
              (paste((nrow(AFR_data[AFR_data$statin_sum_standardized == 1, ]))," (",(round(nrow(AFR_data[AFR_data$statin_sum_standardized == 1, ])/(nrow(AFR_data))*100)), ")",sep="")))

veg_demographics_main_new <- cbind(row_label, EUR_demog, EAS_demog, CSA_demog, AFR_demog)


save(veg_demographics_main_new, file = "veg_demographics_main_new.RData")

write.table(veg_demographics_main_new, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_demographics_main_new.txt", col.names = FALSE, append = TRUE,
            row.names = F, quote = FALSE, na = "",sep='\t')


