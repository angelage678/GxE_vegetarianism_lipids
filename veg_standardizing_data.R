# Standardizing data

# Cleaning up PRS information
load(paste("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/","veg_excluded_data.RData",sep = ""))
load(paste("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/","veg_excluded_data_strict.RData",sep = ""))
load(paste("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/","veg_excluded_data_self.RData",sep = ""))

########### PRS
veg_excluded_data_PRS <- veg_excluded_data[is.na(veg_excluded_data$TC_Graham_EUR_PRS)==F | is.na(veg_excluded_data$TC_Graham_AFR_PRS)==F |
                           is.na(veg_excluded_data$TC_Graham_EAS_PRS)==F | is.na(veg_excluded_data$TC_Graham_SAS_PRS)==F | 
                           is.na(veg_excluded_data$TC_Willer_EUR_PRS )==F | 
                           is.na(veg_excluded_data$TC_Graham_no_UKB_EUR_PRS)==F | is.na(veg_excluded_data$TC_Graham_no_UKB_AFR_PRS)==F |
                           is.na(veg_excluded_data$TC_Graham_no_UKB_SAS_PRS)==F ,]
veg_excluded_data_strict_PRS <- veg_excluded_data_strict[is.na(veg_excluded_data_strict$TC_Graham_EUR_PRS)==F | is.na(veg_excluded_data_strict$TC_Graham_AFR_PRS)==F |
                                                           is.na(veg_excluded_data_strict$TC_Graham_EAS_PRS)==F | is.na(veg_excluded_data_strict$TC_Graham_SAS_PRS)==F | 
                                                           is.na(veg_excluded_data_strict$TC_Willer_EUR_PRS )==F | 
                                                           is.na(veg_excluded_data_strict$TC_Graham_no_UKB_EUR_PRS)==F | is.na(veg_excluded_data_strict$TC_Graham_no_UKB_AFR_PRS)==F |
                                                           is.na(veg_excluded_data_strict$TC_Graham_no_UKB_SAS_PRS)==F ,]
veg_excluded_data_self_PRS <- veg_excluded_data_self[is.na(veg_excluded_data_self$TC_Graham_EUR_PRS)==F | is.na(veg_excluded_data_self$TC_Graham_AFR_PRS)==F |
                                                       is.na(veg_excluded_data_self$TC_Graham_EAS_PRS)==F | is.na(veg_excluded_data_self$TC_Graham_SAS_PRS)==F | 
                                                       is.na(veg_excluded_data_self$TC_Willer_EUR_PRS )==F | 
                                                       is.na(veg_excluded_data_self$TC_Graham_no_UKB_EUR_PRS)==F | is.na(veg_excluded_data_self$TC_Graham_no_UKB_AFR_PRS)==F |
                                                       is.na(veg_excluded_data_self$TC_Graham_no_UKB_SAS_PRS)==F ,]

#Standardize blood lipid phenotypes
veg_excluded_data_PRS$TC <- scale(veg_excluded_data_PRS$Tot_Chol)
veg_excluded_data_PRS$LDL <- scale(veg_excluded_data_PRS$LDL)
veg_excluded_data_PRS$HDL <- scale(veg_excluded_data_PRS$HDL)
veg_excluded_data_PRS$TGs <- scale(veg_excluded_data_PRS$TAGs)

veg_excluded_data_strict_PRS$TC <- scale(veg_excluded_data_strict_PRS$Tot_Chol)
veg_excluded_data_strict_PRS$LDL <- scale(veg_excluded_data_strict_PRS$LDL)
veg_excluded_data_strict_PRS$HDL <- scale(veg_excluded_data_strict_PRS$HDL)
veg_excluded_data_strict_PRS$TGs <- scale(veg_excluded_data_strict_PRS$TAGs)

veg_excluded_data_self_PRS$TC <- scale(veg_excluded_data_self_PRS$Tot_Chol)
veg_excluded_data_self_PRS$LDL <- scale(veg_excluded_data_self_PRS$LDL)
veg_excluded_data_self_PRS$HDL <- scale(veg_excluded_data_self_PRS$HDL)
veg_excluded_data_self_PRS$TGs <- scale(veg_excluded_data_self_PRS$TAGs)

#Standardize all PRS
veg_excluded_data_PRS$TC_Willer_EUR_PRS <- scale(veg_excluded_data_PRS$TC_Willer_EUR_PRS) 
veg_excluded_data_PRS$LDL_Willer_EUR_PRS <- scale(veg_excluded_data_PRS$LDL_Willer_EUR_PRS) 
veg_excluded_data_PRS$HDL_Willer_EUR_PRS <- scale(veg_excluded_data_PRS$HDL_Willer_EUR_PRS) 
veg_excluded_data_PRS$TGs_Willer_EUR_PRS <- scale(veg_excluded_data_PRS$TGs_Willer_EUR_PRS) 
veg_excluded_data_PRS$TC_Graham_EUR_PRS <- scale(veg_excluded_data_PRS$TC_Graham_EUR_PRS) 
veg_excluded_data_PRS$LDL_Graham_EUR_PRS <- scale(veg_excluded_data_PRS$LDL_Graham_EUR_PRS) 
veg_excluded_data_PRS$HDL_Graham_EUR_PRS <- scale(veg_excluded_data_PRS$HDL_Graham_EUR_PRS) 
veg_excluded_data_PRS$TGs_Graham_EUR_PRS <- scale(veg_excluded_data_PRS$TGs_Graham_EUR_PRS) 
veg_excluded_data_PRS$TC_Graham_AFR_PRS <- scale(veg_excluded_data_PRS$TC_Graham_AFR_PRS) 
veg_excluded_data_PRS$LDL_Graham_AFR_PRS <- scale(veg_excluded_data_PRS$LDL_Graham_AFR_PRS) 
veg_excluded_data_PRS$HDL_Graham_AFR_PRS <- scale(veg_excluded_data_PRS$HDL_Graham_AFR_PRS) 
veg_excluded_data_PRS$TGs_Graham_AFR_PRS <- scale(veg_excluded_data_PRS$TGs_Graham_AFR_PRS) 
veg_excluded_data_PRS$TC_Graham_EAS_PRS <- scale(veg_excluded_data_PRS$TC_Graham_EAS_PRS) 
veg_excluded_data_PRS$LDL_Graham_EAS_PRS <- scale(veg_excluded_data_PRS$LDL_Graham_EAS_PRS) 
veg_excluded_data_PRS$HDL_Graham_EAS_PRS <- scale(veg_excluded_data_PRS$HDL_Graham_EAS_PRS) 
veg_excluded_data_PRS$TGs_Graham_EAS_PRS <- scale(veg_excluded_data_PRS$TGs_Graham_EAS_PRS) 
veg_excluded_data_PRS$TC_Graham_SAS_PRS <- scale(veg_excluded_data_PRS$TC_Graham_SAS_PRS) 
veg_excluded_data_PRS$LDL_Graham_SAS_PRS <- scale(veg_excluded_data_PRS$LDL_Graham_SAS_PRS) 
veg_excluded_data_PRS$HDL_Graham_SAS_PRS <- scale(veg_excluded_data_PRS$HDL_Graham_SAS_PRS) 
veg_excluded_data_PRS$TGs_Graham_SAS_PRS <- scale(veg_excluded_data_PRS$TGs_Graham_SAS_PRS) 
veg_excluded_data_PRS$TC_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_PRS$TC_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_PRS$LDL_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_PRS$LDL_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_PRS$HDL_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_PRS$HDL_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_PRS$TGs_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_PRS$TGs_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_PRS$TC_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_PRS$TC_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_PRS$LDL_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_PRS$LDL_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_PRS$HDL_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_PRS$HDL_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_PRS$TGs_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_PRS$TGs_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_PRS$TC_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_PRS$TC_Graham_no_UKB_SAS_PRS) 
veg_excluded_data_PRS$LDL_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_PRS$LDL_Graham_no_UKB_SAS_PRS) 
veg_excluded_data_PRS$HDL_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_PRS$HDL_Graham_no_UKB_SAS_PRS) 
veg_excluded_data_PRS$TGs_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_PRS$TGs_Graham_no_UKB_SAS_PRS) 

veg_excluded_data_strict_PRS$TC_Willer_EUR_PRS <- scale(veg_excluded_data_strict_PRS$TC_Willer_EUR_PRS) 
veg_excluded_data_strict_PRS$LDL_Willer_EUR_PRS <- scale(veg_excluded_data_strict_PRS$LDL_Willer_EUR_PRS) 
veg_excluded_data_strict_PRS$HDL_Willer_EUR_PRS <- scale(veg_excluded_data_strict_PRS$HDL_Willer_EUR_PRS) 
veg_excluded_data_strict_PRS$TGs_Willer_EUR_PRS <- scale(veg_excluded_data_strict_PRS$TGs_Willer_EUR_PRS) 
veg_excluded_data_strict_PRS$TC_Graham_EUR_PRS <- scale(veg_excluded_data_strict_PRS$TC_Graham_EUR_PRS) 
veg_excluded_data_strict_PRS$LDL_Graham_EUR_PRS <- scale(veg_excluded_data_strict_PRS$LDL_Graham_EUR_PRS) 
veg_excluded_data_strict_PRS$HDL_Graham_EUR_PRS <- scale(veg_excluded_data_strict_PRS$HDL_Graham_EUR_PRS) 
veg_excluded_data_strict_PRS$TGs_Graham_EUR_PRS <- scale(veg_excluded_data_strict_PRS$TGs_Graham_EUR_PRS) 
veg_excluded_data_strict_PRS$TC_Graham_AFR_PRS <- scale(veg_excluded_data_strict_PRS$TC_Graham_AFR_PRS) 
veg_excluded_data_strict_PRS$LDL_Graham_AFR_PRS <- scale(veg_excluded_data_strict_PRS$LDL_Graham_AFR_PRS) 
veg_excluded_data_strict_PRS$HDL_Graham_AFR_PRS <- scale(veg_excluded_data_strict_PRS$HDL_Graham_AFR_PRS) 
veg_excluded_data_strict_PRS$TGs_Graham_AFR_PRS <- scale(veg_excluded_data_strict_PRS$TGs_Graham_AFR_PRS) 
veg_excluded_data_strict_PRS$TC_Graham_EAS_PRS <- scale(veg_excluded_data_strict_PRS$TC_Graham_EAS_PRS) 
veg_excluded_data_strict_PRS$LDL_Graham_EAS_PRS <- scale(veg_excluded_data_strict_PRS$LDL_Graham_EAS_PRS) 
veg_excluded_data_strict_PRS$HDL_Graham_EAS_PRS <- scale(veg_excluded_data_strict_PRS$HDL_Graham_EAS_PRS) 
veg_excluded_data_strict_PRS$TGs_Graham_EAS_PRS <- scale(veg_excluded_data_strict_PRS$TGs_Graham_EAS_PRS) 
veg_excluded_data_strict_PRS$TC_Graham_SAS_PRS <- scale(veg_excluded_data_strict_PRS$TC_Graham_SAS_PRS) 
veg_excluded_data_strict_PRS$LDL_Graham_SAS_PRS <- scale(veg_excluded_data_strict_PRS$LDL_Graham_SAS_PRS) 
veg_excluded_data_strict_PRS$HDL_Graham_SAS_PRS <- scale(veg_excluded_data_strict_PRS$HDL_Graham_SAS_PRS) 
veg_excluded_data_strict_PRS$TGs_Graham_SAS_PRS <- scale(veg_excluded_data_strict_PRS$TGs_Graham_SAS_PRS) 
veg_excluded_data_strict_PRS$TC_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_strict_PRS$TC_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_strict_PRS$LDL_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_strict_PRS$LDL_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_strict_PRS$HDL_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_strict_PRS$HDL_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_strict_PRS$TGs_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_strict_PRS$TGs_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_strict_PRS$TC_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_strict_PRS$TC_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_strict_PRS$LDL_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_strict_PRS$LDL_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_strict_PRS$HDL_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_strict_PRS$HDL_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_strict_PRS$TGs_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_strict_PRS$TGs_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_strict_PRS$TC_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_strict_PRS$TC_Graham_no_UKB_SAS_PRS) 
veg_excluded_data_strict_PRS$LDL_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_strict_PRS$LDL_Graham_no_UKB_SAS_PRS) 
veg_excluded_data_strict_PRS$HDL_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_strict_PRS$HDL_Graham_no_UKB_SAS_PRS) 
veg_excluded_data_strict_PRS$TGs_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_strict_PRS$TGs_Graham_no_UKB_SAS_PRS) 

veg_excluded_data_self_PRS$TC_Willer_EUR_PRS <- scale(veg_excluded_data_self_PRS$TC_Willer_EUR_PRS) 
veg_excluded_data_self_PRS$LDL_Willer_EUR_PRS <- scale(veg_excluded_data_self_PRS$LDL_Willer_EUR_PRS) 
veg_excluded_data_self_PRS$HDL_Willer_EUR_PRS <- scale(veg_excluded_data_self_PRS$HDL_Willer_EUR_PRS) 
veg_excluded_data_self_PRS$TGs_Willer_EUR_PRS <- scale(veg_excluded_data_self_PRS$TGs_Willer_EUR_PRS) 
veg_excluded_data_self_PRS$TC_Graham_EUR_PRS <- scale(veg_excluded_data_self_PRS$TC_Graham_EUR_PRS) 
veg_excluded_data_self_PRS$LDL_Graham_EUR_PRS <- scale(veg_excluded_data_self_PRS$LDL_Graham_EUR_PRS) 
veg_excluded_data_self_PRS$HDL_Graham_EUR_PRS <- scale(veg_excluded_data_self_PRS$HDL_Graham_EUR_PRS) 
veg_excluded_data_self_PRS$TGs_Graham_EUR_PRS <- scale(veg_excluded_data_self_PRS$TGs_Graham_EUR_PRS) 
veg_excluded_data_self_PRS$TC_Graham_AFR_PRS <- scale(veg_excluded_data_self_PRS$TC_Graham_AFR_PRS) 
veg_excluded_data_self_PRS$LDL_Graham_AFR_PRS <- scale(veg_excluded_data_self_PRS$LDL_Graham_AFR_PRS) 
veg_excluded_data_self_PRS$HDL_Graham_AFR_PRS <- scale(veg_excluded_data_self_PRS$HDL_Graham_AFR_PRS) 
veg_excluded_data_self_PRS$TGs_Graham_AFR_PRS <- scale(veg_excluded_data_self_PRS$TGs_Graham_AFR_PRS) 
veg_excluded_data_self_PRS$TC_Graham_EAS_PRS <- scale(veg_excluded_data_self_PRS$TC_Graham_EAS_PRS) 
veg_excluded_data_self_PRS$LDL_Graham_EAS_PRS <- scale(veg_excluded_data_self_PRS$LDL_Graham_EAS_PRS) 
veg_excluded_data_self_PRS$HDL_Graham_EAS_PRS <- scale(veg_excluded_data_self_PRS$HDL_Graham_EAS_PRS) 
veg_excluded_data_self_PRS$TGs_Graham_EAS_PRS <- scale(veg_excluded_data_self_PRS$TGs_Graham_EAS_PRS) 
veg_excluded_data_self_PRS$TC_Graham_SAS_PRS <- scale(veg_excluded_data_self_PRS$TC_Graham_SAS_PRS) 
veg_excluded_data_self_PRS$LDL_Graham_SAS_PRS <- scale(veg_excluded_data_self_PRS$LDL_Graham_SAS_PRS) 
veg_excluded_data_self_PRS$HDL_Graham_SAS_PRS <- scale(veg_excluded_data_self_PRS$HDL_Graham_SAS_PRS) 
veg_excluded_data_self_PRS$TGs_Graham_SAS_PRS <- scale(veg_excluded_data_self_PRS$TGs_Graham_SAS_PRS) 
veg_excluded_data_self_PRS$TC_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_self_PRS$TC_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_self_PRS$LDL_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_self_PRS$LDL_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_self_PRS$HDL_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_self_PRS$HDL_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_self_PRS$TGs_Graham_no_UKB_EUR_PRS <- scale(veg_excluded_data_self_PRS$TGs_Graham_no_UKB_EUR_PRS) 
veg_excluded_data_self_PRS$TC_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_self_PRS$TC_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_self_PRS$LDL_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_self_PRS$LDL_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_self_PRS$HDL_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_self_PRS$HDL_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_self_PRS$TGs_Graham_no_UKB_AFR_PRS <- scale(veg_excluded_data_self_PRS$TGs_Graham_no_UKB_AFR_PRS) 
veg_excluded_data_self_PRS$TC_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_self_PRS$TC_Graham_no_UKB_SAS_PRS) 
veg_excluded_data_self_PRS$LDL_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_self_PRS$LDL_Graham_no_UKB_SAS_PRS) 
veg_excluded_data_self_PRS$HDL_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_self_PRS$HDL_Graham_no_UKB_SAS_PRS) 
veg_excluded_data_self_PRS$TGs_Graham_no_UKB_SAS_PRS <- scale(veg_excluded_data_self_PRS$TGs_Graham_no_UKB_SAS_PRS) 

#save data & new factorized variable
save(veg_excluded_data_PRS, file = "veg_excluded_data_PRS.RData")
veg_excluded_data_PRS_factorized <- veg_excluded_data_PRS

save(veg_excluded_data_strict_PRS, file = "veg_excluded_data_strict_PRS.RData")
veg_excluded_data_strict_PRS_factorized <- veg_excluded_data_strict_PRS

save(veg_excluded_data_self_PRS, file = "veg_excluded_data_self_PRS.RData")
veg_excluded_data_self_PRS_factorized <- veg_excluded_data_self_PRS

#factorize covariates
veg_excluded_data_PRS_factorized$Smoking_status_cov <- as.factor(veg_excluded_data_PRS_factorized$smoking_standardized)
veg_excluded_data_PRS_factorized$Smoking_status_cov<- factor(veg_excluded_data_PRS_factorized$Smoking_status_cov, levels = c("Never","Previous","Current"))
veg_excluded_data_PRS_factorized$Alcohol_status_cov <- as.factor(veg_excluded_data_PRS_factorized$alcohol_standardized)
veg_excluded_data_PRS_factorized$Alcohol_status_cov<- factor(veg_excluded_data_PRS_factorized$Alcohol_status_cov, levels = c("Never","Previous","Current"))
veg_excluded_data_PRS_factorized$Physical_activity_cov <- as.factor(veg_excluded_data_PRS_factorized$physical_activity_standardized)
veg_excluded_data_PRS_factorized$Physical_activity_cov<- factor(veg_excluded_data_PRS_factorized$Physical_activity_cov, levels = c("Low","Moderate","High"))
save(veg_excluded_data_PRS_factorized, file = "veg_excluded_data_PRS_factorized.RData")

veg_excluded_data_strict_PRS_factorized$Smoking_status_cov <- as.factor(veg_excluded_data_strict_PRS_factorized$smoking_standardized)
veg_excluded_data_strict_PRS_factorized$Smoking_status_cov<- factor(veg_excluded_data_strict_PRS_factorized$Smoking_status_cov, levels = c("Never","Previous","Current"))
veg_excluded_data_strict_PRS_factorized$Alcohol_status_cov <- as.factor(veg_excluded_data_strict_PRS_factorized$alcohol_standardized)
veg_excluded_data_strict_PRS_factorized$Alcohol_status_cov<- factor(veg_excluded_data_strict_PRS_factorized$Alcohol_status_cov, levels = c("Never","Previous","Current"))
veg_excluded_data_strict_PRS_factorized$Physical_activity_cov <- as.factor(veg_excluded_data_strict_PRS_factorized$physical_activity_standardized)
veg_excluded_data_strict_PRS_factorized$Physical_activity_cov<- factor(veg_excluded_data_strict_PRS_factorized$Physical_activity_cov, levels = c("Low","Moderate","High"))
save(veg_excluded_data_strict_PRS_factorized, file = "veg_excluded_data_strict_PRS_factorized.RData")

veg_excluded_data_self_PRS_factorized$Smoking_status_cov <- as.factor(veg_excluded_data_self_PRS_factorized$smoking_standardized)
veg_excluded_data_self_PRS_factorized$Smoking_status_cov<- factor(veg_excluded_data_self_PRS_factorized$Smoking_status_cov, levels = c("Never","Previous","Current"))
veg_excluded_data_self_PRS_factorized$Alcohol_status_cov <- as.factor(veg_excluded_data_self_PRS_factorized$alcohol_standardized)
veg_excluded_data_self_PRS_factorized$Alcohol_status_cov<- factor(veg_excluded_data_self_PRS_factorized$Alcohol_status_cov, levels = c("Never","Previous","Current"))
veg_excluded_data_self_PRS_factorized$Physical_activity_cov <- as.factor(veg_excluded_data_self_PRS_factorized$physical_activity_standardized)
veg_excluded_data_self_PRS_factorized$Physical_activity_cov<- factor(veg_excluded_data_self_PRS_factorized$Physical_activity_cov, levels = c("Low","Moderate","High"))
save(veg_excluded_data_self_PRS_factorized, file = "veg_excluded_data_self_PRS_factorized.RData")



