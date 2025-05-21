# generates tables for making supplemental plots

'%ni%' <- Negate('%in%') 
library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library("lubridate")


Pathway=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files")
Pathway_out=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files")

load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_excluded_data_strict_PRS_factorized.RData")
load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_excluded_data_self_PRS_factorized.RData")

bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS_factorized)
bd_pheno_self <- data.frame(veg_excluded_data_self_PRS_factorized)

################################### ###########       No UKB   Graham SE et al. Nature (2021)
################################### Strict
for (x in c("Graham")) {
  for (popul in c("EUR")) {
    
    for (e in c("TC","LDL","HDL","TGs")) {
      if(e=="TC") {
        Lipid = "Total cholesterol"
      }
      else if (e=="LDL") {
        Lipid = "LDL cholesterol"
      }
      else if (e=="HDL") {
        Lipid = "HDL cholesterol"
      }
      else if (e=="TGs") {
        Lipid = "Triglycerides"
      }
      Lipids_name=paste(e, "_",x,"_no_UKB_",popul,"_PRS", sep="")
      
      bd_pheno_strict=bd_pheno_strict %>% select(FID,e,Lipids_name, everything())
      
      ########### PRS
      df_tem <- bd_pheno_strict[colnames(bd_pheno_strict)[3]]
      colnames(df_tem)[1]="PRS_lipds"
      df_tem$un_sd_PRS_lipds=df_tem$PRS_lipds
      df_tem$PRS_lipds=scale(df_tem$PRS_lipds)
      bd_pheno_strict=cbind(bd_pheno_strict,df_tem)
      
      #Standardize genotypic lipid
      df_tem <- bd_pheno_strict[colnames(bd_pheno_strict)[2]]
      colnames(df_tem)[1]="Lipds"
      bd_pheno_strict=cbind(bd_pheno_strict,df_tem)
      
      TC_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                        ac10003 + ac11001 + ac11002 + 
                        ac11003 + ac11004 + ac11005 + 
                        ac11006 + ac11007 + ac11008 + 
                        ac11009 + ac11010 + ac11011 + 
                        ac11011 + ac11012 + ac11013 + 
                        ac11014 + ac11016 + ac11017 + 
                        ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                        PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                        PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                        Physical_activity_cov, data = filter(bd_pheno_strict, Strict == 1), na.action = "na.omit")
      
      TC_No_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                           ac10003 + ac11001 + ac11002 + 
                           ac11003 + ac11004 + ac11005 + 
                           ac11006 + ac11007 + ac11008 + 
                           ac11009 + ac11010 + ac11011 + 
                           ac11011 + ac11012 + ac11013 + 
                           ac11014 + ac11016 + ac11017 + 
                           ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                           PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                           PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                           Physical_activity_cov, data = filter(bd_pheno_strict, Strict == 0), na.action = "na.omit")
      
      ########### GxE 
      TC_Veg_lm_int <- lm(Lipds ~ I(PRS_lipds*Strict) + Strict + PRS_lipds + Sex + Age + I(Age^2) + Array + 
                            ac10003 + ac11001 + ac11002 + 
                            ac11003 + ac11004 + ac11005 + 
                            ac11006 + ac11007 + ac11008 + 
                            ac11009 + ac11010 + ac11011 + 
                            ac11011 + ac11012 + ac11013 + 
                            ac11014 + ac11016 + ac11017 + 
                            ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                            PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                            PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                            Physical_activity_cov, data = bd_pheno_strict, na.action = "na.omit")
      
      #Veg
      #beta
      a1=summary(TC_Veg_lm)$coefficients[2,1]
      #lower and upper
      b1=summary(TC_Veg_lm)$coefficients[2,2]
      #pval
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))

      #No_Veg
      a=summary(TC_No_Veg_lm)$coefficients[2,1]
      b=summary(TC_No_Veg_lm)$coefficients[2,2]
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      
      #pvalint
      a_int=summary(TC_Veg_lm_int)$coefficients[2,1]
      b_int=summary(TC_Veg_lm_int)$coefficients[2,2]
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      final_res=data.frame(Lipid, NA,NA,NA, pval_int)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_noukb_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Vegetarian", a1,b1,pval1, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_noukb_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Non-vegetarian", a,b,pval, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_noukb_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS_factorized)
      bd_pheno_self <- data.frame(veg_excluded_data_self_PRS_factorized)      
    }
  }
}

################################### Self
for (x in c("Graham")) {
  for (popul in c("EUR")) {
    
    for (e in c("TC","LDL","HDL","TGs")) {
      
      if(e=="TC") {
        Lipid = "Total cholesterol"
      }
      else if (e=="LDL") {
        Lipid = "LDL cholesterol"
      }
      else if (e=="HDL") {
        Lipid = "HDL cholesterol"
      }
      else if (e=="TGs") {
        Lipid = "Triglycerides"
      }
      
      Lipids_name=paste(e, "_",x,"_no_UKB_",popul,"_PRS", sep="")
      
      bd_pheno_self=bd_pheno_self %>% select(FID,e,Lipids_name, everything())
      
      ########### PRS
      df_tem <- bd_pheno_self[colnames(bd_pheno_self)[3]]
      colnames(df_tem)[1]="PRS_lipds"
      df_tem$un_sd_PRS_lipds=df_tem$PRS_lipds
      df_tem$PRS_lipds=scale(df_tem$PRS_lipds)
      bd_pheno_self=cbind(bd_pheno_self,df_tem)
      
      #Standardize genotypic lipid
      df_tem <- bd_pheno_self[colnames(bd_pheno_self)[2]]
      colnames(df_tem)[1]="Lipds"
      bd_pheno_self=cbind(bd_pheno_self,df_tem)
      
      TC_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                        ac10003 + ac11001 + ac11002 + 
                        ac11003 + ac11004 + ac11005 + 
                        ac11006 + ac11007 + ac11008 + 
                        ac11009 + ac11010 + ac11011 + 
                        ac11011 + ac11012 + ac11013 + 
                        ac11014 + ac11016 + ac11017 + 
                        ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                        PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                        PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                        Physical_activity_cov, data = filter(bd_pheno_self, SelfID == 1), na.action = "na.omit")
      
      TC_No_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                           ac10003 + ac11001 + ac11002 + 
                           ac11003 + ac11004 + ac11005 + 
                           ac11006 + ac11007 + ac11008 + 
                           ac11009 + ac11010 + ac11011 + 
                           ac11011 + ac11012 + ac11013 + 
                           ac11014 + ac11016 + ac11017 + 
                           ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                           PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                           PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                           Physical_activity_cov, data = filter(bd_pheno_self, SelfID == 0), na.action = "na.omit")
      
      ########### GxE 
      TC_Veg_lm_int <- lm(Lipds ~ I(PRS_lipds*SelfID) + SelfID + PRS_lipds + Sex + Age + I(Age^2) + Array + 
                            ac10003 + ac11001 + ac11002 + 
                            ac11003 + ac11004 + ac11005 + 
                            ac11006 + ac11007 + ac11008 + 
                            ac11009 + ac11010 + ac11011 + 
                            ac11011 + ac11012 + ac11013 + 
                            ac11014 + ac11016 + ac11017 + 
                            ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                            PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                            PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                            Physical_activity_cov, data = bd_pheno_self, na.action = "na.omit")
      
      #Veg
      #beta
      a1=summary(TC_Veg_lm)$coefficients[2,1]
      #lower and upper
      b1=summary(TC_Veg_lm)$coefficients[2,2]
      #pval
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))

      #No_Veg
      a=summary(TC_No_Veg_lm)$coefficients[2,1]
      b=summary(TC_No_Veg_lm)$coefficients[2,2]
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      
      #pvalint
      a_int=summary(TC_Veg_lm_int)$coefficients[2,1]
      b_int=summary(TC_Veg_lm_int)$coefficients[2,2]
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      
      final_res=data.frame(Lipid, NA,NA,NA, pval_int)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_noukb_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Vegetarian", a1,b1,pval1, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_noukb_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Non-vegetarian", a,b,pval, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_noukb_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS_factorized)
      bd_pheno_self <- data.frame(veg_excluded_data_self_PRS_factorized)          
      
    }
  }
}

################################### ###########       Graham SE et al. Nature (2021)
################################### Strict
for (x in c("Graham")) {
  for (popul in c("EUR")) {
    
    for (e in c("TC","LDL","HDL","TGs")) {
      
      if(e=="TC") {
        Lipid = "Total cholesterol"
      }
      else if (e=="LDL") {
        Lipid = "LDL cholesterol"
      }
      else if (e=="HDL") {
        Lipid = "HDL cholesterol"
      }
      else if (e=="TGs") {
        Lipid = "Triglycerides"
      }
      
      Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
      
      bd_pheno_strict=bd_pheno_strict %>% select(FID,e,Lipids_name, everything())
      
      ########### PRS
      df_tem <- bd_pheno_strict[colnames(bd_pheno_strict)[3]]
      colnames(df_tem)[1]="PRS_lipds"
      df_tem$un_sd_PRS_lipds=df_tem$PRS_lipds
      df_tem$PRS_lipds=scale(df_tem$PRS_lipds)
      bd_pheno_strict=cbind(bd_pheno_strict,df_tem)
      
      #Standardize genotypic lipid
      df_tem <- bd_pheno_strict[colnames(bd_pheno_strict)[2]]
      colnames(df_tem)[1]="Lipds"
      bd_pheno_strict=cbind(bd_pheno_strict,df_tem)
      
      TC_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                        ac10003 + ac11001 + ac11002 + 
                        ac11003 + ac11004 + ac11005 + 
                        ac11006 + ac11007 + ac11008 + 
                        ac11009 + ac11010 + ac11011 + 
                        ac11011 + ac11012 + ac11013 + 
                        ac11014 + ac11016 + ac11017 + 
                        ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                        PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                        PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                        Physical_activity_cov, data = filter(bd_pheno_strict, Strict == 1), na.action = "na.omit")
      
      TC_No_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                           ac10003 + ac11001 + ac11002 + 
                           ac11003 + ac11004 + ac11005 + 
                           ac11006 + ac11007 + ac11008 + 
                           ac11009 + ac11010 + ac11011 + 
                           ac11011 + ac11012 + ac11013 + 
                           ac11014 + ac11016 + ac11017 + 
                           ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                           PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                           PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                           Physical_activity_cov, data = filter(bd_pheno_strict, Strict == 0), na.action = "na.omit")
      
      ########### GxE 
      TC_Veg_lm_int <- lm(Lipds ~ I(PRS_lipds*Strict) + Strict + PRS_lipds + Sex + Age + I(Age^2) + Array + 
                            ac10003 + ac11001 + ac11002 + 
                            ac11003 + ac11004 + ac11005 + 
                            ac11006 + ac11007 + ac11008 + 
                            ac11009 + ac11010 + ac11011 + 
                            ac11011 + ac11012 + ac11013 + 
                            ac11014 + ac11016 + ac11017 + 
                            ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                            PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                            PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                            Physical_activity_cov, data = bd_pheno_strict, na.action = "na.omit")
      
      #Veg
      #beta
      a1=summary(TC_Veg_lm)$coefficients[2,1]
      #lower and upper
      b1=summary(TC_Veg_lm)$coefficients[2,2]
      #pval
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))

      #No_Veg
      a=summary(TC_No_Veg_lm)$coefficients[2,1]
      b=summary(TC_No_Veg_lm)$coefficients[2,2]
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      
      #pvalint
      a_int=summary(TC_Veg_lm_int)$coefficients[2,1]
      b_int=summary(TC_Veg_lm_int)$coefficients[2,2]
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      
      final_res=data.frame(Lipid, NA,NA,NA, pval_int)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_graham_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Vegetarian", a1,b1,pval1, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_graham_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Non-vegetarian", a,b,pval, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_graham_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS_factorized)
      bd_pheno_self <- data.frame(veg_excluded_data_self_PRS_factorized)         
    }
  }
}

################################### Self
for (x in c("Graham")) {
  for (popul in c("EUR")) {
    
    for (e in c("TC","LDL","HDL","TGs")) {
      
      if(e=="TC") {
        Lipid = "Total cholesterol"
      }
      else if (e=="LDL") {
        Lipid = "LDL cholesterol"
      }
      else if (e=="HDL") {
        Lipid = "HDL cholesterol"
      }
      else if (e=="TGs") {
        Lipid = "Triglycerides"
      }
      
      Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
      
      bd_pheno_self=bd_pheno_self %>% select(FID,e,Lipids_name, everything())
      
      ########### PRS
      df_tem <- bd_pheno_self[colnames(bd_pheno_self)[3]]
      colnames(df_tem)[1]="PRS_lipds"
      df_tem$un_sd_PRS_lipds=df_tem$PRS_lipds
      df_tem$PRS_lipds=scale(df_tem$PRS_lipds)
      bd_pheno_self=cbind(bd_pheno_self,df_tem)
      
      #Standardize genotypic lipid
      df_tem <- bd_pheno_self[colnames(bd_pheno_self)[2]]
      colnames(df_tem)[1]="Lipds"
      bd_pheno_self=cbind(bd_pheno_self,df_tem)
      
      TC_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                        ac10003 + ac11001 + ac11002 + 
                        ac11003 + ac11004 + ac11005 + 
                        ac11006 + ac11007 + ac11008 + 
                        ac11009 + ac11010 + ac11011 + 
                        ac11011 + ac11012 + ac11013 + 
                        ac11014 + ac11016 + ac11017 + 
                        ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                        PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                        PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                        Physical_activity_cov, data = filter(bd_pheno_self, SelfID == 1), na.action = "na.omit")
      
      TC_No_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                           ac10003 + ac11001 + ac11002 + 
                           ac11003 + ac11004 + ac11005 + 
                           ac11006 + ac11007 + ac11008 + 
                           ac11009 + ac11010 + ac11011 + 
                           ac11011 + ac11012 + ac11013 + 
                           ac11014 + ac11016 + ac11017 + 
                           ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                           PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                           PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                           Physical_activity_cov, data = filter(bd_pheno_self, SelfID == 0), na.action = "na.omit")
      
      ########### GxE 
      TC_Veg_lm_int <- lm(Lipds ~ I(PRS_lipds*SelfID) + SelfID + PRS_lipds + Sex + Age + I(Age^2) + Array + 
                            ac10003 + ac11001 + ac11002 + 
                            ac11003 + ac11004 + ac11005 + 
                            ac11006 + ac11007 + ac11008 + 
                            ac11009 + ac11010 + ac11011 + 
                            ac11011 + ac11012 + ac11013 + 
                            ac11014 + ac11016 + ac11017 + 
                            ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                            PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                            PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                            Physical_activity_cov, data = bd_pheno_self, na.action = "na.omit")
      
      #Veg
      #beta
      a1=summary(TC_Veg_lm)$coefficients[2,1]
      #lower and upper
      b1=summary(TC_Veg_lm)$coefficients[2,2]
      #pval
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))

      #No_Veg
      a=summary(TC_No_Veg_lm)$coefficients[2,1]
      b=summary(TC_No_Veg_lm)$coefficients[2,2]
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      
      #pvalint
      a_int=summary(TC_Veg_lm_int)$coefficients[2,1]
      b_int=summary(TC_Veg_lm_int)$coefficients[2,2]
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      
      final_res=data.frame(Lipid, NA,NA,NA, pval_int)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_graham_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Vegetarian", a1,b1,pval1, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_graham_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Non-vegetarian", a,b,pval, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_graham_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS_factorized)
      bd_pheno_self <- data.frame(veg_excluded_data_self_PRS_factorized)      
      
    }
  }
}

################################### ###########       Willer CJ et al. Nat. Genet. 2013
################################### Strict
for (x in c("Willer")) {
  for (popul in c("EUR")) {
    
    for (e in c("TC","LDL","HDL","TGs")) {
      
      if(e=="TC") {
        Lipid = "Total cholesterol"
      }
      else if (e=="LDL") {
        Lipid = "LDL cholesterol"
      }
      else if (e=="HDL") {
        Lipid = "HDL cholesterol"
      }
      else if (e=="TGs") {
        Lipid = "Triglycerides"
      }
      
      Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
      
      bd_pheno_strict=bd_pheno_strict %>% select(FID,e,Lipids_name, everything())
      
      ########### PRS
      df_tem <- bd_pheno_strict[colnames(bd_pheno_strict)[3]]
      colnames(df_tem)[1]="PRS_lipds"
      df_tem$un_sd_PRS_lipds=df_tem$PRS_lipds
      df_tem$PRS_lipds=scale(df_tem$PRS_lipds)
      bd_pheno_strict=cbind(bd_pheno_strict,df_tem)
      
      #Standardize genotypic lipid
      df_tem <- bd_pheno_strict[colnames(bd_pheno_strict)[2]]
      colnames(df_tem)[1]="Lipds"
      bd_pheno_strict=cbind(bd_pheno_strict,df_tem)
      
      TC_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                        ac10003 + ac11001 + ac11002 + 
                        ac11003 + ac11004 + ac11005 + 
                        ac11006 + ac11007 + ac11008 + 
                        ac11009 + ac11010 + ac11011 + 
                        ac11011 + ac11012 + ac11013 + 
                        ac11014 + ac11016 + ac11017 + 
                        ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                        PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                        PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                        Physical_activity_cov, data = filter(bd_pheno_strict, Strict == 1), na.action = "na.omit")
      
      TC_No_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                           ac10003 + ac11001 + ac11002 + 
                           ac11003 + ac11004 + ac11005 + 
                           ac11006 + ac11007 + ac11008 + 
                           ac11009 + ac11010 + ac11011 + 
                           ac11011 + ac11012 + ac11013 + 
                           ac11014 + ac11016 + ac11017 + 
                           ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                           PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                           PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                           Physical_activity_cov, data = filter(bd_pheno_strict, Strict == 0), na.action = "na.omit")
      
      ########### GxE 
      TC_Veg_lm_int <- lm(Lipds ~ I(PRS_lipds*Strict) + Strict + PRS_lipds + Sex + Age + I(Age^2) + Array + 
                            ac10003 + ac11001 + ac11002 + 
                            ac11003 + ac11004 + ac11005 + 
                            ac11006 + ac11007 + ac11008 + 
                            ac11009 + ac11010 + ac11011 + 
                            ac11011 + ac11012 + ac11013 + 
                            ac11014 + ac11016 + ac11017 + 
                            ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                            PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                            PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                            Physical_activity_cov, data = bd_pheno_strict, na.action = "na.omit")
      
      #Veg
      #beta
      a1=summary(TC_Veg_lm)$coefficients[2,1]
      #lower and upper
      b1=summary(TC_Veg_lm)$coefficients[2,2]
      #pval
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))

      #No_Veg
      a=summary(TC_No_Veg_lm)$coefficients[2,1]
      b=summary(TC_No_Veg_lm)$coefficients[2,2]
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      
      #pvalint
      a_int=summary(TC_Veg_lm_int)$coefficients[2,1]
      b_int=summary(TC_Veg_lm_int)$coefficients[2,2]
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      
      final_res=data.frame(Lipid, NA,NA,NA, pval_int)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_willer_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Vegetarian", a1,b1,pval1, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_willer_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Non-vegetarian", a,b,pval, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_willer_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS_factorized)
      bd_pheno_self <- data.frame(veg_excluded_data_self_PRS_factorized)      
    }
  }
}

################################### Self
for (x in c("Willer")) {
  for (popul in c("EUR")) {
    
    for (e in c("TC","LDL","HDL","TGs")) {
      
      if(e=="TC") {
        Lipid = "Total cholesterol"
      }
      else if (e=="LDL") {
        Lipid = "LDL cholesterol"
      }
      else if (e=="HDL") {
        Lipid = "HDL cholesterol"
      }
      else if (e=="TGs") {
        Lipid = "Triglycerides"
      }
      
      Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
      
      bd_pheno_self=bd_pheno_self %>% select(FID,e,Lipids_name, everything())
      
      ########### PRS
      df_tem <- bd_pheno_self[colnames(bd_pheno_self)[3]]
      colnames(df_tem)[1]="PRS_lipds"
      df_tem$un_sd_PRS_lipds=df_tem$PRS_lipds
      df_tem$PRS_lipds=scale(df_tem$PRS_lipds)
      bd_pheno_self=cbind(bd_pheno_self,df_tem)
      
      #Standardize genotypic lipid
      df_tem <- bd_pheno_self[colnames(bd_pheno_self)[2]]
      colnames(df_tem)[1]="Lipds"
      bd_pheno_self=cbind(bd_pheno_self,df_tem)
      
      TC_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                        ac10003 + ac11001 + ac11002 + 
                        ac11003 + ac11004 + ac11005 + 
                        ac11006 + ac11007 + ac11008 + 
                        ac11009 + ac11010 + ac11011 + 
                        ac11011 + ac11012 + ac11013 + 
                        ac11014 + ac11016 + ac11017 + 
                        ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                        PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                        PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                        Physical_activity_cov, data = filter(bd_pheno_self, SelfID == 1), na.action = "na.omit")
      
      TC_No_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                           ac10003 + ac11001 + ac11002 + 
                           ac11003 + ac11004 + ac11005 + 
                           ac11006 + ac11007 + ac11008 + 
                           ac11009 + ac11010 + ac11011 + 
                           ac11011 + ac11012 + ac11013 + 
                           ac11014 + ac11016 + ac11017 + 
                           ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                           PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                           PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                           Physical_activity_cov, data = filter(bd_pheno_self, SelfID == 0), na.action = "na.omit")
      
      ########### GxE 
      TC_Veg_lm_int <- lm(Lipds ~ I(PRS_lipds*SelfID) + SelfID + PRS_lipds + Sex + Age + I(Age^2) + Array + 
                            ac10003 + ac11001 + ac11002 + 
                            ac11003 + ac11004 + ac11005 + 
                            ac11006 + ac11007 + ac11008 + 
                            ac11009 + ac11010 + ac11011 + 
                            ac11011 + ac11012 + ac11013 + 
                            ac11014 + ac11016 + ac11017 + 
                            ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                            PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                            PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov + 
                            Physical_activity_cov, data = bd_pheno_self, na.action = "na.omit")
      
      #Veg
      #beta
      a1=summary(TC_Veg_lm)$coefficients[2,1]
      #lower and upper
      b1=summary(TC_Veg_lm)$coefficients[2,2]
      #pval
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))

      #No_Veg
      a=summary(TC_No_Veg_lm)$coefficients[2,1]
      b=summary(TC_No_Veg_lm)$coefficients[2,2]
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      
      #pvalint
      a_int=summary(TC_Veg_lm_int)$coefficients[2,1]
      b_int=summary(TC_Veg_lm_int)$coefficients[2,2]
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      
      final_res=data.frame(Lipid, NA,NA,NA, pval_int)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_willer_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Vegetarian", a1,b1,pval1, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_willer_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      final_res=data.frame("    Non-vegetarian", a,b,pval, NA)
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_willer_simple.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS_factorized)
      bd_pheno_self <- data.frame(veg_excluded_data_self_PRS_factorized)      
      
    }
  }
}






