# generating tables for plot generation (EUR and CSA interaction plots)

'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr)
library(tidyverse)
library(readr) 
library("lubridate")
library("survival")

Pathway=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/")
Pathway_out=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/")

load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_excluded_data_strict_PRS_factorized.RData")
load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_excluded_data_self_PRS_factorized.RData")

bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS_factorized)
bd_pheno_self <- data.frame(veg_excluded_data_self_PRS_factorized)

bd_pheno_EUR <- bd_pheno_strict[(bd_pheno_strict$pop=="EUR"), ]
bd_pheno_SAS <- bd_pheno_strict[(bd_pheno_strict$pop=="CSA"), ]

bd_pheno_self_EUR <- bd_pheno_self[(bd_pheno_self$pop=="EUR"), ]
bd_pheno_self_SAS <- bd_pheno_self[(bd_pheno_self$pop=="CSA"), ]

#EUR

for (e in c("TC","LDL","HDL","TGs")) {
  for (popul in c("EUR")) {
      if (e == c("TC")) {
        exposure = "Total cholesterol"
        final_res_2=data.frame("Total cholesterol", NA, NA, NA)
      } else if (e == c("LDL")) {
        exposure = "LDL cholesterol"
        final_res_2=data.frame("LDL cholesterol", NA, NA, NA)
      } else if (e == c("HDL")) {
        exposure = "HDL cholesterol"
        final_res_2=data.frame("HDL cholesterol", NA, NA, NA)
      } else if (e == c("TGs")) {
        exposure = "Triglycerides"
        final_res_2=data.frame("Triglycerides", NA, NA, NA)
      }
    write.table(final_res_2, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/new_forest_plot_EUR.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
      
    if (popul == "EUR") {
      bd_pheno_copy <- bd_pheno_EUR
      bd_pheno_tem1=bd_pheno_EUR
      bd_pheno_tem=bd_pheno_EUR
    } else if (popul == "SAS") {
      bd_pheno_copy <- bd_pheno_SAS
      bd_pheno_tem1=bd_pheno_SAS
      bd_pheno_tem=bd_pheno_SAS
    }
    
    bd_pheno_strict=bd_pheno_tem1
    bd_pheno_strict=bd_pheno_strict %>% select(FID,e,everything())
    #Standardize genotypic lipid
    df_tem <- bd_pheno_strict[colnames(bd_pheno_strict)[2]]
    colnames(df_tem)[1]="Lipds"
    bd_pheno_strict=cbind(bd_pheno_strict,df_tem)
    Veg_lipid_strict <- lm(Lipds ~ Strict + Sex + Age + I(Age^2) + Array + 
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
    beta_strict=ifelse(abs(summary(Veg_lipid_strict)$coefficients[2,1])<0.0005,formatC(summary(Veg_lipid_strict)$coefficients[2,1], format = "e", digits = 2),format(round(summary(Veg_lipid_strict)$coefficients[2,1],digit=3), nsmall = 3))
    se_strict=ifelse(summary(Veg_lipid_strict)$coefficients[2,2]<0.0005,formatC(summary(Veg_lipid_strict)$coefficients[2,2], format = "e", digits = 2),format(round(summary(Veg_lipid_strict)$coefficients[2,2],digit=3), nsmall = 3))
    pval_strict=ifelse(summary(Veg_lipid_strict)$coefficients[2,4]<0.0005,formatC(summary(Veg_lipid_strict)$coefficients[2,4], format = "e", digits = 2),format(round(summary(Veg_lipid_strict)$coefficients[2,4],digit=3), nsmall = 3))
    final_res_2=data.frame("    Verified vegetarian",beta_strict,se_strict,pval_strict)
    write.table(final_res_2, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/new_forest_plot_EUR.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    if (popul == "EUR") {
      bd_pheno_copy <- bd_pheno_self_EUR
      bd_pheno_tem1=bd_pheno_self_EUR
      bd_pheno_tem=bd_pheno_self_EUR
    } else if (popul == "SAS") {
      bd_pheno_copy <- bd_pheno_self_SAS
      bd_pheno_tem1=bd_pheno_self_SAS
      bd_pheno_tem=bd_pheno_self_SAS
    }

    bd_pheno_self=bd_pheno_tem1
    bd_pheno_self=bd_pheno_self %>% select(FID,e,everything())
    #Standardize genotypic lipid
    df_tem <- bd_pheno_self[colnames(bd_pheno_self)[2]]
    colnames(df_tem)[1]="Lipds"
    bd_pheno_self=cbind(bd_pheno_self,df_tem)
    Veg_lipid_self <- lm(Lipds ~ SelfID + Sex + Age + I(Age^2) + Array + 
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
    beta_self=ifelse(abs(summary(Veg_lipid_self)$coefficients[2,1])<0.0005,formatC(summary(Veg_lipid_self)$coefficients[2,1], format = "e", digits = 2),format(round(summary(Veg_lipid_self)$coefficients[2,1],digit=3), nsmall = 3))
    se_self=ifelse(summary(Veg_lipid_self)$coefficients[2,2]<0.0005,formatC(summary(Veg_lipid_self)$coefficients[2,2], format = "e", digits = 2),format(round(summary(Veg_lipid_self)$coefficients[2,2],digit=3), nsmall = 3))
    pval_self=ifelse(summary(Veg_lipid_self)$coefficients[2,4]<0.0005,formatC(summary(Veg_lipid_self)$coefficients[2,4], format = "e", digits = 2),format(round(summary(Veg_lipid_self)$coefficients[2,4],digit=3), nsmall = 3))
    final_res_2=data.frame("    Self-reported vegetarian",beta_self,se_self,pval_self)
    write.table(final_res_2, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/new_forest_plot_EUR.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
  }
}






#CSA

for (e in c("TC","LDL","HDL","TGs")) {
  for (popul in c("SAS")) {
    if (e == c("TC")) {
      exposure = "Total cholesterol"
      final_res_2=data.frame("Total cholesterol", NA, NA, NA)
    } else if (e == c("LDL")) {
      exposure = "LDL cholesterol"
      final_res_2=data.frame("LDL cholesterol", NA, NA, NA)
    } else if (e == c("HDL")) {
      exposure = "HDL cholesterol"
      final_res_2=data.frame("HDL cholesterol", NA, NA, NA)
    } else if (e == c("TGs")) {
      exposure = "Triglycerides"
      final_res_2=data.frame("Triglycerides", NA, NA, NA)
    }
    write.table(final_res_2, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/new_forest_plot_CSA.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    if (popul == "EUR") {
      bd_pheno_copy <- bd_pheno_EUR
      bd_pheno_tem1=bd_pheno_EUR
      bd_pheno_tem=bd_pheno_EUR
    } else if (popul == "SAS") {
      bd_pheno_copy <- bd_pheno_SAS
      bd_pheno_tem1=bd_pheno_SAS
      bd_pheno_tem=bd_pheno_SAS
    }
    
    bd_pheno_strict=bd_pheno_tem1
    bd_pheno_strict=bd_pheno_strict %>% select(FID,e,everything())
    #Standardize genotypic lipid
    df_tem <- bd_pheno_strict[colnames(bd_pheno_strict)[2]]
    colnames(df_tem)[1]="Lipds"
    bd_pheno_strict=cbind(bd_pheno_strict,df_tem)
    Veg_lipid_strict <- lm(Lipds ~ Strict + Sex + Age + I(Age^2) + Array + 
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
    beta_strict=ifelse(abs(summary(Veg_lipid_strict)$coefficients[2,1])<0.0005,formatC(summary(Veg_lipid_strict)$coefficients[2,1], format = "e", digits = 2),format(round(summary(Veg_lipid_strict)$coefficients[2,1],digit=3), nsmall = 3))
    se_strict=ifelse(summary(Veg_lipid_strict)$coefficients[2,2]<0.0005,formatC(summary(Veg_lipid_strict)$coefficients[2,2], format = "e", digits = 2),format(round(summary(Veg_lipid_strict)$coefficients[2,2],digit=3), nsmall = 3))
    pval_strict=ifelse(summary(Veg_lipid_strict)$coefficients[2,4]<0.0005,formatC(summary(Veg_lipid_strict)$coefficients[2,4], format = "e", digits = 2),format(round(summary(Veg_lipid_strict)$coefficients[2,4],digit=3), nsmall = 3))
    final_res_2=data.frame("    Verified vegetarian",beta_strict,se_strict,pval_strict)
    write.table(final_res_2, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/new_forest_plot_CSA.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    if (popul == "EUR") {
      bd_pheno_copy <- bd_pheno_self_EUR
      bd_pheno_tem1=bd_pheno_self_EUR
      bd_pheno_tem=bd_pheno_self_EUR
    } else if (popul == "SAS") {
      bd_pheno_copy <- bd_pheno_self_SAS
      bd_pheno_tem1=bd_pheno_self_SAS
      bd_pheno_tem=bd_pheno_self_SAS
    }
    
    bd_pheno_self=bd_pheno_tem1
    bd_pheno_self=bd_pheno_self %>% select(FID,e,everything())
    #Standardize genotypic lipid
    df_tem <- bd_pheno_self[colnames(bd_pheno_self)[2]]
    colnames(df_tem)[1]="Lipds"
    bd_pheno_self=cbind(bd_pheno_self,df_tem)
    Veg_lipid_self <- lm(Lipds ~ SelfID + Sex + Age + I(Age^2) + Array + 
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
    beta_self=ifelse(abs(summary(Veg_lipid_self)$coefficients[2,1])<0.0005,formatC(summary(Veg_lipid_self)$coefficients[2,1], format = "e", digits = 2),format(round(summary(Veg_lipid_self)$coefficients[2,1],digit=3), nsmall = 3))
    se_self=ifelse(summary(Veg_lipid_self)$coefficients[2,2]<0.0005,formatC(summary(Veg_lipid_self)$coefficients[2,2], format = "e", digits = 2),format(round(summary(Veg_lipid_self)$coefficients[2,2],digit=3), nsmall = 3))
    pval_self=ifelse(summary(Veg_lipid_self)$coefficients[2,4]<0.0005,formatC(summary(Veg_lipid_self)$coefficients[2,4], format = "e", digits = 2),format(round(summary(Veg_lipid_self)$coefficients[2,4],digit=3), nsmall = 3))
    final_res_2=data.frame("    Self-reported vegetarian",beta_self,se_self,pval_self)
    write.table(final_res_2, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/new_forest_plot_CSA.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
  }
}
