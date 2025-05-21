# veg-lipid effect and interaction effect analysis

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

################################### Graham SE et al. Nature (2021)
################################### Strict
for (x in c("Graham")) {
  final_res=data.frame(NA,NA,"Vegetarian",NA,NA,"Not vegetarian",NA,NA)
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  final_res=data.frame("Lipids","n","β (95% CI)","P-value",
                       "n","β (95% CI)","P-value","Pinteractionb")
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  final_res=data.frame("Graham SE et al. (included UKB data)",NA,NA,NA,NA,NA,NA,NA)
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  for (popul in c("EUR", "AFR", "EAS", "SAS")) {
    
    if(popul == "SAS"){
      final_res=data.frame("CSA",NA,NA,NA,NA,NA,NA,NA)
    }
    else{
      final_res=data.frame(popul,NA,NA,NA,NA,NA,NA,NA)
    }
    
    write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    for (e in c("TC","LDL","HDL","TGs")) {
      
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
      
      if(popul == "EAS"){
        TC_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                          ac10003 + ac11001 + ac11002 + 
                          ac11003 + ac11004 + ac11005 + 
                          ac11006 + ac11007 + ac11008 + 
                          ac11009 + ac11010 + ac11011 + 
                          ac11011 + ac11012 + ac11013 + 
                          ac11014 + ac11016 + ac11017 + 
                          ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                          PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                          PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Alcohol_status_cov +
                          Physical_activity_cov, data = filter(bd_pheno_strict, Strict == 1), na.action = "na.omit")
      }
      else{
        TC_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                          ac10003 + ac11001 + ac11002 + 
                          ac11003 + ac11004 + ac11005 + 
                          ac11006 + ac11007 + ac11008 + 
                          ac11009 + ac11010 + ac11011 + 
                          ac11011 + ac11012 + ac11013 + 
                          ac11014 + ac11016 + ac11017 + 
                          ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                          PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                          PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Alcohol_status_cov + Smoking_status_cov +
                          Physical_activity_cov, data = filter(bd_pheno_strict, Strict == 1), na.action = "na.omit")
      }
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
      
      print(sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      (summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])/sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      print(pnorm(abs(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]) / sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2), lower.tail=FALSE) * 2)
      
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      
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
    
      a1=format(round(summary(TC_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new1=paste(a1," (",b1,", ",c1,")",sep="")
      
      a=format(round(summary(TC_No_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new=paste(a," (",b,", ",c,")",sep="")
      
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      if (e == c("TC")) {
        final_res=data.frame("Total cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("LDL")) {
        final_res=data.frame("LDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("HDL")) {
        final_res=data.frame("HDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("TGs")) {
        final_res=data.frame("Triglycerides, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      }
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_strict <- subset (bd_pheno_strict, select = -c(PRS_lipds,Lipds,un_sd_PRS_lipds))
      
    }
  }
}

################################### Self
for (x in c("Graham")) {
  final_res=data.frame(NA,NA,"Vegetarian",NA,NA,"Not vegetarian",NA,NA)
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  final_res=data.frame("Lipids","n","β (95% CI)","P-value",
                       "n","β (95% CI)","P-value","Pinteractionb")
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  final_res=data.frame("Graham SE et al. (included UKB data)",NA,NA,NA,NA,NA,NA,NA)
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  for (popul in c("EUR", "AFR", "EAS", "SAS")) {
    
    if(popul == "SAS"){
      final_res=data.frame("CSA",NA,NA,NA,NA,NA,NA,NA)
    }
    else{
      final_res=data.frame(popul,NA,NA,NA,NA,NA,NA,NA)
    }
    
    write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    for (e in c("TC","LDL","HDL","TGs")) {
      
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
      
      if(popul == "EAS"){
        TC_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                          ac10003 + ac11001 + ac11002 + 
                          ac11003 + ac11004 + ac11005 + 
                          ac11006 + ac11007 + ac11008 + 
                          ac11009 + ac11010 + ac11011 + 
                          ac11011 + ac11012 + ac11013 + 
                          ac11014 + ac11016 + ac11017 + 
                          ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                          PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                          PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Alcohol_status_cov +
                          Physical_activity_cov, data = filter(bd_pheno_self, SelfID == 1), na.action = "na.omit")
      }
      else{
        TC_Veg_lm <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                          ac10003 + ac11001 + ac11002 + 
                          ac11003 + ac11004 + ac11005 + 
                          ac11006 + ac11007 + ac11008 + 
                          ac11009 + ac11010 + ac11011 + 
                          ac11011 + ac11012 + ac11013 + 
                          ac11014 + ac11016 + ac11017 + 
                          ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                          PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                          PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Alcohol_status_cov + Smoking_status_cov +
                          Physical_activity_cov, data = filter(bd_pheno_self, SelfID == 1), na.action = "na.omit")
      }
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
      
      #SE_d_TC=
      print(sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      #Lower_CI_d_TC=
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #Higher_CI_d_TC=
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #z_TC=
      (summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])/sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #p_TC=
      print(pnorm(abs(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]) / sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2), lower.tail=FALSE) * 2)
      
      #RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])
      #Lower_CI_RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      #Higher_CI_RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      
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
      
      a1=format(round(summary(TC_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new1=paste(a1," (",b1,", ",c1,")",sep="")
      
      a=format(round(summary(TC_No_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new=paste(a," (",b,", ",c,")",sep="")
      
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      if (e == c("TC")) {
        final_res=data.frame("Total cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("LDL")) {
        final_res=data.frame("LDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("HDL")) {
        final_res=data.frame("HDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("TGs")) {
        final_res=data.frame("Triglycerides, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      }
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_self <- subset (bd_pheno_self, select = -c(PRS_lipds,Lipds,un_sd_PRS_lipds))
      
    }
  }
}








































################################### ###########       No UKB   Graham SE et al. Nature (2021)
################################### Strict
for (x in c("Graham")) {
  final_res=data.frame(NA,NA,"Vegetarian",NA,NA,"Not vegetarian",NA,NA)
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully_main.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  final_res=data.frame("Lipids","n","β (95% CI)","P-value",
                       "n","β (95% CI)","P-value","Pinteractionb")
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully_main.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  final_res=data.frame("Graham SE et al. (excluded UKB data)",NA,NA,NA,NA,NA,NA,NA)
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully_main.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  for (popul in c("EUR", "SAS", "AFR")) {
    
    if(popul == "SAS"){
      final_res=data.frame("CSA",NA,NA,NA,NA,NA,NA,NA)
    }
    else{
      final_res=data.frame(popul,NA,NA,NA,NA,NA,NA,NA)
    }
    
    write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully_main.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    for (e in c("TC","LDL","HDL","TGs")) {
      
      Lipids_name=paste(e, "_",x,"_no_UKB_",popul,"_PRS", sep="")
      
      bd_pheno_strict=bd_pheno_strict %>% select(FID,e,Lipids_name, everything())
      
      ########### PRS
      df_tem <- bd_pheno_strict[colnames(bd_pheno_strict)[3]]
      colnames(df_tem)[1]="PRS_lipds"
      df_tem$un_sd_PRS_lipds=df_tem$PRS_lipds
      df_tem$PRS_lipds=scale(df_tem$PRS_lipds)
      bd_pheno_strict=cbind(bd_pheno_strict,df_tem)
      
      #Standardize genotypic lipid
      # bd_pheno$df_TC <- scale(bd_pheno$TC_Willer_EUR_PRS)
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
      
      #SE_d_TC=
      print(sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      #Lower_CI_d_TC=
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #Higher_CI_d_TC=
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #z_TC=
      (summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])/sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #p_TC=
      print(pnorm(abs(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]) / sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2), lower.tail=FALSE) * 2)
      
      #RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])
      #Lower_CI_RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      #Higher_CI_RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      
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
      
      a1=format(round(summary(TC_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new1=paste(a1," (",b1,", ",c1,")",sep="")
      
      a=format(round(summary(TC_No_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new=paste(a," (",b,", ",c,")",sep="")
      
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      if (e == c("TC")) {
        final_res=data.frame("Total cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("LDL")) {
        final_res=data.frame("LDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("HDL")) {
        final_res=data.frame("HDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("TGs")) {
        final_res=data.frame("Triglycerides, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      }
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully_main.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_strict <- subset (bd_pheno_strict, select = -c(PRS_lipds,Lipds,un_sd_PRS_lipds))
      
    }
  }
}

################################### Self
for (x in c("Graham")) {
  final_res=data.frame(NA,NA,"Vegetarian",NA,NA,"Not vegetarian",NA,NA)
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully_main.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  final_res=data.frame("Lipids","n","β (95% CI)","P-value",
                       "n","β (95% CI)","P-value","Pinteractionb")
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully_main.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  final_res=data.frame("Graham SE et al. (excluded UKB data)",NA,NA,NA,NA,NA,NA,NA)
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully_main.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  for (popul in c("EUR", "SAS", "AFR")) {
    
    if(popul == "SAS"){
      final_res=data.frame("CSA",NA,NA,NA,NA,NA,NA,NA)
    }
    else{
      final_res=data.frame(popul,NA,NA,NA,NA,NA,NA,NA)
    }
    
    write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully_main.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    for (e in c("TC","LDL","HDL","TGs")) {
      
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
      
      #SE_d_TC=
      print(sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      #Lower_CI_d_TC=
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #Higher_CI_d_TC=
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #z_TC=
      (summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])/sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #p_TC=
      print(pnorm(abs(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]) / sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2), lower.tail=FALSE) * 2)
      
      #RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])
      #Lower_CI_RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      #Higher_CI_RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      
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
      
      a1=format(round(summary(TC_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new1=paste(a1," (",b1,", ",c1,")",sep="")
      
      a=format(round(summary(TC_No_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new=paste(a," (",b,", ",c,")",sep="")
      
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      if (e == c("TC")) {
        final_res=data.frame("Total cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("LDL")) {
        final_res=data.frame("LDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("HDL")) {
        final_res=data.frame("HDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("TGs")) {
        final_res=data.frame("Triglycerides, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      }
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully_main.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_self <- subset (bd_pheno_self, select = -c(PRS_lipds,Lipds,un_sd_PRS_lipds))
      
    }
  }
}













































































































################################### Willer CJ et al. Nat. Genet. 2013
################################### Strict
for (x in c("Willer")) {
  final_res=data.frame("Willer CJ et al.",NA,NA,NA,NA,NA,NA,NA)
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')  
  for (popul in c("EUR")) {
    
    final_res=data.frame(popul,NA,NA,NA,NA,NA,NA,NA)
    
    write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    for (e in c("TC","LDL","HDL","TGs")) {
      
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
      
      #SE_d_TC=
      print(sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      #Lower_CI_d_TC=
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #Higher_CI_d_TC=
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #z_TC=
      (summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])/sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #p_TC=
      print(pnorm(abs(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]) / sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2), lower.tail=FALSE) * 2)
      
      #RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])
      #Lower_CI_RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      #Higher_CI_RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      
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
      
      a1=format(round(summary(TC_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new1=paste(a1," (",b1,", ",c1,")",sep="")
      
      a=format(round(summary(TC_No_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new=paste(a," (",b,", ",c,")",sep="")
      
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      if (e == c("TC")) {
        final_res=data.frame("Total cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("LDL")) {
        final_res=data.frame("LDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("HDL")) {
        final_res=data.frame("HDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("TGs")) {
        final_res=data.frame("Triglycerides, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      }
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_fully.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_strict <- subset (bd_pheno_strict, select = -c(PRS_lipds,Lipds,un_sd_PRS_lipds))
      
    }
  }
}

################################### Self
for (x in c("Willer")) {
  final_res=data.frame("Willer CJ et al.",NA,NA,NA,NA,NA,NA,NA)
  
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  for (popul in c("EUR")) {
    
    final_res=data.frame(popul,NA,NA,NA,NA,NA,NA,NA)
    
    write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    for (e in c("TC","LDL","HDL","TGs")) {
      
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
      
      #SE_d_TC=
      print(sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      #Lower_CI_d_TC=
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #Higher_CI_d_TC=
      summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #z_TC=
      (summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])/sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2)
      #p_TC=
      print(pnorm(abs(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]) / sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2), lower.tail=FALSE) * 2)
      
      #RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1])
      #Lower_CI_RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      #Higher_CI_RRR=
      exp(summary(TC_Veg_lm)$coefficients[2,1]-summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*sqrt(summary(TC_Veg_lm)$coefficients[2,2]^2+summary(TC_No_Veg_lm)$coefficients[2,2]^2))
      
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
      
      a1=format(round(summary(TC_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c1=format(round(summary(TC_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval1=ifelse(summary(TC_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new1=paste(a1," (",b1,", ",c1,")",sep="")
      
      a=format(round(summary(TC_No_Veg_lm)$coefficients[2,1],digit=3), nsmall = 3)
      b=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      c=format(round(summary(TC_No_Veg_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_No_Veg_lm)$coefficients[2,2],digit=3), nsmall = 3)
      pval=ifelse(summary(TC_No_Veg_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_No_Veg_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_No_Veg_lm)$coefficients[2,4],digit=3), nsmall = 3))
      new=paste(a," (",b,", ",c,")",sep="")
      
      pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
      
      if (e == c("TC")) {
        final_res=data.frame("Total cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("LDL")) {
        final_res=data.frame("LDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("HDL")) {
        final_res=data.frame("HDL cholesterol, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      } else if (e == c("TGs")) {
        final_res=data.frame("Triglycerides, SD",nobs(TC_Veg_lm),new1,pval1,
                             nobs(TC_No_Veg_lm),new,pval,pval_int)
      }
      
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/final_res_self_fully.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      bd_pheno_self <- subset (bd_pheno_self, select = -c(PRS_lipds,Lipds,un_sd_PRS_lipds))
      
    }
  }
}











































