# supp table 3
'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr)
library(tidyverse)
library(readr)
library("lubridate")
library("survival")


# Pathway=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/")
# Pathway_out=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/")

# load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_excluded_data_strict_PRS_factorized.RData")
# load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_excluded_data_self_PRS_factorized.RData")

bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS_factorized)
bd_pheno_self <- data.frame(veg_excluded_data_self_PRS_factorized)

bd_pheno_EUR <- bd_pheno_strict[(bd_pheno_strict$pop=="EUR"), ]
bd_pheno_SAS <- bd_pheno_strict[(bd_pheno_strict$pop=="CSA"), ]
bd_pheno_AFR <- bd_pheno_strict[(bd_pheno_strict$pop=="AFR"), ]
bd_pheno_EAS <- bd_pheno_strict[(bd_pheno_strict$pop=="EAS"), ]

bd_pheno_self_EUR <- bd_pheno_self[(bd_pheno_self$pop=="EUR"), ]
bd_pheno_self_SAS <- bd_pheno_self[(bd_pheno_self$pop=="CSA"), ]
bd_pheno_self_AFR <- bd_pheno_self[(bd_pheno_self$pop=="AFR"), ]
bd_pheno_self_EAS <- bd_pheno_self[(bd_pheno_self$pop=="EAS"), ]

#Strict


final_res_2=data.frame(NA,"Total population","Vegetarian population","β","SE","P")

write.table(final_res_2, file= "/Users/angelage/Desktop/veg_lipid_beta_table_2.txt", col.names = FALSE, append = TRUE,
            row.names = F, quote = FALSE, na = "",sep='\t')
for (e in c("TC","LDL","HDL","TGs")) {
  for (popul in c("EUR", "AFR", "SAS", "EAS")) {
    if(popul == "EUR"){
      if (e == c("TC")) {
        final_res=data.frame("Total cholesterol, SD",NA,NA,NA,NA,NA)
      } else if (e == c("LDL")) {
        final_res=data.frame("LDL cholesterol, SD",NA,NA,NA,NA,NA)
      } else if (e == c("HDL")) {
        final_res=data.frame("HDL cholesterol, SD",NA,NA,NA,NA,NA)
      } else if (e == c("TGs")) {
        final_res=data.frame("Triglycerides, SD",NA,NA,NA,NA,NA)
      }
      write.table(final_res, file= "/Users/angelage/Desktop/veg_lipid_beta_table_2.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
    }
    if (popul == "SAS") {
      final_res_2=data.frame("CSA",NA,NA,NA,NA,NA)
    }
    else{
      final_res_2=data.frame(popul,NA,NA,NA,NA,NA)
    }
    write.table(final_res_2, file= "/Users/angelage/Desktop/veg_lipid_beta_table_2.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    if (popul == "EUR") {
      bd_pheno_copy <- bd_pheno_EUR
      bd_pheno_tem1=bd_pheno_EUR
      bd_pheno_tem=bd_pheno_EUR
    } else if (popul == "AFR") {
      bd_pheno_copy <- bd_pheno_AFR
      bd_pheno_tem1=bd_pheno_AFR
      bd_pheno_tem=bd_pheno_AFR
    } else if (popul == "EAS") {
      bd_pheno_copy <- bd_pheno_EAS
      bd_pheno_tem1=bd_pheno_EAS
      bd_pheno_tem=bd_pheno_EAS
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
    
    TC_Veg_lm_int <- lm(Lipds ~ Strict + Sex + Age + I(Age^2) + Array + 
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
    
    veg_count <- lm(Lipds ~ Strict + Sex + Age + I(Age^2) + Array + 
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
    
    beta_Veg=summary(TC_Veg_lm_int)$coefficients[2,1]
    se_Veg=summary(TC_Veg_lm_int)$coefficients[2,2]
    pval_Veg=summary(TC_Veg_lm_int)$coefficients[2,4]
    hr_Veg="—"
    a1=format(summary(TC_Veg_lm_int)$coefficients[2,1])
    b1=format(summary(TC_Veg_lm_int)$coefficients[2,1]+qnorm(0.025)*summary(TC_Veg_lm_int)$coefficients[2,2])
    c1=format(summary(TC_Veg_lm_int)$coefficients[2,1]+qnorm(0.975)*summary(TC_Veg_lm_int)$coefficients[2,2])
    new1=paste(a1," (",b1,", ",c1,")",sep="")
    
    final_res_2=data.frame("Vegetarian",nobs(TC_Veg_lm_int),nobs(veg_count),new1,se_Veg,pval_Veg)
    
    write.table(final_res_2, file= "/Users/angelage/Desktop/veg_lipid_beta_table_2.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
  }
}



#Self


final_res_2=data.frame(NA,"Total population","Vegetarian population","β","SE","P")

write.table(final_res_2, file= "/Users/angelage/Desktop/veg_lipid_beta_table_2_self.txt", col.names = FALSE, append = TRUE,
            row.names = F, quote = FALSE, na = "",sep='\t')


for (e in c("TC","LDL","HDL","TGs")) {
  #NO EAS FOR NO UKB
  for (popul in c("EUR", "AFR", "SAS", "EAS")) {
    
    if(popul == "EUR"){
      if (e == c("TC")) {
        final_res=data.frame("Total cholesterol, SD",NA,NA,NA,NA,NA)
      } else if (e == c("LDL")) {
        final_res=data.frame("LDL cholesterol, SD",NA,NA,NA,NA,NA)
      } else if (e == c("HDL")) {
        final_res=data.frame("HDL cholesterol, SD",NA,NA,NA,NA,NA)
      } else if (e == c("TGs")) {
        final_res=data.frame("Triglycerides, SD",NA,NA,NA,NA,NA)
      }
      write.table(final_res, file= "/Users/angelage/Desktop/veg_lipid_beta_table_2_self.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
    }
    if (popul == "SAS") {
      final_res_2=data.frame("CSA",NA,NA,NA,NA,NA)
    }
    else{
      final_res_2=data.frame(popul,NA,NA,NA,NA,NA)
    }
    write.table(final_res_2, file= "/Users/angelage/Desktop/veg_lipid_beta_table_2_self.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    if (popul == "EUR") {
      bd_pheno_copy <- bd_pheno_self_EUR
      bd_pheno_tem1=bd_pheno_self_EUR
      bd_pheno_tem=bd_pheno_self_EUR
    } else if (popul == "AFR") {
      bd_pheno_copy <- bd_pheno_self_AFR
      bd_pheno_tem1=bd_pheno_self_AFR
      bd_pheno_tem=bd_pheno_self_AFR
    } else if (popul == "EAS") {
      bd_pheno_copy <- bd_pheno_self_EAS
      bd_pheno_tem1=bd_pheno_self_EAS
      bd_pheno_tem=bd_pheno_self_EAS
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
    
    TC_Veg_lm_int <- lm(Lipds ~ SelfID + Sex + Age + I(Age^2) + Array + 
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
    veg_count <- lm(Lipds ~ SelfID + Sex + Age + I(Age^2) + Array + 
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
    
    beta_Veg=summary(TC_Veg_lm_int)$coefficients[2,1]
    se_Veg=summary(TC_Veg_lm_int)$coefficients[2,2]
    pval_Veg=summary(TC_Veg_lm_int)$coefficients[2,4]
    hr_Veg="—"
    a1=summary(TC_Veg_lm_int)$coefficients[2,1]
    b1=summary(TC_Veg_lm_int)$coefficients[2,1]+qnorm(0.025)*summary(TC_Veg_lm_int)$coefficients[2,2]
    c1=summary(TC_Veg_lm_int)$coefficients[2,1]+qnorm(0.975)*summary(TC_Veg_lm_int)$coefficients[2,2]
    new1=paste(a1," (",b1,", ",c1,")",sep="")
    
    final_res_2=data.frame("Vegetarian",nobs(TC_Veg_lm_int),nobs(veg_count),new1,se_Veg,pval_Veg)
    
    write.table(final_res_2, file= "/Users/angelage/Desktop/veg_lipid_beta_table_2_self.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
  }
}













