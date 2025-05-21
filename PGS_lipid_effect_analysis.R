# PGS-lipid effect analysis

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
load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_excluded_data_PRS_factorized.RData")


bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS_factorized)
bd_pheno_self <- data.frame(veg_excluded_data_self_PRS_factorized)
bd_pheno_full <- data.frame(veg_excluded_data_PRS_factorized)

bd_pheno_full_EUR <- bd_pheno_full[(bd_pheno_full$pop=="EUR"), ]
bd_pheno_full_SAS <- bd_pheno_full[(bd_pheno_full$pop=="CSA"), ]
bd_pheno_full_AFR <- bd_pheno_full[(bd_pheno_full$pop=="AFR"), ]
bd_pheno_full_EAS <- bd_pheno_full[(bd_pheno_full$pop=="EAS"), ]

bd_pheno_EUR <- bd_pheno_strict[(bd_pheno_strict$pop=="EUR"), ]
bd_pheno_SAS <- bd_pheno_strict[(bd_pheno_strict$pop=="CSA"), ]
bd_pheno_AFR <- bd_pheno_strict[(bd_pheno_strict$pop=="AFR"), ]
bd_pheno_EAS <- bd_pheno_strict[(bd_pheno_strict$pop=="EAS"), ]

bd_pheno_self_EUR <- bd_pheno_self[(bd_pheno_self$pop=="EUR"), ]
bd_pheno_self_SAS <- bd_pheno_self[(bd_pheno_self$pop=="CSA"), ]
bd_pheno_self_AFR <- bd_pheno_self[(bd_pheno_self$pop=="AFR"), ]
bd_pheno_self_EAS <- bd_pheno_self[(bd_pheno_self$pop=="EAS"), ]

############################################################################################################################################################### FULL

bd_pheno_full_tem1 = bd_pheno_full
for (percent_num in c(33)) {
  final_res=data.frame(NA,"Total population","Vegetarian population","β","SE","P")
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_full_120324.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  for (popul in c("EUR", "SAS", "AFR", "EAS")) {
    if (popul == c("EUR")) {
      final_res=data.frame("European",NA,NA,NA,NA,NA,NA,NA)
    } else if (popul == c("SAS")) {
      final_res=data.frame("Central/South Asian",NA,NA,NA,NA,NA,NA,NA)
    } else if (popul == c("AFR")) {
      final_res=data.frame("African",NA,NA,NA,NA,NA,NA,NA)
    } else if (popul == c("EAS")) {
      final_res=data.frame("East Asian",NA,NA,NA,NA,NA,NA,NA)
    }
    write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_full_120324.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    for (e in c("TC","LDL","HDL","TGs")) {
      if (e == c("TC")) {
        final_res=data.frame("Total cholesterol",NA,NA,NA,NA,NA,NA,NA)
      } else if (e == c("LDL")) {
        final_res=data.frame("LDL cholesterol",NA,NA,NA,NA,NA,NA,NA)
      } else if (e == c("HDL")) {
        final_res=data.frame("HDL cholesterol",NA,NA,NA,NA,NA,NA,NA)
      } else if (e == c("TGs")) {
        final_res=data.frame("Triglycerides",NA,NA,NA,NA,NA,NA,NA)
      }
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_full_120324.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      for (x in c("Graha","Graham","Willer")) {
        
        if (x == c("Graha")) {
          Lipids_name=paste(e, "_",x,"m_no_UKB_",popul,"_PRS", sep="")
        } else if (x == c("Graham")) {
          Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
        } else if (x == c("Willer")) {
          Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
        } 
        
        bd_pheno_full=bd_pheno_full_tem1
        
        if (Lipids_name %in% c(colnames(bd_pheno_full))) {

          bd_pheno_full=bd_pheno_full %>% select(FID,e,Lipids_name, everything())

          df_tem <- bd_pheno_full[colnames(bd_pheno_full)[3]]
          colnames(df_tem)[1]="PRS_lipds"
          df_tem$un_sd_PRS_lipds=df_tem$PRS_lipds
          df_tem$PRS_lipds=scale(df_tem$PRS_lipds)
          bd_pheno_full=cbind(bd_pheno_full,df_tem)
          
          df_tem <- bd_pheno_full[colnames(bd_pheno_full)[2]]
          colnames(df_tem)[1]="Lipds"
          bd_pheno_full=cbind(bd_pheno_full,df_tem)
          
          bd_pheno_full <- bd_pheno_full[is.na(bd_pheno_full$PRS_lipds)==F,]
          
          if (percent_num == c(33)) {
            bd_pheno_full$tertile <- ntile(bd_pheno_full$PRS_lipds, 3)
            
            bd_pheno_full$PRS_lipds_tem=ifelse(is.na(bd_pheno_full$tertile)==F & bd_pheno_full$tertile==1, "Low", NA)
            bd_pheno_full$PRS_lipds_tem=ifelse(is.na(bd_pheno_full$tertile)==F & bd_pheno_full$tertile==3, "High", bd_pheno_full$PRS_lipds_tem)
            bd_pheno_full$PRS_lipds_tem=ifelse(is.na(bd_pheno_full$PRS_lipds_tem)==T, "Intermediate", bd_pheno_full$PRS_lipds_tem)
            
            bd_pheno_full$PRS_lipds_tem <- as.factor(bd_pheno_full$PRS_lipds_tem)
            bd_pheno_full$PRS_lipds_tem<- factor(bd_pheno_full$PRS_lipds_tem, levels = c("Low","Intermediate","High"))
            
          } 
          
          bd_pheno_full$Strict <- as.factor(bd_pheno_full$Strict)
          bd_pheno_full$Strict<- factor(bd_pheno_full$Strict, levels = c(0,1))
          
          bd_pheno_full_tem <- bd_pheno_full %>% select(Lipds,PRS_lipds,Sex,Age,Array,
                                                            ac10003,ac11001,ac11002,
                                                            ac11003,ac11004,ac11005,
                                                            ac11006,ac11007,ac11008,
                                                            ac11009,ac11010,ac11011,
                                                            ac11011,ac11012,ac11013,
                                                            ac11014,ac11016,ac11017,
                                                            ac11018,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,
                                                            PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,
                                                            PC19,PC20,BMI,Townsend,statin_sum_standardized,Smoking_status_cov,Alcohol_status_cov,
                                                            Physical_activity_cov,Strict,PRS_lipds_tem)
          
          TC_Veg_lm_int <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                                ac10003 + ac11001 + ac11002 +
                                ac11003 + ac11004 + ac11005 +
                                ac11006 + ac11007 + ac11008 +
                                ac11009 + ac11010 + ac11011 +
                                ac11011 + ac11012 + ac11013 +
                                ac11014 + ac11016 + ac11017 +
                                ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +
                                PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                                PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov +
                                Physical_activity_cov, data = bd_pheno_full, na.action = "na.omit")
          
          veg_count <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                            ac10003 + ac11001 + ac11002 +
                            ac11003 + ac11004 + ac11005 +
                            ac11006 + ac11007 + ac11008 +
                            ac11009 + ac11010 + ac11011 +
                            ac11011 + ac11012 + ac11013 +
                            ac11014 + ac11016 + ac11017 +
                            ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +
                            PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                            PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Alcohol_status_cov +
                            Physical_activity_cov, data = filter(bd_pheno_full, Strict == 1), na.action = "na.omit")
          
          beta_Veg=ifelse(abs(summary(TC_Veg_lm_int)$coefficients[2,1])<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,1], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,1],digit=3), nsmall = 3))
          se_Veg=ifelse(summary(TC_Veg_lm_int)$coefficients[2,2]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,2], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,2],digit=3), nsmall = 3))
          pval_Veg=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
          hr_Veg="—"
          a1=format(round(summary(TC_Veg_lm_int)$coefficients[2,1],digit=3), nsmall = 3)
          b1=format(round(summary(TC_Veg_lm_int)$coefficients[2,1]+qnorm(0.025)*summary(TC_Veg_lm_int)$coefficients[2,2],digit=3), nsmall = 3)
          c1=format(round(summary(TC_Veg_lm_int)$coefficients[2,1]+qnorm(0.975)*summary(TC_Veg_lm_int)$coefficients[2,2],digit=3), nsmall = 3)
          new1=paste(a1," (",b1,", ",c1,")",sep="")
          if (x == c("Graha")) {
            final_res_2=data.frame("     Graham SE et al. (excluded UKB data)",nobs(TC_Veg_lm_int),nobs(veg_count),new1,se_Veg,pval_Veg)
          } else if (x == c("Graham")) {
            final_res_2=data.frame("     Graham SE et al. (included UKB data)",nobs(TC_Veg_lm_int),nobs(veg_count),new1,se_Veg,pval_Veg)
          } else if (x == c("Willer")) {
            final_res_2=data.frame("     Willer CJ et a",nobs(TC_Veg_lm_int),nobs(veg_count),new1,se_Veg,pval_Veg)
          } 
          
          write.table(final_res_2, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_full_120324.txt", col.names = FALSE, append = TRUE,
                      row.names = F, quote = FALSE, na = "",sep='\t')
          
          bd_pheno_full <- subset (bd_pheno_full, select = -c(PRS_lipds,Lipds,un_sd_PRS_lipds))
        }
      }
    }
  }
}






############################################################################################################################################################### STRICT 
bd_pheno_strict_tem1 = bd_pheno_strict
for (percent_num in c(33)) {
  final_res=data.frame(NA,"Total population","Vegetarian population","β","SE","P")
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_120324.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')

      for (popul in c("EUR", "SAS", "AFR", "EAS")) {
        if (popul == c("EUR")) {
          final_res=data.frame("European",NA,NA,NA,NA,NA,NA,NA)
        } else if (popul == c("SAS")) {
          final_res=data.frame("Central/South Asian",NA,NA,NA,NA,NA,NA,NA)
        } else if (popul == c("AFR")) {
          final_res=data.frame("African",NA,NA,NA,NA,NA,NA,NA)
        } else if (popul == c("EAS")) {
          final_res=data.frame("East Asian",NA,NA,NA,NA,NA,NA,NA)
        }
        write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_120324.txt", col.names = FALSE, append = TRUE,
                    row.names = F, quote = FALSE, na = "",sep='\t')
        
        for (e in c("TC","LDL","HDL","TGs")) {
          if (e == c("TC")) {
            final_res=data.frame("Total cholesterol",NA,NA,NA,NA,NA,NA,NA)
          } else if (e == c("LDL")) {
            final_res=data.frame("LDL cholesterol",NA,NA,NA,NA,NA,NA,NA)
          } else if (e == c("HDL")) {
            final_res=data.frame("HDL cholesterol",NA,NA,NA,NA,NA,NA,NA)
          } else if (e == c("TGs")) {
            final_res=data.frame("Triglycerides",NA,NA,NA,NA,NA,NA,NA)
          }
          write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_120324.txt", col.names = FALSE, append = TRUE,
                      row.names = F, quote = FALSE, na = "",sep='\t')
          
          for (x in c("Graha","Graham","Willer")) {
            
        if (x == c("Graha")) {
          Lipids_name=paste(e, "_",x,"m_no_UKB_",popul,"_PRS", sep="")
        } else if (x == c("Graham")) {
          Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
        } else if (x == c("Willer")) {
          Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
        } 
        
        bd_pheno_strict=bd_pheno_strict_tem1
        
        
        if (Lipids_name %in% c(colnames(bd_pheno_strict))) {
          
          bd_pheno_strict=bd_pheno_strict %>% select(FID,e,Lipids_name, everything())

          df_tem <- bd_pheno_strict[colnames(bd_pheno_strict)[3]]
          colnames(df_tem)[1]="PRS_lipds"
          df_tem$un_sd_PRS_lipds=df_tem$PRS_lipds
          df_tem$PRS_lipds=scale(df_tem$PRS_lipds)
          bd_pheno_strict=cbind(bd_pheno_strict,df_tem)
          
          #Standardize genotypic lipid
          df_tem <- bd_pheno_strict[colnames(bd_pheno_strict)[2]]
          colnames(df_tem)[1]="Lipds"
          bd_pheno_strict=cbind(bd_pheno_strict,df_tem)
          
          bd_pheno_strict <- bd_pheno_strict[is.na(bd_pheno_strict$PRS_lipds)==F,]
          
          if (percent_num == c(33)) {
            bd_pheno_strict$tertile <- ntile(bd_pheno_strict$PRS_lipds, 3)
            
            bd_pheno_strict$PRS_lipds_tem=ifelse(is.na(bd_pheno_strict$tertile)==F & bd_pheno_strict$tertile==1, "Low", NA)
            bd_pheno_strict$PRS_lipds_tem=ifelse(is.na(bd_pheno_strict$tertile)==F & bd_pheno_strict$tertile==3, "High", bd_pheno_strict$PRS_lipds_tem)
            bd_pheno_strict$PRS_lipds_tem=ifelse(is.na(bd_pheno_strict$PRS_lipds_tem)==T, "Intermediate", bd_pheno_strict$PRS_lipds_tem)
            bd_pheno_strict$PRS_lipds_tem <- as.factor(bd_pheno_strict$PRS_lipds_tem)
            bd_pheno_strict$PRS_lipds_tem<- factor(bd_pheno_strict$PRS_lipds_tem, levels = c("Low","Intermediate","High"))
            
          } 
          
          bd_pheno_strict$Strict <- as.factor(bd_pheno_strict$Strict)
          bd_pheno_strict$Strict<- factor(bd_pheno_strict$Strict, levels = c(0,1))
          
          bd_pheno_strict_tem <- bd_pheno_strict %>% select(Lipds,PRS_lipds,Sex,Age,Array,
                                                            ac10003,ac11001,ac11002,
                                                            ac11003,ac11004,ac11005,
                                                            ac11006,ac11007,ac11008,
                                                            ac11009,ac11010,ac11011,
                                                            ac11011,ac11012,ac11013,
                                                            ac11014,ac11016,ac11017,
                                                            ac11018,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,
                                                            PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,
                                                            PC19,PC20,BMI,Townsend,statin_sum_standardized,Smoking_status_cov,Alcohol_status_cov,
                                                            Physical_activity_cov,Strict,PRS_lipds_tem)
          
          TC_Veg_lm_int <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
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
          
          veg_count <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
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
     
          beta_Veg=ifelse(abs(summary(TC_Veg_lm_int)$coefficients[2,1])<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,1], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,1],digit=3), nsmall = 3))
          se_Veg=ifelse(summary(TC_Veg_lm_int)$coefficients[2,2]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,2], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,2],digit=3), nsmall = 3))
          pval_Veg=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
          hr_Veg="—"
          a1=format(round(summary(TC_Veg_lm_int)$coefficients[2,1],digit=3), nsmall = 3)
          b1=format(round(summary(TC_Veg_lm_int)$coefficients[2,1]+qnorm(0.025)*summary(TC_Veg_lm_int)$coefficients[2,2],digit=3), nsmall = 3)
          c1=format(round(summary(TC_Veg_lm_int)$coefficients[2,1]+qnorm(0.975)*summary(TC_Veg_lm_int)$coefficients[2,2],digit=3), nsmall = 3)
          new1=paste(a1," (",b1,", ",c1,")",sep="")
          if (x == c("Graha")) {
            final_res_2=data.frame("     Graham SE et al. (excluded UKB data)",nobs(TC_Veg_lm_int),nobs(veg_count),new1,se_Veg,pval_Veg)
          } else if (x == c("Graham")) {
            final_res_2=data.frame("     Graham SE et al. (included UKB data)",nobs(TC_Veg_lm_int),nobs(veg_count),new1,se_Veg,pval_Veg)
          } else if (x == c("Willer")) {
            final_res_2=data.frame("     Willer CJ et a",nobs(TC_Veg_lm_int),nobs(veg_count),new1,se_Veg,pval_Veg)
          } 
          
          write.table(final_res_2, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_120324.txt", col.names = FALSE, append = TRUE,
                      row.names = F, quote = FALSE, na = "",sep='\t')
          
          bd_pheno_strict <- subset (bd_pheno_strict, select = -c(PRS_lipds,Lipds,un_sd_PRS_lipds))
        }
      }
    }
  }
}




############################################################################################################################################################### SELF 
bd_pheno_self_tem1 = bd_pheno_self
for (percent_num in c(33)) {
  final_res=data.frame(NA,"Total population","Vegetarian population","β","SE","P")
  write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_self_120324.txt", col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  
  for (popul in c("EUR", "SAS", "AFR", "EAS")) {
    if (popul == c("EUR")) {
      final_res=data.frame("European",NA,NA,NA,NA,NA,NA,NA)
    } else if (popul == c("SAS")) {
      final_res=data.frame("Central/South Asian",NA,NA,NA,NA,NA,NA,NA)
    } else if (popul == c("AFR")) {
      final_res=data.frame("African",NA,NA,NA,NA,NA,NA,NA)
    } else if (popul == c("EAS")) {
      final_res=data.frame("East Asian",NA,NA,NA,NA,NA,NA,NA)
    }
    write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_self_120324.txt", col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    
    for (e in c("TC","LDL","HDL","TGs")) {
      if (e == c("TC")) {
        final_res=data.frame("Total cholesterol",NA,NA,NA,NA,NA,NA,NA)
      } else if (e == c("LDL")) {
        final_res=data.frame("LDL cholesterol",NA,NA,NA,NA,NA,NA,NA)
      } else if (e == c("HDL")) {
        final_res=data.frame("HDL cholesterol",NA,NA,NA,NA,NA,NA,NA)
      } else if (e == c("TGs")) {
        final_res=data.frame("Triglycerides",NA,NA,NA,NA,NA,NA,NA)
      }
      write.table(final_res, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_self_120324.txt", col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      
      for (x in c("Graha","Graham","Willer")) {
        
        if (x == c("Graha")) {
          Lipids_name=paste(e, "_",x,"m_no_UKB_",popul,"_PRS", sep="")
        } else if (x == c("Graham")) {
          Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
        } else if (x == c("Willer")) {
          Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
        } 
        
        bd_pheno_self=bd_pheno_self_tem1
        
        
        if (Lipids_name %in% c(colnames(bd_pheno_self))) {
          
          bd_pheno_self=bd_pheno_self %>% select(FID,e,Lipids_name, everything())

          df_tem <- bd_pheno_self[colnames(bd_pheno_self)[3]]
          colnames(df_tem)[1]="PRS_lipds"
          df_tem$un_sd_PRS_lipds=df_tem$PRS_lipds
          df_tem$PRS_lipds=scale(df_tem$PRS_lipds)
          bd_pheno_self=cbind(bd_pheno_self,df_tem)
          
          #Standardize genotypic lipid
          df_tem <- bd_pheno_self[colnames(bd_pheno_self)[2]]
          colnames(df_tem)[1]="Lipds"
          bd_pheno_self=cbind(bd_pheno_self,df_tem)
          
          bd_pheno_self <- bd_pheno_self[is.na(bd_pheno_self$PRS_lipds)==F,]
          
          if (percent_num == c(33)) {
            bd_pheno_self$tertile <- ntile(bd_pheno_self$PRS_lipds, 3)
            
            bd_pheno_self$PRS_lipds_tem=ifelse(is.na(bd_pheno_self$tertile)==F & bd_pheno_self$tertile==1, "Low", NA)
            bd_pheno_self$PRS_lipds_tem=ifelse(is.na(bd_pheno_self$tertile)==F & bd_pheno_self$tertile==3, "High", bd_pheno_self$PRS_lipds_tem)
            bd_pheno_self$PRS_lipds_tem=ifelse(is.na(bd_pheno_self$PRS_lipds_tem)==T, "Intermediate", bd_pheno_self$PRS_lipds_tem)
            
            bd_pheno_self$PRS_lipds_tem <- as.factor(bd_pheno_self$PRS_lipds_tem)
            bd_pheno_self$PRS_lipds_tem<- factor(bd_pheno_self$PRS_lipds_tem, levels = c("Low","Intermediate","High"))
            
          } 
          
          bd_pheno_self$Self <- as.factor(bd_pheno_self$Self)
          bd_pheno_self$Self<- factor(bd_pheno_self$Self, levels = c(0,1))
          
          bd_pheno_self_tem <- bd_pheno_self %>% select(Lipds,PRS_lipds,Sex,Age,Array,
                                                            ac10003,ac11001,ac11002,
                                                            ac11003,ac11004,ac11005,
                                                            ac11006,ac11007,ac11008,
                                                            ac11009,ac11010,ac11011,
                                                            ac11011,ac11012,ac11013,
                                                            ac11014,ac11016,ac11017,
                                                            ac11018,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,
                                                            PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,
                                                            PC19,PC20,BMI,Townsend,statin_sum_standardized,Smoking_status_cov,Alcohol_status_cov,
                                                            Physical_activity_cov,Self,PRS_lipds_tem)
          
          TC_Veg_lm_int <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
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
          
          veg_count <- lm(Lipds ~ PRS_lipds + Sex + Age + I(Age^2) + Array + 
                            ac10003 + ac11001 + ac11002 +
                            ac11003 + ac11004 + ac11005 +
                            ac11006 + ac11007 + ac11008 +
                            ac11009 + ac11010 + ac11011 +
                            ac11011 + ac11012 + ac11013 +
                            ac11014 + ac11016 + ac11017 +
                            ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +
                            PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                            PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Alcohol_status_cov +
                            Physical_activity_cov, data = filter(bd_pheno_self, Self == 1), na.action = "na.omit")
          
          beta_Veg=ifelse(abs(summary(TC_Veg_lm_int)$coefficients[2,1])<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,1], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,1],digit=3), nsmall = 3))
          se_Veg=ifelse(summary(TC_Veg_lm_int)$coefficients[2,2]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,2], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,2],digit=3), nsmall = 3))
          pval_Veg=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
          hr_Veg="—"
          a1=format(round(summary(TC_Veg_lm_int)$coefficients[2,1],digit=3), nsmall = 3)
          b1=format(round(summary(TC_Veg_lm_int)$coefficients[2,1]+qnorm(0.025)*summary(TC_Veg_lm_int)$coefficients[2,2],digit=3), nsmall = 3)
          c1=format(round(summary(TC_Veg_lm_int)$coefficients[2,1]+qnorm(0.975)*summary(TC_Veg_lm_int)$coefficients[2,2],digit=3), nsmall = 3)
          new1=paste(a1," (",b1,", ",c1,")",sep="")
          if (x == c("Graha")) {
            final_res_2=data.frame("     Graham SE et al. (excluded UKB data)",nobs(TC_Veg_lm_int),nobs(veg_count),new1,se_Veg,pval_Veg)
          } else if (x == c("Graham")) {
            final_res_2=data.frame("     Graham SE et al. (included UKB data)",nobs(TC_Veg_lm_int),nobs(veg_count),new1,se_Veg,pval_Veg)
          } else if (x == c("Willer")) {
            final_res_2=data.frame("     Willer CJ et a",nobs(TC_Veg_lm_int),nobs(veg_count),new1,se_Veg,pval_Veg)
          } 
          
          write.table(final_res_2, file= "/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/lipid_PGS_beta_self_120324.txt", col.names = FALSE, append = TRUE,
                      row.names = F, quote = FALSE, na = "",sep='\t')
          
          bd_pheno_self <- subset (bd_pheno_self, select = -c(PRS_lipds,Lipds,un_sd_PRS_lipds))
        }
      }
    }
  }
}

