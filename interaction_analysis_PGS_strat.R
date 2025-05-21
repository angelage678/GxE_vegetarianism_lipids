# stratified by PGS interaction analysis

'%ni%' <- Negate('%in%')
library(plyr)
library(dplyr) 
library(tidyverse)
library(readr)
library("lubridate")
library("survival")

Pathway=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/")
Pathway_out=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/strat_results")
load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_excluded_data_strict_PRS_factorized.RData")
load("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/veg_excluded_data_self_PRS_factorized.RData")

bd_pheno_strict <- data.frame(veg_excluded_data_strict_PRS_factorized)
bd_pheno_self <- data.frame(veg_excluded_data_self_PRS_factorized)

bd_pheno_strict <- bd_pheno_strict[is.na(bd_pheno_strict$TC_Graham_EUR_PRS)==F | is.na(bd_pheno_strict$TC_Graham_AFR_PRS)==F | 
                                     is.na(bd_pheno_strict$TC_Graham_EAS_PRS)==F | is.na(bd_pheno_strict$TC_Graham_SAS_PRS)==F ,]
bd_pheno_self <- bd_pheno_self[is.na(bd_pheno_self$TC_Graham_EUR_PRS)==F | is.na(bd_pheno_self$TC_Graham_AFR_PRS)==F | 
                                 is.na(bd_pheno_self$TC_Graham_EAS_PRS)==F | is.na(bd_pheno_self$TC_Graham_SAS_PRS)==F ,]

bd_pheno_strict_copy <- bd_pheno_strict
bd_pheno_self_copy <- bd_pheno_self

########################################################################################################################## STRICT 
bd_pheno_strict_tem1 = bd_pheno_strict
for (percent_num in c(33)) {
  final_res=data.frame(NA,"Low PGS",NA,"Intermediate PGS",NA,"High PGS",NA,NA)
  write.table(final_res, file= paste(Pathway_out,"strat1204.txt", sep=""), col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  final_res=data.frame(NA,"Beta (95% CI)","P-value","Beta (95% CI)","P-value",
                       "Beta (95% CI)","P-value","Pinteractionb")
  write.table(final_res, file= paste(Pathway_out,"strat1204.txt", sep=""), col.names = FALSE, append = TRUE,
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
    write.table(final_res, file= paste(Pathway_out,"strat1204.txt", sep=""), col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    for (x in c("Graha","Graham","Willer")) {
      if (x == c("Graha")) {
        final_res=data.frame("     Graham SE et al. (excluded UKB data)",NA,NA,NA,NA,NA,NA,NA)
      } else if (x == c("Graham")) {
        final_res=data.frame("     Graham SE et al. (included UKB data)",NA,NA,NA,NA,NA,NA,NA)
      } else if (x == c("Willer")) {
        final_res=data.frame("     Willer CJ et al.",NA,NA,NA,NA,NA,NA,NA)
      } 
      write.table(final_res, file= paste(Pathway_out,"strat1204.txt", sep=""), col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      for (popul in c("EUR")) {
        if (x == c("Graha")) {
          final_res=data.frame("Graham SE et al. (excluded UKB data)",NA,NA,NA,NA,NA,NA,NA)
          Lipids_name=paste(e, "_",x,"m_no_UKB_",popul,"_PRS", sep="")
        } else if (x == c("Graham")) {
          final_res=data.frame("Graham SE et al. (included UKB data)",NA,NA,NA,NA,NA,NA,NA)
          Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
        } else if (x == c("Willer")) {
          final_res=data.frame("Willer CJ et al.",NA,NA,NA,NA,NA,NA,NA)
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
          
          bd_pheno_strict$Strict_factor <- as.factor(bd_pheno_strict$Strict)
          bd_pheno_strict$Strict_factor<- factor(bd_pheno_strict$Strict, levels = c(0,1))
          
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
                                                            Physical_activity_cov,Strict_factor,Strict,PRS_lipds_tem)
          
          TC_Low_lm <- lm(Lipds ~ Strict + Sex + Age + I(Age^2) + Array + 
                            ac10003 + ac11001 + ac11002 + 
                            ac11003 + ac11004 + ac11005 + 
                            ac11006 + ac11007 + ac11008 + 
                            ac11009 + ac11010 + ac11011 + 
                            ac11011 + ac11012 + ac11013 + 
                            ac11014 + ac11016 + ac11017 + 
                            ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                            PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                            PC19 + PC20 + BMI + Townsend +
                            statin_sum_standardized + 
                            Smoking_status_cov + 
                            Alcohol_status_cov + 
                            Physical_activity_cov, data = filter(bd_pheno_strict_tem, PRS_lipds_tem == "Low"))
          TC_Intermediate_lm <- lm(Lipds ~ Strict + Sex + Age + I(Age^2) + Array +
                                     ac10003 + ac11001 + ac11002 +
                                     ac11003 + ac11004 + ac11005 +
                                     ac11006 + ac11007 + ac11008 +
                                     ac11009 + ac11010 + ac11011 +
                                     ac11011 + ac11012 + ac11013 +
                                     ac11014 + ac11016 + ac11017 +
                                     ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +
                                     PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                                     PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov +
                                     Physical_activity_cov, data = filter(bd_pheno_strict_tem, PRS_lipds_tem == "Intermediate"))
          TC_High_lm <- lm(Lipds ~ Strict + Sex + Age + I(Age^2) + Array +
                             ac10003 + ac11001 + ac11002 +
                             ac11003 + ac11004 + ac11005 +
                             ac11006 + ac11007 + ac11008 +
                             ac11009 + ac11010 + ac11011 +
                             ac11011 + ac11012 + ac11013 +
                             ac11014 + ac11016 + ac11017 +
                             ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +
                             PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                             PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov +
                             Physical_activity_cov, data = filter(bd_pheno_strict_tem, PRS_lipds_tem == "High"))
        
          if (percent_num == 33) {
            bd_pheno_strict$PRS_lipds_num <- ifelse(bd_pheno_strict$tertile == 1, 0,
                                                        ifelse(bd_pheno_strict$tertile == 3, 1, NA))
          }
          
          TC_Veg_lm_int <- lm(Lipds ~ I(PRS_lipds_num * Strict) + Strict + PRS_lipds_num + Sex + Age + I(Age^2) + Array +
                                ac10003 + ac11001 + ac11002 +
                                ac11003 + ac11004 + ac11005 +
                                ac11006 + ac11007 + ac11008 +
                                ac11009 + ac11010 + ac11011 +
                                ac11011 + ac11012 + ac11013 +
                                ac11014 + ac11016 + ac11017 +
                                ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +
                                PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                                PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov +
                                Physical_activity_cov, data = bd_pheno_strict)

         
          bd_pheno_strict_tem_no_NA <- na.omit(bd_pheno_strict_tem)
          
          a_Low=format(round(summary(TC_Low_lm)$coefficients[2,1],digit=3), nsmall = 3)
          b_Low=format(round(summary(TC_Low_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_Low_lm)$coefficients[2,2],digit=3), nsmall = 3)
          c_Low=format(round(summary(TC_Low_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_Low_lm)$coefficients[2,2],digit=3), nsmall = 3)
          pval_Low=ifelse(summary(TC_Low_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Low_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Low_lm)$coefficients[2,4],digit=3), nsmall = 3))
          new_low=paste(a_Low," (",b_Low,", ",c_Low,")",sep="")
          
          a_Intermediate=format(round(summary(TC_Intermediate_lm)$coefficients[2,1],digit=3), nsmall = 3)
          b_Intermediate=format(round(summary(TC_Intermediate_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_Intermediate_lm)$coefficients[2,2],digit=3), nsmall = 3)
          c_Intermediate=format(round(summary(TC_Intermediate_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_Intermediate_lm)$coefficients[2,2],digit=3), nsmall = 3)
          pval_Intermediate=ifelse(summary(TC_Intermediate_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Intermediate_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Intermediate_lm)$coefficients[2,4],digit=3), nsmall = 3))
          new_Intermediate=paste(a_Intermediate," (",b_Intermediate,", ",c_Intermediate,")",sep="")
          
          a_High=format(round(summary(TC_High_lm)$coefficients[2,1],digit=3), nsmall = 3)
          b_High=format(round(summary(TC_High_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_High_lm)$coefficients[2,2],digit=3), nsmall = 3)
          c_High=format(round(summary(TC_High_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_High_lm)$coefficients[2,2],digit=3), nsmall = 3)
          pval_High=ifelse(summary(TC_High_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_High_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_High_lm)$coefficients[2,4],digit=3), nsmall = 3))
          new_High=paste(a_High," (",b_High,", ",c_High,")",sep="")
          
          pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
          
          
          final_res=data.frame(paste("          ",popul),new_low,pval_Low,new_Intermediate,pval_Intermediate,
                               new_High,pval_High,pval_int)
          
          
          
          write.table(final_res, file= paste(Pathway_out,"strat1204.txt", sep=""), col.names = FALSE, append = TRUE,
                      row.names = F, quote = FALSE, na = "",sep='\t')
          
          bd_pheno_strict <- subset (bd_pheno_strict, select = -c(PRS_lipds,Lipds,un_sd_PRS_lipds))
        }
      }
    }
  }
}































############################################################################################################################# SELF 
bd_pheno_self_tem1 = bd_pheno_self
for (percent_num in c(33)) {
  final_res=data.frame(NA,"Low PGS",NA,"Intermediate PGS",NA,"High PGS",NA,NA)
  write.table(final_res, file= paste(Pathway_out,"strat1204self.txt", sep=""), col.names = FALSE, append = TRUE,
              row.names = F, quote = FALSE, na = "",sep='\t')
  final_res=data.frame(NA,"Beta (95% CI)","P-value","Beta (95% CI)","P-value",
                       "Beta (95% CI)","P-value","Pinteractionb")
  write.table(final_res, file= paste(Pathway_out,"strat1204self.txt", sep=""), col.names = FALSE, append = TRUE,
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
    write.table(final_res, file= paste(Pathway_out,"strat1204self.txt", sep=""), col.names = FALSE, append = TRUE,
                row.names = F, quote = FALSE, na = "",sep='\t')
    for (x in c("Graha","Graham","Willer")) {
      if (x == c("Graha")) {
        final_res=data.frame("     Graham SE et al. (excluded UKB data)",NA,NA,NA,NA,NA,NA,NA)
      } else if (x == c("Graham")) {
        final_res=data.frame("     Graham SE et al. (included UKB data)",NA,NA,NA,NA,NA,NA,NA)
      } else if (x == c("Willer")) {
        final_res=data.frame("     Willer CJ et al.",NA,NA,NA,NA,NA,NA,NA)
      } 
      write.table(final_res, file= paste(Pathway_out,"strat1204self.txt", sep=""), col.names = FALSE, append = TRUE,
                  row.names = F, quote = FALSE, na = "",sep='\t')
      for (popul in c("EUR")) {
        
        if (x == c("Graha")) {
          final_res=data.frame("Graham SE et al. (excluded UKB data)",NA,NA,NA,NA,NA,NA,NA)
          Lipids_name=paste(e, "_",x,"m_no_UKB_",popul,"_PRS", sep="")
        } else if (x == c("Graham")) {
          final_res=data.frame("Graham SE et al. (included UKB data)",NA,NA,NA,NA,NA,NA,NA)
          Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
        } else if (x == c("Willer")) {
          final_res=data.frame("Willer CJ et al.",NA,NA,NA,NA,NA,NA,NA)
          Lipids_name=paste(e, "_",x,"_",popul,"_PRS", sep="")
        } 
        
        bd_pheno_self=bd_pheno_self_tem1
        
        
        if (Lipids_name %in% c(colnames(bd_pheno_self))) {
          
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
          
          bd_pheno_self <- bd_pheno_self[is.na(bd_pheno_self$PRS_lipds)==F,]
          
          if (percent_num == c(33)) {
            bd_pheno_self$tertile <- ntile(bd_pheno_self$PRS_lipds, 3)
            
            bd_pheno_self$PRS_lipds_tem=ifelse(is.na(bd_pheno_self$tertile)==F & bd_pheno_self$tertile==1, "Low", NA)
            bd_pheno_self$PRS_lipds_tem=ifelse(is.na(bd_pheno_self$tertile)==F & bd_pheno_self$tertile==3, "High", bd_pheno_self$PRS_lipds_tem)
            bd_pheno_self$PRS_lipds_tem=ifelse(is.na(bd_pheno_self$PRS_lipds_tem)==T, "Intermediate", bd_pheno_self$PRS_lipds_tem)
            bd_pheno_self$PRS_lipds_tem <- as.factor(bd_pheno_self$PRS_lipds_tem)
            bd_pheno_self$PRS_lipds_tem<- factor(bd_pheno_self$PRS_lipds_tem, levels = c("Low","Intermediate","High"))
            
          } 
          
          bd_pheno_self$SelfID_factor <- as.factor(bd_pheno_self$SelfID)
          bd_pheno_self$SelfID_factor<- factor(bd_pheno_self$SelfID, levels = c(0,1))
          
          
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
                                                        Physical_activity_cov,SelfID,PRS_lipds_tem)
          
          TC_Low_lm <- lm(Lipds ~ SelfID + Sex + Age + I(Age^2) + Array + 
                            ac10003 + ac11001 + ac11002 + 
                            ac11003 + ac11004 + ac11005 + 
                            ac11006 + ac11007 + ac11008 + 
                            ac11009 + ac11010 + ac11011 + 
                            ac11011 + ac11012 + ac11013 + 
                            ac11014 + ac11016 + ac11017 + 
                            ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                            PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + 
                            PC19 + PC20 + BMI + Townsend +
                            statin_sum_standardized + 
                            Smoking_status_cov + 
                            Alcohol_status_cov + 
                            Physical_activity_cov, data = filter(bd_pheno_self_tem, PRS_lipds_tem == "Low"))
          TC_Intermediate_lm <- lm(Lipds ~ SelfID + Sex + Age + I(Age^2) + Array +
                                     ac10003 + ac11001 + ac11002 +
                                     ac11003 + ac11004 + ac11005 +
                                     ac11006 + ac11007 + ac11008 +
                                     ac11009 + ac11010 + ac11011 +
                                     ac11011 + ac11012 + ac11013 +
                                     ac11014 + ac11016 + ac11017 +
                                     ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +
                                     PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                                     PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov +
                                     Physical_activity_cov, data = filter(bd_pheno_self_tem, PRS_lipds_tem == "Intermediate"))
          TC_High_lm <- lm(Lipds ~ SelfID + Sex + Age + I(Age^2) + Array +
                             ac10003 + ac11001 + ac11002 +
                             ac11003 + ac11004 + ac11005 +
                             ac11006 + ac11007 + ac11008 +
                             ac11009 + ac11010 + ac11011 +
                             ac11011 + ac11012 + ac11013 +
                             ac11014 + ac11016 + ac11017 +
                             ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +
                             PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                             PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov +
                             Physical_activity_cov, data = filter(bd_pheno_self_tem, PRS_lipds_tem == "High"))
          if (percent_num == 33) {
           bd_pheno_self$PRS_lipds_num <- ifelse(bd_pheno_self$tertile == 1, 0,
                                                    ifelse(bd_pheno_self$tertile == 3, 1, NA))
            }
          
          TC_Veg_lm_int <- lm(Lipds ~ I(PRS_lipds_num * SelfID) + SelfID + PRS_lipds_num + Sex + Age + I(Age^2) + Array +
                                ac10003 + ac11001 + ac11002 +
                                ac11003 + ac11004 + ac11005 +
                                ac11006 + ac11007 + ac11008 +
                                ac11009 + ac11010 + ac11011 +
                                ac11011 + ac11012 + ac11013 +
                                ac11014 + ac11016 + ac11017 +
                                ac11018 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 +
                                PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +
                                PC19 + PC20 + BMI + Townsend + statin_sum_standardized + Smoking_status_cov + Alcohol_status_cov +
                                Physical_activity_cov, data = bd_pheno_self)
          
          bd_pheno_self_tem_no_NA <- na.omit(bd_pheno_self_tem)
          
          a_Low=format(round(summary(TC_Low_lm)$coefficients[2,1],digit=3), nsmall = 3)
          b_Low=format(round(summary(TC_Low_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_Low_lm)$coefficients[2,2],digit=3), nsmall = 3)
          c_Low=format(round(summary(TC_Low_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_Low_lm)$coefficients[2,2],digit=3), nsmall = 3)
          pval_Low=ifelse(summary(TC_Low_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Low_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Low_lm)$coefficients[2,4],digit=3), nsmall = 3))
          new_low=paste(a_Low," (",b_Low,", ",c_Low,")",sep="")
          
          a_Intermediate=format(round(summary(TC_Intermediate_lm)$coefficients[2,1],digit=3), nsmall = 3)
          b_Intermediate=format(round(summary(TC_Intermediate_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_Intermediate_lm)$coefficients[2,2],digit=3), nsmall = 3)
          c_Intermediate=format(round(summary(TC_Intermediate_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_Intermediate_lm)$coefficients[2,2],digit=3), nsmall = 3)
          pval_Intermediate=ifelse(summary(TC_Intermediate_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_Intermediate_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Intermediate_lm)$coefficients[2,4],digit=3), nsmall = 3))
          new_Intermediate=paste(a_Intermediate," (",b_Intermediate,", ",c_Intermediate,")",sep="")
          
          a_High=format(round(summary(TC_High_lm)$coefficients[2,1],digit=3), nsmall = 3)
          b_High=format(round(summary(TC_High_lm)$coefficients[2,1]+qnorm(0.025)*summary(TC_High_lm)$coefficients[2,2],digit=3), nsmall = 3)
          c_High=format(round(summary(TC_High_lm)$coefficients[2,1]+qnorm(0.975)*summary(TC_High_lm)$coefficients[2,2],digit=3), nsmall = 3)
          pval_High=ifelse(summary(TC_High_lm)$coefficients[2,4]<0.0005,formatC(summary(TC_High_lm)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_High_lm)$coefficients[2,4],digit=3), nsmall = 3))
          new_High=paste(a_High," (",b_High,", ",c_High,")",sep="")
          
          pval_int=ifelse(summary(TC_Veg_lm_int)$coefficients[2,4]<0.0005,formatC(summary(TC_Veg_lm_int)$coefficients[2,4], format = "e", digits = 2),format(round(summary(TC_Veg_lm_int)$coefficients[2,4],digit=3), nsmall = 3))
          
          
          final_res=data.frame(paste("          ",popul),new_low,pval_Low,new_Intermediate,pval_Intermediate,
                               new_High,pval_High,pval_int)
          
          
          write.table(final_res, file= paste(Pathway_out,"strat1204self.txt", sep=""), col.names = FALSE, append = TRUE,
                      row.names = F, quote = FALSE, na = "",sep='\t')
          
          bd_pheno_self <- subset (bd_pheno_self, select = -c(PRS_lipds,Lipds,un_sd_PRS_lipds))
        }
      }
    }
  }
}









