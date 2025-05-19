# Adding PRS data

'%ni%' <- Negate('%in%')
library(dplyr)

Pathway=c("/scratch/ag42790/vegetarianPRSLipids/All_PRS/PRS/")
Pathway_out=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/")

load(paste("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/","veg_merged_data.RData",sep = ""))
veg_bd_pheno <- veg_merged_data

################################### Willer CJ et al. Nat. Genet. 2013
# TC EUR
# LDL EUR
# HDLcholesterol HDL EUR
# Triglycerides TGs EUR

for (x in c("Willer")) {
  for (popul in c("EUR")) {
    for (e in c("TC","LDL","HDL","TGs")) {
      dat=read.table(paste(Pathway,x,"_chr","1","_",e,"_",popul,".profile", sep=""),header=T,sep = "")
      dat<-dat%>%select(IID, SCORESUM)
      colnames(dat)[2]=c("SCORE_1")
      for (chr in (c("2" ,"3", "4" ,"5" ,"6" ,"7", "8" ,"9" ,"10" ,"11", "12", "13", "14", "15", "16" ,"17", "18" ,"19", "20" ,"21", "22", "X" ,"XY"))) {
        Position=try(read.table(paste(Pathway,x,"_chr",chr,"_",e,"_",popul,".profile", sep=""),header=T,sep = ""), silent=TRUE)
        if (length(colnames(Position))>0) {
          Position<-Position%>%select(IID, SCORESUM)
          colnames(Position)[2]=paste("SCORE_",chr, sep="")
          dat<- dat %>% left_join(Position, by= "IID")
        }
      }
      
      dat$Total_PRS <- rowSums(dat, na.rm=TRUE)-dat$IID
      colnames(dat)[1]=c("FID")
      dat<-dat%>%select(FID, Total_PRS)
      print(summary(dat$Total_PRS))
      colnames(dat)[length(colnames(dat))]=paste(e, "_",x,"_",popul,"_PRS", sep="")
      print(colnames(dat)[length(colnames(dat))])
      
      veg_bd_pheno <-veg_bd_pheno %>% left_join(dat, by= "FID")
    }
  }
}

################################### Graham SE et al. Nature (2021)
############# EUR
############# AFR
############# EAS
############# SAS
for (x in c("Graham")) {
  for (popul in c("EUR","AFR","EAS","SAS")) {
    for (e in c("TC","LDL","HDL","TGs")) {
      dat=read.table(paste(Pathway,x,"_chr","1","_",e,"_",popul,".profile", sep=""),header=T,sep = "")
      dat<-dat%>%select(IID, SCORESUM)
      colnames(dat)[2]=c("SCORE_1")
      for (chr in (c("2" ,"3", "4" ,"5" ,"6" ,"7", "8" ,"9" ,"10" ,"11", "12", "13", "14", "15", "16" ,"17", "18" ,"19", "20" ,"21", "22", "X" ,"XY"))) {
        Position=try(read.table(paste(Pathway,x,"_chr",chr,"_",e,"_",popul,".profile", sep=""),header=T,sep = ""), silent=TRUE)
        if (length(colnames(Position))>0) {
          Position<-Position%>%select(IID, SCORESUM)
          colnames(Position)[2]=paste("SCORE_",chr, sep="")
          dat<- dat %>% left_join(Position, by= "IID")
        }
      }
      
      dat$Total_PRS <- rowSums(dat, na.rm=TRUE)-dat$IID
      colnames(dat)[1]=c("FID")
      dat<-dat%>%select(FID, Total_PRS)
      print(summary(dat$Total_PRS))
      colnames(dat)[length(colnames(dat))]=paste(e, "_",x,"_",popul,"_PRS", sep="")
      print(colnames(dat)[length(colnames(dat))])
      
      veg_bd_pheno <-veg_bd_pheno %>% left_join(dat, by= "FID")
    }
  }
}

################################### ###########       No UKB   Graham SE et al. Nature (2021)
############# EUR
############# AFR
############# SAS
for (x in c("Graham")) {
  for (popul in c("EUR","AFR","SAS")) {
    for (e in c("TC","LDL","HDL","TGs")) {
      dat=try(read.table(paste(Pathway,x,"_no_UKB_chr","1","_",e,"_",popul,".profile", sep=""),header=T,sep = ""), silent=TRUE)
      dat<-try(dat%>%select(IID, SCORESUM), silent=TRUE)
      if (length(colnames(dat))>0) {
        colnames(dat)[2]=c("SCORE_1")
        
        for (chr in (c("2" ,"3", "4" ,"5" ,"6" ,"7", "8" ,"9" ,"10" ,"11", "12", "13", "14", "15", "16" ,"17", "18" ,"19", "20" ,"21", "22", "X" ,"XY"))) {
          Position=try(read.table(paste(Pathway,x,"_no_UKB_chr",chr,"_",e,"_",popul,".profile", sep=""),header=T,sep = ""), silent=TRUE)
          if (length(colnames(Position))>0) {
            Position<-try(Position%>%select(IID, SCORESUM), silent=TRUE)
            colnames(Position)[2]=paste("SCORE_",chr, sep="")
            dat<- try(dat %>% left_join(Position, by= "IID"), silent=TRUE)
          }
        }
      } else {
        dat=try(read.table(paste(Pathway,x,"_no_UKB_chr","2","_",e,"_",popul,".profile", sep=""),header=T,sep = ""), silent=TRUE)
        dat<-try(dat%>%select(IID, SCORESUM), silent=TRUE)
        colnames(dat)[2]=c("SCORE_2")
        for (chr in (c("3", "4" ,"5" ,"6" ,"7", "8" ,"9" ,"10" ,"11", "12", "13", "14", "15", "16" ,"17", "18" ,"19", "20" ,"21", "22", "X" ,"XY"))) {
          Position=try(read.table(paste(Pathway,x,"_no_UKB_chr",chr,"_",e,"_",popul,".profile", sep=""),header=T,sep = ""), silent=TRUE)
          if (length(colnames(Position))>0) {
            Position<-try(Position%>%select(IID, SCORESUM), silent=TRUE)
            colnames(Position)[2]=paste("SCORE_",chr, sep="")
            dat<- try(dat %>% left_join(Position, by= "IID"), silent=TRUE)
          }
        }
      }
      
      dat$Total_PRS <- rowSums(dat, na.rm=TRUE)-dat$IID
      colnames(dat)[1]=c("FID")
      dat<-dat%>%select(FID, Total_PRS)
      print(summary(dat$Total_PRS))
      colnames(dat)[length(colnames(dat))]=paste(e, "_",x,"_no_UKB_",popul,"_PRS", sep="")
      print(colnames(dat)[length(colnames(dat))])
      
      veg_bd_pheno <-veg_bd_pheno %>% left_join(dat, by= "FID")
    }
  }
}

save(veg_bd_pheno, file = "veg_bd_pheno.RData")

