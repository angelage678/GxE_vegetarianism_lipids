# generating supplemental plots for manuscript (veg-lipid effect plots) 

'%ni%' <- Negate('%in%')
library(plyr) 
library(dplyr)
library(tidyverse)
library(readr)
library(dplyr)
library(forestplot)
library("RColorBrewer")

library(ggplot2)
library(RColorBrewer)
library("ggsci")
library(Hmisc)

Pathway=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/")
Pathway_out=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_OR_plots_04012025/")

####################################################################################################
#SELF
####################################################################################################
########## Graham_no_UKB_EUR_PRS			
data <- read.csv(paste(Pathway,"final_res_self_noukb_simple.txt", sep=""),as.is=T, header=F, sep="\t")
names(data)[1] <- "Exposure"
names(data)[2] <- "Beta"
names(data)[3] <- "SE"
names(data)[4] <- "pval"
names(data)[5] <- "pval_int"

data1=data%>% add_row()
data1=data[13,]
data4=rbind(data1,data)
data=data4

data$a1=data$Beta
data$b1=data$Beta+(qnorm(0.025))*data$SE
data$c1=data$Beta+qnorm(0.975)*data$SE
data$a=format(round(data$Beta,digit=3), nsmall = 3)
data$b=format(round(data$Beta+qnorm(0.025)*data$SE,digit=3), nsmall = 3)
data$c=format(round(data$Beta+qnorm(0.975)*data$SE,digit=3), nsmall = 3)
data$pval=ifelse(data$pval<0.0005,formatC(data$pval, format = "e", digits = 3),format(round(data$pval,digit=3), nsmall = 3))
data$pval_int=ifelse(data$pval_int<0.0005,formatC(data$pval_int, format = "e", digits = 3),format(round(data$pval_int,digit=3), nsmall = 3))
data$new=paste(data$a," [",data$b,"-",data$c,"]",sep="")

data$new[2]=NA
data$new[5]=NA
data$new[8]=NA
data$new[11]=NA
data$Beta[2]=NA
data$SE[2]=NA
data$pval[2]=NA
data$Beta[5]=NA
data$SE[5]=NA
data$pval[5]=NA
data$Beta[8]=NA
data$SE[8]=NA
data$pval[8]=NA
data$Beta[11]=NA
data$SE[11]=NA
data$pval[11]=NA
data$Exposure[1]="Lipids"
data$new[1]="Beta [95% CI]"
data$pval[1]="P-value"
data$pval_int[1]="P interaction" 

base_data <- tibble::tibble(mean  = data$a1[2:13],
                            lower = data$b1[2:13],
                            upper = data$c1[2:13],
                            Exposure = data$Exposure[2:13],
                            pval_int= data$pval_int[2:13],
                            new = data$new[2:13])

pdf(paste0(Pathway_out,"Graham_no_UKB_EUR_PRS_self.OR_plots.0312.pdf"), pointsize =11,  width =6.8, height  = 3.4,family="Times", onefile=F)


base_data |>
  forestplot(labeltext = c(Exposure, new,pval_int),
             boxsize = .25, # We set the box size to better visualize the type
             xticks = c(0.20, 0.24,0.28,0.32),
             zero = NA,
             txt_gp = fpTxtGp(label=gpar(fontfamily="Times"),ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.9), cex = 0.9), 
             fn.ci_norm=matrix(c("fpDrawCircleCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI"), 
                               nrow = 13, ncol=1, byrow=T),
             is.summary=c(rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2),
                          rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2)),
  )|>
  fp_set_style(
    box = "#006e00",
    line = "#006e00",
    summary = "#006e00"
  ) |> 
  fp_add_header(Exposure = c("Lipids"),
                new = c("Beta [95% CI]"),
                pval_int = c("P interaction")) 

dev.off()
########## Graham_no_UKB_EUR_PRS			
data <- read.csv(paste(Pathway,"final_res_self_graham_simple.txt", sep=""),as.is=T, header=F, sep="\t")
names(data)[1] <- "Exposure"
names(data)[2] <- "Beta"
names(data)[3] <- "SE"
names(data)[4] <- "pval"
names(data)[5] <- "pval_int"

data1=data%>% add_row()
data1=data[13,]
data4=rbind(data1,data)
data=data4

data$a1=data$Beta
data$b1=data$Beta+(qnorm(0.025))*data$SE
data$c1=data$Beta+qnorm(0.975)*data$SE
data$a=format(round(data$Beta,digit=3), nsmall = 3)
data$b=format(round(data$Beta+qnorm(0.025)*data$SE,digit=3), nsmall = 3)
data$c=format(round(data$Beta+qnorm(0.975)*data$SE,digit=3), nsmall = 3)
data$pval=ifelse(data$pval<0.0005,formatC(data$pval, format = "e", digits = 3),format(round(data$pval,digit=3), nsmall = 3))
data$pval_int=ifelse(data$pval_int<0.0005,formatC(data$pval_int, format = "e", digits = 3),format(round(data$pval_int,digit=3), nsmall = 3))
data$new=paste(data$a," [",data$b,"-",data$c,"]",sep="")

data$new[2]=NA
data$new[5]=NA
data$new[8]=NA
data$new[11]=NA
data$Beta[2]=NA
data$SE[2]=NA
data$pval[2]=NA
data$Beta[5]=NA
data$SE[5]=NA
data$pval[5]=NA
data$Beta[8]=NA
data$SE[8]=NA
data$pval[8]=NA
data$Beta[11]=NA
data$SE[11]=NA
data$pval[11]=NA
data$Exposure[1]="Lipids"
data$new[1]="Beta [95% CI]"
data$pval[1]="P-value"
data$pval_int[1]="P interaction" 

base_data <- tibble::tibble(mean  = data$a1[2:13],
                            lower = data$b1[2:13],
                            upper = data$c1[2:13],
                            Exposure = data$Exposure[2:13],
                            pval_int= data$pval_int[2:13],
                            new = data$new[2:13])

pdf(paste0(Pathway_out,"Graham_EUR_PRS_self.OR_plots.0312.pdf"), pointsize =11,  width =6.8, height  = 3.4,family="Times", onefile=F)


base_data |>
  forestplot(labeltext = c(Exposure, new,pval_int),
             boxsize = .25, # We set the box size to better visualize the type
             xticks = c(0.22, 0.26,0.30,0.34),
             zero = NA,
             txt_gp = fpTxtGp(label=gpar(fontfamily="Times"),ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.9), cex = 0.9), 
             fn.ci_norm=matrix(c("fpDrawCircleCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI"), 
                               nrow = 13, ncol=1, byrow=T),
             is.summary=c(rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2),
                          rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2)),
  )|>
  fp_set_style(
    box = "#006e00",
    line = "#006e00",
    summary = "#006e00"
  ) |> 
  fp_add_header(Exposure = c("Lipids"),
                new = c("Beta [95% CI]"),
                pval_int = c("P interaction")) 

dev.off()
########## Graham_no_UKB_EUR_PRS			
data <- read.csv(paste(Pathway,"final_res_self_willer_simple.txt", sep=""),as.is=T, header=F, sep="\t")
names(data)[1] <- "Exposure"
names(data)[2] <- "Beta"
names(data)[3] <- "SE"
names(data)[4] <- "pval"
names(data)[5] <- "pval_int"

data1=data%>% add_row()
data1=data[13,]
data4=rbind(data1,data)
data=data4

data$a1=data$Beta
data$b1=data$Beta+(qnorm(0.025))*data$SE
data$c1=data$Beta+qnorm(0.975)*data$SE
data$a=format(round(data$Beta,digit=3), nsmall = 3)
data$b=format(round(data$Beta+qnorm(0.025)*data$SE,digit=3), nsmall = 3)
data$c=format(round(data$Beta+qnorm(0.975)*data$SE,digit=3), nsmall = 3)
data$pval=ifelse(data$pval<0.0005,formatC(data$pval, format = "e", digits = 3),format(round(data$pval,digit=3), nsmall = 3))
data$pval_int=ifelse(data$pval_int<0.0005,formatC(data$pval_int, format = "e", digits = 3),format(round(data$pval_int,digit=3), nsmall = 3))
data$new=paste(data$a," [",data$b,"-",data$c,"]",sep="")

data$new[2]=NA
data$new[5]=NA
data$new[8]=NA
data$new[11]=NA
data$Beta[2]=NA
data$SE[2]=NA
data$pval[2]=NA
data$Beta[5]=NA
data$SE[5]=NA
data$pval[5]=NA
data$Beta[8]=NA
data$SE[8]=NA
data$pval[8]=NA
data$Beta[11]=NA
data$SE[11]=NA
data$pval[11]=NA
data$Exposure[1]="Lipids"
data$new[1]="Beta [95% CI]"
data$pval[1]="P-value"
data$pval_int[1]="P interaction" 

base_data <- tibble::tibble(mean  = data$a1[2:13],
                            lower = data$b1[2:13],
                            upper = data$c1[2:13],
                            Exposure = data$Exposure[2:13],
                            pval_int= data$pval_int[2:13],
                            new = data$new[2:13])

pdf(paste0(Pathway_out,"Willer_EUR_PRS_self.OR_plots.0312.pdf"), pointsize =11,  width =6.8, height  = 3.4,family="Times", onefile=F)


base_data |>
  forestplot(labeltext = c(Exposure, new,pval_int),
             boxsize = .25, # We set the box size to better visualize the type
             xticks = c(0.14, 0.20, 0.26,0.32),
             zero = NA,
             txt_gp = fpTxtGp(label=gpar(fontfamily="Times"),ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.9), cex = 0.9), 
             fn.ci_norm=matrix(c("fpDrawCircleCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI"), 
                               nrow = 13, ncol=1, byrow=T),
             is.summary=c(rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2),
                          rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2)),
  )|>
  fp_set_style(
    box = "#006e00",
    line = "#006e00",
    summary = "#006e00"
  ) |> 
  fp_add_header(Exposure = c("Lipids"),
                new = c("Beta [95% CI]"),
                pval_int = c("P interaction")) 

dev.off()



























####################################################################################################
#STRICT
####################################################################################################
########## Graham_no_UKB_EUR_PRS			
data <- read.csv(paste(Pathway,"final_res_noukb_simple.txt", sep=""),as.is=T, header=F, sep="\t")
names(data)[1] <- "Exposure"
names(data)[2] <- "Beta"
names(data)[3] <- "SE"
names(data)[4] <- "pval"
names(data)[5] <- "pval_int"

data1=data%>% add_row()
data1=data[13,]
data4=rbind(data1,data)
data=data4

data$a1=data$Beta
data$b1=data$Beta+(qnorm(0.025))*data$SE
data$c1=data$Beta+qnorm(0.975)*data$SE
data$a=format(round(data$Beta,digit=3), nsmall = 3)
data$b=format(round(data$Beta+qnorm(0.025)*data$SE,digit=3), nsmall = 3)
data$c=format(round(data$Beta+qnorm(0.975)*data$SE,digit=3), nsmall = 3)
data$pval=ifelse(data$pval<0.0005,formatC(data$pval, format = "e", digits = 3),format(round(data$pval,digit=3), nsmall = 3))
data$pval_int=ifelse(data$pval_int<0.0005,formatC(data$pval_int, format = "e", digits = 3),format(round(data$pval_int,digit=3), nsmall = 3))
data$new=paste(data$a," [",data$b,"-",data$c,"]",sep="")

data$new[2]=NA
data$new[5]=NA
data$new[8]=NA
data$new[11]=NA
data$Beta[2]=NA
data$SE[2]=NA
data$pval[2]=NA
data$Beta[5]=NA
data$SE[5]=NA
data$pval[5]=NA
data$Beta[8]=NA
data$SE[8]=NA
data$pval[8]=NA
data$Beta[11]=NA
data$SE[11]=NA
data$pval[11]=NA
data$Exposure[1]="Lipids"
data$new[1]="Beta [95% CI]"
data$pval[1]="P-value"
data$pval_int[1]="P interaction" 

base_data <- tibble::tibble(mean  = data$a1[2:13],
                            lower = data$b1[2:13],
                            upper = data$c1[2:13],
                            Exposure = data$Exposure[2:13],
                            pval_int= data$pval_int[2:13],
                            new = data$new[2:13])

pdf(paste0(Pathway_out,"Graham_no_UKB_EUR_PRS_strict.OR_plots.0312.pdf"), pointsize =11,  width =6.8, height  = 3.4,family="Times", onefile=F)


base_data |>
  forestplot(labeltext = c(Exposure, new,pval_int),
             boxsize = .25, # We set the box size to better visualize the type
             xticks = c(0.18, 0.24,0.30,0.36),
             zero = NA,
             txt_gp = fpTxtGp(label=gpar(fontfamily="Times"),ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.9), cex = 0.9), 
             fn.ci_norm=matrix(c("fpDrawCircleCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI"), 
                               nrow = 13, ncol=1, byrow=T),
             is.summary=c(rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2),
                          rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2)),
  )|>
  fp_set_style(
    box = "#006e00",
    line = "#006e00",
    summary = "#006e00"
  ) |> 
  fp_add_header(Exposure = c("Lipids"),
                new = c("Beta [95% CI]"),
                pval_int = c("P interaction")) 

dev.off()
########## Graham_no_UKB_EUR_PRS			
data <- read.csv(paste(Pathway,"final_res_graham_simple.txt", sep=""),as.is=T, header=F, sep="\t")
names(data)[1] <- "Exposure"
names(data)[2] <- "Beta"
names(data)[3] <- "SE"
names(data)[4] <- "pval"
names(data)[5] <- "pval_int"

data1=data%>% add_row()
data1=data[13,]
data4=rbind(data1,data)
data=data4

data$a1=data$Beta
data$b1=data$Beta+(qnorm(0.025))*data$SE
data$c1=data$Beta+qnorm(0.975)*data$SE
data$a=format(round(data$Beta,digit=3), nsmall = 3)
data$b=format(round(data$Beta+qnorm(0.025)*data$SE,digit=3), nsmall = 3)
data$c=format(round(data$Beta+qnorm(0.975)*data$SE,digit=3), nsmall = 3)
data$pval=ifelse(data$pval<0.0005,formatC(data$pval, format = "e", digits = 3),format(round(data$pval,digit=3), nsmall = 3))
data$pval_int=ifelse(data$pval_int<0.0005,formatC(data$pval_int, format = "e", digits = 3),format(round(data$pval_int,digit=3), nsmall = 3))
data$new=paste(data$a," [",data$b,"-",data$c,"]",sep="")

data$new[2]=NA
data$new[5]=NA
data$new[8]=NA
data$new[11]=NA
data$Beta[2]=NA
data$SE[2]=NA
data$pval[2]=NA
data$Beta[5]=NA
data$SE[5]=NA
data$pval[5]=NA
data$Beta[8]=NA
data$SE[8]=NA
data$pval[8]=NA
data$Beta[11]=NA
data$SE[11]=NA
data$pval[11]=NA
data$Exposure[1]="Lipids"
data$new[1]="Beta [95% CI]"
data$pval[1]="P-value"
data$pval_int[1]="P interaction" 

base_data <- tibble::tibble(mean  = data$a1[2:13],
                            lower = data$b1[2:13],
                            upper = data$c1[2:13],
                            Exposure = data$Exposure[2:13],
                            pval_int= data$pval_int[2:13],
                            new = data$new[2:13])

pdf(paste0(Pathway_out,"Graham_EUR_PRS_strict.OR_plots.0312.pdf"), pointsize =11,  width =6.8, height  = 3.4,family="Times", onefile=F)


base_data |>
  forestplot(labeltext = c(Exposure, new,pval_int),
             boxsize = .25, # We set the box size to better visualize the type
             xticks = c(0.20, 0.26,0.32,0.38),
             zero = NA,
             txt_gp = fpTxtGp(label=gpar(fontfamily="Times"),ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.9), cex = 0.9), 
             fn.ci_norm=matrix(c("fpDrawCircleCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI"), 
                               nrow = 13, ncol=1, byrow=T),
             is.summary=c(rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2),
                          rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2)),
  )|>
  fp_set_style(
    box = "#006e00",
    line = "#006e00",
    summary = "#006e00"
  ) |> 
  fp_add_header(Exposure = c("Lipids"),
                new = c("Beta [95% CI]"),
                pval_int = c("P interaction")) 

dev.off()
########## Graham_no_UKB_EUR_PRS			
data <- read.csv(paste(Pathway,"final_res_willer_simple.txt", sep=""),as.is=T, header=F, sep="\t")
names(data)[1] <- "Exposure"
names(data)[2] <- "Beta"
names(data)[3] <- "SE"
names(data)[4] <- "pval"
names(data)[5] <- "pval_int"

data1=data%>% add_row()
data1=data[13,]
data4=rbind(data1,data)
data=data4

data$a1=data$Beta
data$b1=data$Beta+(qnorm(0.025))*data$SE
data$c1=data$Beta+qnorm(0.975)*data$SE
data$a=format(round(data$Beta,digit=3), nsmall = 3)
data$b=format(round(data$Beta+qnorm(0.025)*data$SE,digit=3), nsmall = 3)
data$c=format(round(data$Beta+qnorm(0.975)*data$SE,digit=3), nsmall = 3)
data$pval=ifelse(data$pval<0.0005,formatC(data$pval, format = "e", digits = 3),format(round(data$pval,digit=3), nsmall = 3))
data$pval_int=ifelse(data$pval_int<0.0005,formatC(data$pval_int, format = "e", digits = 3),format(round(data$pval_int,digit=3), nsmall = 3))
data$new=paste(data$a," [",data$b,"-",data$c,"]",sep="")

data$new[2]=NA
data$new[5]=NA
data$new[8]=NA
data$new[11]=NA
data$Beta[2]=NA
data$SE[2]=NA
data$pval[2]=NA
data$Beta[5]=NA
data$SE[5]=NA
data$pval[5]=NA
data$Beta[8]=NA
data$SE[8]=NA
data$pval[8]=NA
data$Beta[11]=NA
data$SE[11]=NA
data$pval[11]=NA
data$Exposure[1]="Lipids"
data$new[1]="Beta [95% CI]"
data$pval[1]="P-value"
data$pval_int[1]="P interaction" 

base_data <- tibble::tibble(mean  = data$a1[2:13],
                            lower = data$b1[2:13],
                            upper = data$c1[2:13],
                            Exposure = data$Exposure[2:13],
                            pval_int= data$pval_int[2:13],
                            new = data$new[2:13])

pdf(paste0(Pathway_out,"Willer_EUR_PRS_strict.OR_plots.0312.pdf"), pointsize =11,  width =6.8, height  = 3.4,family="Times", onefile=F)


base_data |>
  forestplot(labeltext = c(Exposure, new,pval_int),
             boxsize = .25, # We set the box size to better visualize the type
             xticks = c(0.08, 0.16, 0.24,0.32),
             zero = NA,
             txt_gp = fpTxtGp(label=gpar(fontfamily="Times"),ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.9), cex = 0.9),
             fn.ci_norm=matrix(c("fpDrawCircleCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI"), 
                               nrow = 13, ncol=1, byrow=T),
             is.summary=c(rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2),
                          rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2)),
  )|>
  fp_set_style(
    box = "#006e00",
    line = "#006e00",
    summary = "#006e00"
  ) |> 
  fp_add_header(Exposure = c("Lipids"),
                new = c("Beta [95% CI]"),
                pval_int = c("P interaction")) 

dev.off()
