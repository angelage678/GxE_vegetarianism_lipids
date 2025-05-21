# generating forest plots for manuscript

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

Pathway=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_analyses/")
Pathway_out=c("/scratch/ag42790/vegetarianPRSLipids/vegetarian_files/new_OR_plots_04012025/")

#EUR		
data <- read.csv(paste(Pathway,"new_forest_plot_EUR.txt", sep=""),as.is=T, header=F, sep="\t")
names(data)[1] <- "Exposure"
names(data)[2] <- "Beta"
names(data)[3] <- "SE"
names(data)[4] <- "pval"

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

base_data <- tibble::tibble(mean  = data$a1[2:13],
                            lower = data$b1[2:13],
                            upper = data$c1[2:13],
                            Exposure = data$Exposure[2:13],
                            new = data$new[2:13])

pdf(paste0(Pathway_out,"EUR_plot.pdf"), pointsize =13,  width =6.8, height  = 3.4,family="Times", onefile=F)


base_data |>
  forestplot(labeltext = c(Exposure, new),
             boxsize = .25, # We set the box size to better visualize the type
             #line.margin = .1, # We need to add this to avoid crowding
             #clip = c(0.22, 0.34),
             xticks = c(-0.3, -0.2, -0.1, 0 ,0.1, 0.2, 0.3),
             #colgap=unit(4, "mm"),
             #graph.pos = 2,
             zero = NA,
             #zero.lwd=gpar(lwd=1, col="white"),
             txt_gp = fpTxtGp(label=gpar(fontfamily="Times"),ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.9), cex = 0.9), #文本大小设置
             fn.ci_norm=matrix(c("fpDrawCircleCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI"), 
                               nrow = 13, ncol=1, byrow=T),
             is.summary=c(rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2),
                          rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2)),
             #xlog = F
  )|>
  #fp_add_lines(h_2 = gpar(lwd = 1, col = "#000044")) |> 
  fp_set_style(
    #box=c(rep(pal_npg("nrc", alpha = 1)(9)[1],3),pal_npg("nrc", alpha = 1)(9)[2],
    #rep(pal_npg("nrc", alpha = 1)(9)[1],2),pal_npg("nrc", alpha = 1)(9)[2],
    #rep(pal_npg("nrc", alpha = 1)(9)[1],2),pal_npg("nrc", alpha = 1)(9)[2],
    #rep(pal_npg("nrc", alpha = 1)(9)[1],2),pal_npg("nrc", alpha = 1)(9)[2]),
    
    box = "#006e00",
    line = "#006e00",
    summary = "#006e00"
  ) |> 
  fp_add_header(Exposure = c("Lipids"),
                new = c("Beta [95% CI]")) 

dev.off()



#CSA		
data <- read.csv(paste(Pathway,"new_forest_plot_CSA.txt", sep=""),as.is=T, header=F, sep="\t")
names(data)[1] <- "Exposure"
names(data)[2] <- "Beta"
names(data)[3] <- "SE"
names(data)[4] <- "pval"

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

base_data <- tibble::tibble(mean  = data$a1[2:13],
                            lower = data$b1[2:13],
                            upper = data$c1[2:13],
                            Exposure = data$Exposure[2:13],
                            new = data$new[2:13])

pdf(paste0(Pathway_out,"CSA_plot.pdf"), pointsize =13,  width =6.8, height  = 3.4,family="Times", onefile=F)


base_data |>
  forestplot(labeltext = c(Exposure, new),
             boxsize = .25, # We set the box size to better visualize the type
             #line.margin = .1, # We need to add this to avoid crowding
             #clip = c(0.22, 0.34),
             xticks = c(-0.5, -0.4, -0.3, -0.2, -0.1, 0 ,0.1, 0.2, 0.3, 0.4, 0.5),
             #colgap=unit(4, "mm"),
             #graph.pos = 2,
             zero = NA,
             #zero.lwd=gpar(lwd=1, col="white"),
             txt_gp = fpTxtGp(label=gpar(fontfamily="Times"),ticks = gpar(cex = 0.9), xlab = gpar(cex = 0.9), cex = 0.9), #文本大小设置
             fn.ci_norm=matrix(c("fpDrawCircleCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI",  "fpDrawNormalCI","fpDrawCircleCI","fpDrawDiamondCI"), 
                               nrow = 13, ncol=1, byrow=T),
             is.summary=c(rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2),
                          rep(TRUE,1),rep(FALSE,2),rep(TRUE,1),rep(FALSE,2)),
             #xlog = F
  )|>
  #fp_add_lines(h_2 = gpar(lwd = 1, col = "#000044")) |> 
  fp_set_style(
    #box=c(rep(pal_npg("nrc", alpha = 1)(9)[1],3),pal_npg("nrc", alpha = 1)(9)[2],
    #rep(pal_npg("nrc", alpha = 1)(9)[1],2),pal_npg("nrc", alpha = 1)(9)[2],
    #rep(pal_npg("nrc", alpha = 1)(9)[1],2),pal_npg("nrc", alpha = 1)(9)[2],
    #rep(pal_npg("nrc", alpha = 1)(9)[1],2),pal_npg("nrc", alpha = 1)(9)[2]),
    
    box = "#006e00",
    line = "#006e00",
    summary = "#006e00"
  ) |> 
  fp_add_header(Exposure = c("Lipids"),
                new = c("Beta [95% CI]")) 

dev.off()
