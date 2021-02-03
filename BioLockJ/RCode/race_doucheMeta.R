#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-08-20
#Description: 

#Libraries
library(reshape2)
library(vegan)
library(ggpubr)
library(ggrepel)

rm(list = ls())

# setwd("~/BioLockJ_pipelines/vaginalMicrobiome_2021Feb03/04_Douche/script/")

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")
douchePath <- file.path(dir(pipeRoot, full.names = TRUE, pattern = "_Douche"), "output/")
clusterPath <- file.path(dir(pipeRoot, full.names = TRUE, pattern = "Cluster"), "output/")

meta2=read.csv(file=paste0(douchePath, "metadata2.csv"),row.names=1)

col_d=c("darkblue","steelblue","darkred","tomato")
#metadata for race and douching
gplots=list()
gplots[[1]]=ggboxplot(meta2, x = "race_douche", y = "C_MOMEDU", color = "race_douche", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20)))+labs( x = "", y = "Years of Maternal Education")+rotate_x_text(45)
gplots[[2]]=ggboxplot(meta2, x = "race_douche", y = "C_MAGE24", color = "race_douche", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20)))+labs( x = "", y = "Maternal Age at 24 weeks' gestation")+rotate_x_text(45)
gplots[[3]]=ggboxplot(meta2, x = "race_douche", y = "C_POVERT", color = "race_douche", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "", y = "Percent of Poverty Level")+rotate_x_text(45)
gplots[[4]]=ggboxplot(meta2, x = "race_douche", y = "C_BMI", color = "race_douche", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "", y = "Pre-Pregnancy BMI")+rotate_x_text(45)
gplots[[5]]=ggboxplot(meta2, x = "race_douche", y = "LENEG_S", color = "race_douche", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "", y = "Life events: Sum of negative values")+rotate_x_text(45)
gplots[[6]]=ggboxplot(meta2, x = "race_douche", y = "CES_SCOR", color = "race_douche", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "", y = "CES-D summary score")+rotate_x_text(45)


## ** Figure S1
ggexport(
  plotlist = gplots, filename = paste0(output, "meta_cor_race_douche.pdf"),
  ncol = 2, nrow = 3,width=15,height=25
)

meta1=read.csv(file=paste0(clusterPath, "metadata1.csv"),row.names=1)

aggregate(meta1$C_MAGE24~meta1$race_case,FUN=mean)
aggregate(meta1$C_MAGE24~meta1$race_case,FUN=sd)
by(meta1$C_MAGE24,meta1$race_case,function(i){length(which(is.na(i)))})

aggregate(meta1$C_MOMEDU~meta1$race_case,FUN=mean)
aggregate(meta1$C_MOMEDU~meta1$race_case,FUN=sd)
by(meta1$C_MOMEDU,meta1$race_case,function(i){length(which(is.na(i)))})

aggregate(meta1$C_BMI~meta1$race_case,FUN=mean)
aggregate(meta1$C_BMI~meta1$race_case,FUN=sd)
by(meta1$C_BMI,meta1$race_case,function(i){length(which(is.na(i)))})

a=data.frame(table(paste(meta1$race_case,meta1$CIGS_46A)))
b=rep(table(meta1$race_case),each=3)
a$Perc=a$Freq/b

a=data.frame(table(paste(meta1$race_case,meta1$C_PARITY)))
b=c(rep(table(meta1$race_case)[1],3),rep(table(meta1$race_case)[2],3),rep(table(meta1$race_case)[3],2),rep(table(meta1$race_case)[4],3))
a$Perc=a$Freq/b

a=data.frame(table(paste(meta1$race_case,meta1$C_MARITA)))
b=c(rep(table(meta1$race_case)[1],5),rep(table(meta1$race_case)[2],6),rep(table(meta1$race_case)[3],4),rep(table(meta1$race_case)[4],4))
a$Perc=a$Freq/b

aggregate(meta1$C_POVERT~meta1$race_case,FUN=mean)
aggregate(meta1$C_POVERT~meta1$race_case,FUN=sd)
by(meta1$C_POVERT,meta1$race_case,function(i){length(which(is.na(i)))})

a=data.frame(table(paste(meta1$race_case,meta1$douche)))
b=rep(table(meta1$race_case),each=3)
a$Perc=a$Freq/b

aggregate(meta1$CES_SCOR~meta1$race_case,FUN=mean)
aggregate(meta1$CES_SCOR~meta1$race_case,FUN=sd)
by(meta1$CES_SCOR,meta1$race_case,function(i){length(which(is.na(i)))})

aggregate(meta1$LENEG_S~meta1$race_case,FUN=mean)
aggregate(meta1$LENEG_S~meta1$race_case,FUN=sd)
by(meta1$LENEG_S,meta1$race_case,function(i){length(which(is.na(i)))})
