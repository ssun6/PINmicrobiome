#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-08-20
#Description: 

#Libraries
library(reshape2)
library(vegan)
library(ggpubr)

rm(list = ls())

# setwd("~/Documents/BioLockJ_pipelines/vaginalMicrobiome_2020Dec08/02_MetaRaceCase/script/")

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")
racecasePath <- file.path(dir(pipeRoot, full.names = TRUE, pattern = "_RaceCase"), "output/")

meta1=read.csv(file=paste0(racecasePath, "metadata.csv"),row.names=1)

#metadata for race and case
col_d=c("darkblue","steelblue","darkred","tomato")
gplots=list()

gplots[[1]]=ggboxplot(meta1, x = "race_case", y = "C_MOMEDU", color = "race_case", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20)))+labs( x = "", y = "Years of Maternal Education")+rotate_x_text(45)
gplots[[2]]=ggboxplot(meta1, x = "race_case", y = "C_MAGE24", color = "race_case", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20)))+labs( x = "", y = "Maternal Age at 24 weeks' gestation")+rotate_x_text(45)
gplots[[3]]=ggboxplot(meta1, x = "race_case", y = "C_POVERT", color = "race_case", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "", y = "Percent of Poverty Level")+rotate_x_text(45)
gplots[[4]]=ggboxplot(meta1, x = "race_case", y = "C_BMI", color = "race_case", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "", y = "Pre-Pregnancy BMI")+rotate_x_text(45)
gplots[[5]]=ggboxplot(meta1, x = "race_case", y = "LENEG_S", color = "race_case", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "", y = "Life events: Sum of negative values")+rotate_x_text(45)
gplots[[6]]=ggboxplot(meta1, x = "race_case", y = "CES_SCOR", color = "race_case", add = "jitter", show.legend=FALSE,palette =col_d,
                      ggtheme = theme_bw()+theme(legend.position="none",panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                                 panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                                 axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "", y = "CES-D summary score")+rotate_x_text(45)

## ** Figure S2
ggexport(
  plotlist = gplots, filename = paste0(output, "meta_cor_race_case.pdf"),
  ncol = 2, nrow = 3,width=15,height=25
)

