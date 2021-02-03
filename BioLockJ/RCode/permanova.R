#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-08-20
#Description: 

#Libraries
library(reshape2)
library(vegan)
library(ggpubr)

rm(list = ls())

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")
racecasePath <- file.path(dir(pipeRoot, full.names = TRUE, pattern = "_RaceCase"), "output/")

meta1=read.csv(file=paste0(racecasePath, "metadata.csv"),row.names=1)
stir_norm1=read.csv(file=paste0(racecasePath, "stir_norm1.csv"),row.names=1)

#permanova
meta_t=c("C_MRACE","C_MARITA","edu_cat", "C_MAGE24","C_POVERT","C_PARITY", "C_BMI","douche","LENEG_S","ORRNEG_S","CES_SCOR","CIGS_46A","C_CASE2")
meta1_sn=match(meta_t,colnames(meta1))
ado_R2_vag=matrix(nrow=length(meta1_sn),ncol=2)
n=1
for (i in meta1_sn){
  tab1=na.omit(data.frame(cbind(stir_norm1,meta1[,i]), stringsAsFactors=FALSE))#remove na rows
  tab1[,1:337]=apply(tab1[,1:337],2,as.character)
  tab1[,1:337]=apply(tab1[,1:337],2,as.numeric)
  dim(tab1) 
  if(length(table(tab1[,338]))>15){
    #BMI is significant when log transformed (more changes in rare species)
    fit=adonis(tab1[,1:337] ~ tab1[,338], data=tab1, permutations=999)
  }else{
    fit=adonis(tab1[,1:337] ~ factor(tab1[,338]), data=tab1, permutations=999)
  }
  ado_R2_vag[n,1]=fit$aov.tab[1,5]
  ado_R2_vag[n,2]=fit$aov.tab[1,6]
  n=n+1
}
ado_R2_vag=na.omit(ado_R2_vag)
rownames(ado_R2_vag)=meta_t
colnames(ado_R2_vag)=c("R2","P")
ado_R2_vag=data.frame(ado_R2_vag)
ado_R2_vag$FDR=p.adjust(ado_R2_vag$P,method="fdr")
write.csv(ado_R2_vag,file=paste0(output, "ado_R2_vag_tax_new.csv"))
#ado_R2_vag=read.csv(file="ado_R2_vag_tax_new.csv",row.names=1)
ado_R2_vag=ado_R2_vag[ado_R2_vag$FDR<0.1,]
ado_R2_vag=ado_R2_vag[order(ado_R2_vag[,1],decreasing = T),]

bar_lab=rownames(ado_R2_vag)
bar_lab=c("Poverty level","Education","Marital status" ,"Age","Douching","Race","Parity","Depression","Negative life events","SPTB","pre-pregnancy BMI")
#col_bar=c("pink","yellow","lightgreen","lightblue")[factor(rev(c(1,1,2,2,3,3,2,3,3,2,1,3,2,4,3,3,3,3,3,4,3,4,rep(3,3))))]


## ** Figure 1c
pdf(paste0(output, "adonis_r2_bar.pdf"),width=6,height=8)
par(mar=c(5,15,5,5),mfrow=c(1,1))
barplot(rev(as.numeric(ado_R2_vag[,1])),xlim=c(0,0.1),col="palegreen",names.arg=rev(bar_lab[1:nrow(ado_R2_vag)]),las=2,main="Overall",horiz=T,cex.names=1.5,cex.axis=1.5)
dev.off()
