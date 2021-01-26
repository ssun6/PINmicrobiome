library(vegan)
library(ggplot2)
library(tidyr)
library(ggthemes)
library(reshape2)
library(gridExtra)
library(RGraphics)
library(pheatmap)
library(ggpubr)
library(ggrepel)
library(farver)

rm(list=ls())
setwd("/Users/shansun/Google\ Drive/engel/test/test/")
stir_count_melt=read.table(file="Engel_16S_stirrups_summary_97_070819.txt",sep="|",quote="")
stir_count_melt1=stir_count_melt[,c(1,2,4)]
stir_count_melt1[,3]=as.numeric(as.character(stir_count_melt1[,3]))
stir_count2=acast(stir_count_melt1, V1~ V2,sum)
stir_count2[is.na(stir_count2)]=0
stir_count2=stir_count2[1:1068,]

stir_count_melt_p1=stir_count_melt[stir_count_melt[,3]=="AT",c(1,2,5)]
stir_count_p=acast(stir_count_melt_p1, V1~ V2)
stir_count_p[is.na(stir_count_p)]=0
stir_count_p=stir_count_p[1:1068,]

stir_count_p=cbind(stir_count_p,100-rowSums(stir_count_p))
colnames(stir_count_p)[337]="other"
stir_count_p=stir_count_p[which(rowSums(stir_count2)>1000),]
stir_count_p=stir_count_p
stir_count=stir_count_p*mean(rowSums(stir_count2))/100

stir_count2=stir_count2[which(rowSums(stir_count2)>1000),]
stir_count2=stir_count2/rowSums(stir_count2)*mean(rowSums(stir_count2))

rownames(stir_count)=sapply(strsplit(rownames(stir_count),"\\-"), "[[", 1)
which(table(sapply(strsplit(rownames(stir_count),"\\-"), "[[", 1))>1)

#1152 1068
#1098 1061
#19 duplicates


meta_or=read.csv(file="PIN0108_MicrobCov_NoSNP_20190607.csv",row.names=1,header=T,quote="")
#new
meta_new=read.csv(file="PIN0113_MbCov_AllDiet_20200706.csv",row.names=1,header=T,quote="")
meta=cbind(meta_or,meta_new[,setdiff(colnames(meta_new),colnames(meta_or))])

meta$rc=paste(meta$C_MRACE,meta$C_CASE2,sep="_")
meta$rcn=paste(meta$C_MRACE,meta$C_CASCN2,sep="_")
meta$edu_cat=as.numeric(meta$C_MOMEDU)
meta$edu_cat[meta$edu_cat>16]="16+"
meta$edu_cat[meta$edu_cat%in%c(6:11)]="<12"
meta$edu_cat[meta$edu_cat%in%c(12:15)]="12-15"
meta$douche=meta$MQE17C
meta$douche[meta$douche>0]=1
meta$rd=paste(meta$C_MRACE,meta$douche,sep="_")
meta$rdc=paste(meta$C_MRACE,meta$douche,meta$C_CASE2,sep="_")
meta$rdc=factor(meta$rdc,levels=c("1_0_0","1_0_1","1_1_0","1_1_1","2_0_0" ,"2_0_1" , "2_1_0", "2_1_1"))
meta$race1=c("White","Black")[factor(meta$C_MRACE)]
meta$case1=c("Term","SPTB")[factor(meta$C_CASE2)]
meta$douche=factor(c("No douching","Douching")[factor(meta$douche)],levels=c("No douching","Douching"))
meta$C_PARITY[meta$C_PARITY>0]=1

stir_norm=stir_count/rowSums(stir_count)*mean(rowSums(stir_count))
meta=meta[intersect(rownames(stir_norm),rownames(meta)),]
stir_norm1=stir_norm[intersect(rownames(stir_norm),rownames(meta)),]
stir_norm1=stir_norm1[which(meta$C_MRACE%in%c(1,2)&meta$WHYPTB2%in%c(0,1,2)),]
meta1=meta[which(meta$C_MRACE%in%c(1,2)&meta$WHYPTB2%in%c(0,1,2)),]
meta1$Case=c("T","S")[factor(meta1$C_CASE2)]
meta1$Race=c("W","B")[factor(meta1$C_MRACE)]
meta1$Douche=factor(c("N","D")[factor(meta1$douche)],levels=c("N","D"))
meta1$race_case=factor(paste(meta1$race1,meta1$case1,sep="_"))
meta1$Race_Case=factor(paste(meta1$Race,meta1$Case,sep="_"))
# stir_norm1=stir_norm1[which(meta$C_MRACE%in%c(1,2)),]
# meta1=meta[which(meta$C_MRACE%in%c(1,2)),]
write.csv(stir_norm1,file="stir_norm1.csv")
filename=read.table(file="/Users/shansun/filenames.txt",sep=" ",quote="")
filename=data.frame(filename[1:2136,])
filename$sample=sapply(strsplit(as.character(filename[,1]),"-"),"[[",1)
filename$R=sapply(strsplit(as.character(filename[,1]),"_"),"[[",2)

filename1=filename[filename$R=="R1.hg19unmapped.fastq",]
filename2=filename[filename$R=="R2.hg19unmapped.fastq",]
filename3=merge(filename1,filename2,by="sample",all=T)
filename4=filename3[match(rownames(stir_norm1),filename3$sample),]
write.csv(filename4,file="filename_format.csv")

stir_shan=vegan::diversity(stir_norm1,index = "shannon", MARGIN = 1, base = exp(1))
stir_simp=vegan::diversity(stir_norm1,index = "simpson", MARGIN = 1, base = exp(1))
stir_invsimp=vegan::diversity(stir_norm1,index = "invsimpson", MARGIN = 1, base = exp(1))
meta1$shannon=stir_shan
meta1$simp=stir_simp
meta1$invsimp=stir_invsimp
meta1$lcrisp=log10(stir_norm1[,177]+1)
meta1$liners=log10(stir_norm1[,183]+1)

#Fig1 race
stir_norm_d1=log10(stir_norm1+1)
gen_pcoa=capscale(stir_norm_d1~1,distance="bray")
var_per=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1:6]*100,2)
pcoa_p=paste("PCoA",c(1:6)," (",var_per,"%)",sep="")

col_d=c("steelblue","tomato")
col2=adjustcolor(col_d[factor(meta1$race1)], alpha.f = 1)
pdf("pcoa_race.pdf",width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,1))
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),display="sites",type="none",cex.lab=2,xlab=pcoa_p[1],ylab=pcoa_p[2])
points(pcoa12,"sites",col=adjustcolor(col2, alpha.f = 0.5),pch=16,cex=1.5)
for (n in 1:2){
  ordiellipse(pcoa12, meta1$race1, kind="se", conf=0.95, lwd=4, draw = "lines", col=col_d[n],show.groups=levels(factor(meta1$race1))[n],label=T,font=2,cex=1) 
}
legend("topright",c("Black","White"),col=col_d,cex=1.5,pch=16,bty = "n")
dev.off()


pdf("shannon_race_box.pdf",width=8,height=6)
ggboxplot(meta1, x = "race1", y = "shannon", color = "race1", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race", y = "Shannon Index")+rotate_x_text(0)
dev.off()

pdf("lcrisp_race_box.pdf",width=6,height=6)
ggboxplot(meta1, x = "race1", y = "lcrisp", color = "race1", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race", y = "L.crispatus")+rotate_x_text(0)
dev.off()

pdf("liners_race_box.pdf",width=6,height=6)
ggboxplot(meta1, x = "race1", y = "liners", color = "race1", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race", y = "L.iners")+rotate_x_text(0)
dev.off()

wilcox.test(meta1$shannon~meta1$race1)
#W = 97923, p-value = 2.137e-05
wilcox.test(meta1$lcrisp~meta1$race1)
#W = 72505, p-value = 0.001153
wilcox.test(meta1$liners~meta1$race1)
#W = 98934, p-value = 5.407e-06

#case
col_d=c("purple","turquoise3")
col2=adjustcolor(col_d[factor(meta1$case1)], alpha.f = 1)
pdf("pcoa_case.pdf",width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,1))
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),display="sites",type="none",cex.lab=2,xlab=pcoa_p[1],ylab=pcoa_p[2])
points(pcoa12,"sites",col=adjustcolor(col2, alpha.f = 0.5),pch=16,cex=1.5)
for (n in 1:2){
  ordiellipse(pcoa12, meta1$case1, kind="se", conf=0.95, lwd=4, draw = "lines", col=col_d[n],show.groups=levels(factor(meta1$case1))[n],label=T,font=2,cex=1) 
}
legend("topright",c("SPTB","Term"),col=col_d,cex=1.5,pch=16,bty = "n")
dev.off()


pdf("shannon_case_box.pdf",width=8,height=6)
ggboxplot(meta1, x = "case1", y = "shannon", color = "case1", add = "jitter", show.legend=TRUE,palette =rev(col_d),
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Case", y = "Shannon Index")+rotate_x_text(0)
dev.off()

pdf("lcrisp_case_box.pdf",width=6,height=6)
ggboxplot(meta1, x = "case1", y = "lcrisp", color = "case1", add = "jitter", show.legend=TRUE,palette =rev(col_d),
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Case", y = "L.crispatus")+rotate_x_text(0)
dev.off()

pdf("liners_case_box.pdf",width=6,height=6)
ggboxplot(meta1, x = "case1", y = "liners", color = "case1", add = "jitter", show.legend=TRUE,palette =rev(col_d),
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Case", y = "L.iners")+rotate_x_text(0)
dev.off()

wilcox.test(meta1$shannon~meta1$case1)
#W = 60691, p-value = 0.1156
wilcox.test(meta1$lcrisp~meta1$case1)
#W = 47946, p-value = 0.002645
wilcox.test(meta1$liners~meta1$case1)
#W = 60547, p-value = 0.128

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
write.csv(ado_R2_vag,file="ado_R2_vag_tax_new.csv")
ado_R2_vag=read.csv(file="ado_R2_vag_tax_new.csv",row.names=1)
ado_R2_vag=ado_R2_vag[ado_R2_vag$FDR<0.1,]
ado_R2_vag=ado_R2_vag[order(ado_R2_vag[,1],decreasing = T),]

bar_lab=rownames(ado_R2_vag)
bar_lab=c("Poverty level","Education","Marital status" ,"Age","Douching","Race","Depression","Negative life events","SPTB","Parity")
#col_bar=c("pink","yellow","lightgreen","lightblue")[factor(rev(c(1,1,2,2,3,3,2,3,3,2,1,3,2,4,3,3,3,3,3,4,3,4,rep(3,3))))]
pdf("adonis_r2_bar.pdf",width=6,height=8)
par(mar=c(5,15,5,5),mfrow=c(1,1))
barplot(rev(as.numeric(ado_R2_vag[,1])),xlim=c(0,0.1),col="palegreen",names.arg=rev(bar_lab[1:nrow(ado_R2_vag)]),las=2,main="Overall",horiz=T,cex.names=1.5,cex.axis=1.5)
dev.off()

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

ggexport(
  plotlist = gplots, filename = paste0("meta_cor_race_case.pdf"),
  ncol = 2, nrow = 3,width=15,height=25
)





#cluster
stir_norm1p=stir_norm1/rowSums(stir_norm1)
vagtypes=vector()
for (i in 1:dim(stir_norm1p)[1]){
  if (max(stir_norm1p[i,])<0.3){
    vagtypes[i]="notype"
  }else{
    vagtypes[i]=names(which(stir_norm1p[i,]==max(stir_norm1p[i,])))
  }
}
sample_cluster=cbind(rownames(stir_norm1p),vagtypes)
colnames(sample_cluster)=c("SampleID","cluster")
#write.csv(sample_cluster,file="JF_cluster_group.csv",row.names = F)

meta1$cluster_JF=vagtypes
meta1$cluster_JF_maj1=meta1$cluster_JF
meta1$cluster_JF_maj1[meta1$cluster_JF_maj1%in%names(which(table(meta1$cluster_JF_maj1)<30))]="Others"
meta1$cluster_JF_maj1=factor(meta1$cluster_JF_maj1,levels=names(sort(table(meta1$cluster_JF_maj1),decreasing = T)))

meta1$cluster_JF_maj=meta1$cluster_JF
meta1$cluster_JF_maj=gsub("Lactobacillus_iners","L.iners",meta1$cluster_JF_maj)
meta1$cluster_JF_maj=gsub("Lactobacillus_crispatus_cluster","L.crispatus",meta1$cluster_JF_maj)
meta1$cluster_JF_maj[grep("Lactobacillus",meta1$cluster_JF_maj)]="Lacto_other"
meta1$cluster_JF_maj[meta1$cluster_JF_maj%in%names(which(table(meta1$cluster_JF_maj)<50))]="Others"
meta1$cluster_JF_maj=factor(meta1$cluster_JF_maj,levels=c("L.crispatus", "L.iners", "Lacto_other", "Others"))
meta1$cluster_JF_maj_case=paste(meta1$cluster_JF_maj,meta1$C_CASE2)
meta1$cluster_JF_maj_race=paste(meta1$cluster_JF_maj,meta1$C_MRACE)
meta1$cluster_JF_maj_race_case=paste(meta1$cluster_JF_maj_race,meta1$C_CASE2)


cl_case=table(paste(meta1$cluster_JF_maj,meta1$case1))
fisher_case=matrix(nrow=4,ncol=8)
for (m in 1:4){
  for (n in 1:4){
    mat1=rbind(cl_case[(m*2-1):(m*2)],cl_case[(n*2-1):(n*2)])
    fisher_case[m,n]=fisher.test(mat1)$estimate
    fisher_case[m,n+4]=fisher.test(mat1)$p.value
  }
}
rownames(fisher_case)=levels(meta1$cluster_JF_maj)
colnames(fisher_case)=paste(rep(levels(meta1$cluster_JF_maj),2),c(rep("r2",4),rep("P",4)))
write.csv(fisher_case,file="fisher_cluster_case.csv")

stir_norm1p=stir_norm1/rowSums(stir_norm1)
stir_norm1pm1=stir_norm1p[,apply(stir_norm1p,2,function(i){length(which(i>0.3))})/829>0.01]
Other=1-rowSums(stir_norm1pm1)
stir_norm1pm=cbind(stir_norm1pm1,Other)
stir_norm1pm=stir_norm1pm[order(meta1$cluster_JF_maj1),]
meta1pm=meta1[order(meta1$cluster_JF_maj1),]
sortnames=c(names(sort(stir_norm1pm[meta1pm$cluster_JF_maj1=="Lactobacillus_crispatus_cluster","Lactobacillus_crispatus_cluster"],decreasing = T)),
            names(sort(stir_norm1pm[meta1pm$cluster_JF_maj1=="Lactobacillus_iners","Lactobacillus_iners"],decreasing = T)),
            names(sort(stir_norm1pm[meta1pm$cluster_JF_maj1=="Lactobacillus_gasseri_cluster","Lactobacillus_gasseri_cluster"],decreasing = T)),
            names(sort(stir_norm1pm[meta1pm$cluster_JF_maj1=="Lactobacillus_jensenii/fornicalis/psittaci","Lactobacillus_jensenii/fornicalis/psittaci"],decreasing = T)),
            names(sort(stir_norm1pm[meta1pm$cluster_JF_maj1=="Lachnospiraceae_BVAB1","Lachnospiraceae_BVAB1"],decreasing = T)),
            names(sort(stir_norm1pm[meta1pm$cluster_JF_maj1=="Gardnerella_vaginalis","Gardnerella_vaginalis"],decreasing = T)))
stir_norm1pm =stir_norm1pm[sortnames,]
stir_norm1pm=stir_norm1pm[,order(colSums(stir_norm1pm))]

colnames(stir_norm1pm)
col_b=c("purple","pink","grey","yellow","orange","lightgreen","cyan","red","blue")
pdf(file="class_bar.pdf",height=12,width=12)
par(mfrow=c(2,1),mar=c(3,5,3,3))
barplot(t(stir_norm1pm),col=col_b,space=0, border = NA,axes=F,axisnames = FALSE,ylab="Percentage",cex.lab=2)
axis(side = 2,cex.axis=1)
plot.new()
legend("left",rev(colnames(stir_norm1pm)),fill=rev(col_b),cex=2,bty="n")
dev.off()

#PCoA plot colored by vagitype 
col12=c("blue","red","grey","lightgreen","yellow","cyan","orange")
col2=adjustcolor(col12[factor(meta1$cluster_JF_maj1)], alpha.f = 1)
pdf("pcoa_cluster_JF_maj.pdf",width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,1))
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),display="sites",type="none",cex.lab=2,xlab=pcoa_p[1],ylab=pcoa_p[2])
points(pcoa12,"sites",col=adjustcolor(col2, alpha.f = 0.5),pch=16,cex=1.5)
# for (i in 1:7){
#   ordiellipse(pcoa12, meta1$cluster_JF_maj1, kind="se", conf=0.95, lwd=4, draw = "lines", col=col12[i],show.groups=levels(factor(meta1$cluster_JF_maj1))[i],label=T,font=2,cex=1)
# }
legend("topright",paste("cluster",names(sort(table(meta1$cluster_JF_maj1),decreasing = T))),col=col12,cex=0.8,pch=16,bty = "n")
dev.off()

#pcoa_race_case
col_r=c("green3","darkgreen","orchid","purple")
col2=adjustcolor(col_r[factor(meta1$Race_Case)], alpha.f = 1)
pdf("pcoa_race_case.pdf",width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,1))
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),display="sites",type="none",cex.lab=2,xlab=pcoa_p[1],ylab=pcoa_p[2])
points(pcoa12,"sites",col=adjustcolor(col2, alpha.f = 0.5),pch=16,cex=1.5)
for (n in 1:4){
  ordiellipse(pcoa12, meta1$Race_Case, kind="se", conf=0.95, lwd=4, draw = "lines", col=col_r[n],show.groups=names(table(meta1$Race_Case))[n],label=T,font=2,cex=1) 
}
legend("topright",paste0(names(table(meta1$race_case))," (",names(table(meta1$Race_Case)),")"),col=col_r,cex=1,pch=16,bty = "n")
dev.off()

adonis(stir_norm_d1~meta1$Race,permutations=999)
"             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
meta1$Race   1     2.067 2.06684  12.746 0.01527  0.001 ***
Residuals  822   133.288 0.16215         0.98473           
Total      823   135.355                 1.00000   "
adonis(stir_norm_d1[meta1$C_MRACE==1,]~meta1$C_CASE2[meta1$C_MRACE==1],permutations=999)
"                                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
meta1$C_CASE2[meta1$C_MRACE == 1]   1     0.410 0.41042  2.4981 0.00538   0.01 **
Residuals                         462    75.905 0.16430         0.99462          
Total                             463    76.316                 1.00000           "
adonis(stir_norm_d1[meta1$C_MRACE==2,]~meta1$C_CASE2[meta1$C_MRACE==2],permutations=999)
"                                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
meta1$C_CASE2[meta1$C_MRACE == 2]   1     0.192 0.19168  1.2085 0.00336   0.24
Residuals                         358    56.780 0.15860         0.99664       
Total                             359    56.972                 1.00000"


shan_tukey=TukeyHSD(aov(stir_shan~meta1$cluster_JF_maj))
pdf("cluster_shannon_box.pdf",width=8,height=5)
ggboxplot(meta1, x = "cluster_JF_maj", y = "shannon", color = "cluster_JF_maj", palette = c("red","blue","green","orange"), add = "jitter", lwd = 0.8, xlab = "Vagitype",ylab = "Shannon Index",
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),title=element_text(size=5),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20)))+rotate_x_text(45)
dev.off()


#stacked cluster composition for race
meta1$cluster_JF_maj_race_case=paste(meta1$cluster_JF_maj,meta1$Race,meta1$case1)
cluster_com=data.frame(table(meta1$cluster_JF_maj_race_case))
cluster_com$perc=cluster_com$Freq/829*100
cluster_com$cluster=sapply(strsplit(as.character(cluster_com$Var1)," "),"[[",1)
cluster_com$race=sapply(strsplit(as.character(cluster_com$Var1)," "),"[[",2)
cluster_com$case=sapply(strsplit(as.character(cluster_com$Var1)," "),"[[",3)
cluster_com$perc_race=cluster_com$Freq/rep(c(360,360,464,464),4)*100
perc_race_case=rep(table(meta1$cluster_JF_maj_race_case)[seq(1,16,2)]/(table(meta1$cluster_JF_maj_race_case)[seq(1,16,2)]+table(meta1$cluster_JF_maj_race_case)[seq(2,16,2)]),each=2)
cluster_com$perc_race_case=perc_race_case*100
perc_case=(cluster_com$Freq[seq(1,16,4)]+cluster_com$Freq[seq(3,16,4)])/by(cluster_com$Freq,cluster_com$cluster,sum)
cluster_com$perc_case=rep(perc_case,each=4)*100

pdf("cluster_perc_bar.pdf",width=8,height=5)
ggplot(data=cluster_com, aes(x=cluster, y=perc_case, fill=cluster)) +geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values=c("red","blue","green","orange"))+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=20),  
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs( y = "SPTB (%)")+rotate_x_text(45)+ylim(0, 35)
dev.off()

pdf("cluster_race_bar_stack.pdf",width=6,height=7)
ggplot(cluster_com, aes(fill=cluster, y=perc_race, x=race)) + 
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values=c("red","blue","green","orange"))+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=20),  
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs( y = "Percentage (%)")+rotate_x_text(45)+ylim(0, 100)
dev.off()

cluster_com1=cluster_com[1:8,]
pdf("cluster_race_perc_bar.pdf",width=5,height=7)
ggplot(data=cluster_com1, aes(x=cluster, y=perc_race_case, fill=race)) +geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values=c("green3","orchid"))+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=20),  
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs( x="Vagitype",y = "SPTB (%)")+rotate_x_text(45)+ylim(0, 30)
dev.off()

#overall percentage of SPTB
race_case=data.frame(table(meta1$race_case))
race_case$Percentage=rep(race_case[c(1,3),2]/(race_case[c(1,3),2]+race_case[c(2,4),2])*100,each=2)
race_case$Race=c("Black","Black","White","White")
race_case$type=rep("overall",4)
pdf("race_perc_bar.pdf",width=4,height=7)
ggplot(data=race_case, aes(x=type, y=Percentage, fill=Race)) +geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values=c("green3","orchid"))+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=20),  
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs( x="",y = "SPTB (%)")+rotate_x_text(45)+ylim(0, 30)
dev.off()


fisher_cluster_race=matrix(nrow=2,ncol=4)
for (i in 1:4){
  a=cluster_com[seq(1,16,2),2]+cluster_com[seq(2,16,2),2]
  b=c(sum(a[seq(1,8,2)]),sum(a[seq(2,8,2)]))
  fisher_cluster_race[1,i]=fisher.test(matrix(c(a[c(2*i-1,2*i)],b-a[c(2*i-1,2*i)]),nrow=2))$estimate
  fisher_cluster_race[2,i]=fisher.test(matrix(c(a[c(2*i-1,2*i)],b-a[c(2*i-1,2*i)]),nrow=2))$p.value
}
fisher_cluster_race

fisher_cluster=matrix(nrow=2,ncol=4)
for (i in 1:4){
  fisher_cluster[1,i]=fisher.test(matrix(cluster_com[c((i*4-3):(i*4)),2],nrow=2))$estimate
  fisher_cluster[2,i]=fisher.test(matrix(cluster_com[c((i*4-3):(i*4)),2],nrow=2))$p.value
}
fisher_cluster


#douche
meta2=meta1[!is.na(meta1$douche),]
meta2$douche_case=factor(paste(meta2$douche,meta2$case1,sep="_"))
meta2$Race_Douche=factor(paste(meta2$Race,meta2$Douche,sep="_"))
meta2$race_douche=factor(paste(meta2$race1,meta2$douche,sep="_"))
meta2$rdc=paste(meta2$race_douche,meta2$Case,sep="_")

stir_norm_d=stir_norm1[match(rownames(meta2),rownames(stir_norm1)),]
stir_norm_d2=log10(stir_norm_d+1)

gen_pcoa=capscale(stir_norm_d2~1,distance="bray")
var_per=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1:6]*100,2)
pcoa_p=paste("PCoA",c(1:6)," (",var_per,"%)",sep="")

col_d=c("darkblue","steelblue","darkred","tomato")
col2=adjustcolor(col_d[factor(meta2$Race_Douche)], alpha.f = 1)
pdf("pcoa_race_douche.pdf",width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,1))
# stir_norm_d2=rbind(stir_norm_d2[-1,],stir_norm_d2[1,])
# meta2=rbind(meta2[-1,],meta2[1,])
# gen_pcoa=capscale(stir_norm_d2~1,distance="bray")
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),display="sites",type="none",cex.lab=2,xlab=pcoa_p[1],ylab=pcoa_p[2])
points(pcoa12,"sites",col=adjustcolor(col2, alpha.f = 0.5),pch=16,cex=1.5)
for (n in 1:4){
  ordiellipse(pcoa12, meta2$Race_Douche, kind="se", conf=0.95, lwd=4, draw = "lines", col=col_d[n],show.groups=levels(factor(meta2$Race_Douche))[n],label=T,font=2,cex=1) 
}
legend("topright",c("Black_douching (B_D)",
                   "Black_no douching (B_N)",
                   "White_douching (W_D)",
                   "White_no douching (W_N)"),col=col_d,cex=1.2,pch=16,bty = "n")
dev.off()

pdf("shannon_race_douche_box.pdf",width=10,height=7)
ggboxplot(meta2, x = "race_douche", y = "shannon", color = "race_douche", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race_Douching", y = "Shannon Index")+rotate_x_text(45)
dev.off()

pdf("lcrisp_race_douche_box.pdf",width=10,height=7)
ggboxplot(meta2, x = "race_douche", y = "lcrisp", color = "race_douche", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race_Douching", y = "L.crispatus")+rotate_x_text(45)
dev.off()

pdf("liners_race_douche_box.pdf",width=10,height=7)
ggboxplot(meta2, x = "race_douche", y = "liners", color = "race_douche", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race_Douching", y = "L.iners")+rotate_x_text(45)
dev.off()

TukeyHSD(aov(meta2$shannon~meta2$race_douche))
"                                           diff        lwr         upr     p adj
Black_No douching-Black_Douching     0.03179357 -0.1784374  0.24202450 0.9798636
White_Douching-Black_Douching       -0.03734413 -0.2325680  0.15787971 0.9606065
White_No douching-Black_Douching    -0.22047352 -0.3892134 -0.05173366 0.0045223
White_Douching-Black_No douching    -0.06913770 -0.2827618  0.14448638 0.8381125
White_No douching-Black_No douching -0.25226709 -0.4419933 -0.06254086 0.0036795
White_No douching-White_Douching    -0.18312939 -0.3560784 -0.01018043 0.0331546"

TukeyHSD(aov(meta2$lcrisp~meta2$race_douche))
"Black_No douching-Black_Douching    -0.15800964 -0.6697084 0.3536891 0.8562224
White_Douching-Black_Douching       -0.19802549 -0.6731972 0.2771462 0.7053621
White_No douching-Black_Douching     0.56172168  0.1510116 0.9724318 0.0025982
White_Douching-Black_No douching    -0.04001585 -0.5599734 0.4799417 0.9972321
White_No douching-Black_No douching  0.71973132  0.2579407 1.1815219 0.0003949
White_No douching-White_Douching     0.75974717  0.3387922 1.1807022 0.0000251"

TukeyHSD(aov(meta2$liners~meta2$race_douche))
"                                           diff        lwr        upr     p adj
Black_No douching-Black_Douching     0.05999801 -0.3543325  0.4743285 0.9822389
White_Douching-Black_Douching        0.04532855 -0.3394254  0.4300826 0.9902655
White_No douching-Black_Douching    -0.74229540 -1.0748539 -0.4097369 0.0000001
White_Douching-Black_No douching    -0.01466946 -0.4356873  0.4063484 0.9997401
White_No douching-Black_No douching -0.80229341 -1.1762125 -0.4283743 0.0000003
White_No douching-White_Douching    -0.78762395 -1.1284779 -0.4467701 0.0000000
"


#taxa associated with race and douche
#wilcoxon
stir_norm_dm=stir_norm_d[,apply(stir_norm_d,2,function(i){length(which(i>38))})>10]
pmat=matrix(nrow=ncol(stir_norm_dm),ncol=4)
pmat_dir=matrix(nrow=ncol(stir_norm_dm),ncol=4)
for (i in 1:ncol(stir_norm_dm)){
  pmat[i,1]=try(wilcox.test(stir_norm_dm[meta2$Douche=="D",i]~meta2$race1[meta2$Douche=="D"])$p.value)
  pmat[i,2]=try(wilcox.test(stir_norm_dm[meta2$Douche=="N",i]~meta2$race1[meta2$Douche=="N"])$p.value)
  pmat[i,3]=try(wilcox.test(stir_norm_dm[meta2$Race=="B",i]~meta2$douche[meta2$Race=="B"])$p.value)
  pmat[i,4]=try(wilcox.test(stir_norm_dm[meta2$Race=="W",i]~meta2$douche[meta2$Race=="W"])$p.value)
  
  pmat_dir[i,1]=try(t.test(stir_norm_dm[meta2$Douche=="D",i]~meta2$race1[meta2$Douche=="D"])$statistic)
  pmat_dir[i,2]=try(t.test(stir_norm_dm[meta2$Douche=="N",i]~meta2$race1[meta2$Douche=="N"])$statistic)
  pmat_dir[i,3]=try(t.test(stir_norm_dm[meta2$Race=="B",i]~meta2$douche[meta2$Race=="B"])$statistic)
  pmat_dir[i,4]=try(t.test(stir_norm_dm[meta2$Race=="W",i]~meta2$douche[meta2$Race=="W"])$statistic)
}
rownames(pmat)=colnames(stir_norm_dm)
fdrs=matrix(p.adjust(pmat,method="fdr"),ncol=4)
fdrs[,1]=p.adjust(pmat[,1],method="fdr")
fdrs[,2]=p.adjust(pmat[,2],method="fdr")
fdrs[,3]=p.adjust(pmat[,3],method="fdr")
fdrs[,4]=p.adjust(pmat[,4],method="fdr")

length(which(fdrs[,1]<0.1)) #0
length(which(fdrs[,2]<0.1)) #23
length(which(fdrs[,3]<0.1)) #0
length(which(fdrs[,4]<0.1)) #22

pmat1=data.frame(-log10(pmat))*sign(pmat_dir)

pdf("race_cor_douche.pdf",width=5,height=5)
p=ggplot(pmat1, mapping=aes(x=X1, y=X2)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P)*direction(change) by race in douching group" , y = "-log10(P)*direction(change) by race in no douching group")+xlim(-5.5, 5.5)+ylim(-5.5, 5.5)
tax_lab=rownames(pmat)
tax_lab[fdrs[,1]>0.1&fdrs[,2]>0.1]=""
p+geom_text_repel(aes(label =tax_lab),size = 3.5)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
dev.off()

pdf("douche_cor_race.pdf",width=5,height=5)
p=ggplot(pmat1, mapping=aes(x=X3, y=X4)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P)*direction(change) by douching in Black women" , y = "-log10(P)*direction(change) by douching in White women")+xlim(-7.5, 7.5)+ylim(-7.5, 7.5)
tax_lab=rownames(pmat)
tax_lab[fdrs[,3]>0.1&fdrs[,4]>0.1]=""
p+geom_text_repel(aes(label =tax_lab),size = 3.5)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
dev.off()

meta2$cluster_JF_maj_race_douche=paste(meta2$cluster_JF_maj,meta2$race1,meta2$douche,sep="-")
cluster_com=data.frame(table(meta2$cluster_JF_maj_race_douche))
cluster_com$perc=cluster_com$Freq/sum(cluster_com$Freq)*100
cluster_com$cluster=sapply(strsplit(as.character(cluster_com$Var1),"-"),"[[",1)
cluster_com$race=sapply(strsplit(as.character(cluster_com$Var1),"-"),"[[",2)
cluster_com$douche=sapply(strsplit(as.character(cluster_com$Var1),"-"),"[[",3)
cluster_com$race_douche=factor(paste(cluster_com$race,cluster_com$douche))
cluster_com$perc_race_douche=cluster_com$Freq/rep(table(paste(meta2$race1,meta2$douche)),4)*100

pdf("cluster_douche_race_bar_stack.pdf",width=8,height=7)
par(mar=c(5,8,5,5))
ggplot(cluster_com, aes(fill=cluster, y=perc_race_douche, x=race_douche)) + 
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values=c("red","blue","green","orange"))+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=15),  
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs(x="Race_Douching", y = "Percentage (%)")+rotate_x_text(45)+ylim(0, 101)
dev.off()

fisher_cluster=matrix(nrow=4,ncol=4)
for (i in 1:4){
  fisher_cluster[1,i]=fisher.test(rbind(cluster_com[c((i*4-3):(i*4-2)),2],table(meta2$race_douche)[1:2]-cluster_com[c((i*4-3):(i*4-2)),2]))$estimate
  fisher_cluster[2,i]=fisher.test(rbind(cluster_com[c((i*4-3):(i*4-2)),2],table(meta2$race_douche)[1:2]-cluster_com[c((i*4-3):(i*4-2)),2]))$p.value
  fisher_cluster[3,i]=fisher.test(rbind(cluster_com[c((i*4-1):(i*4)),2],table(meta2$race_douche)[3:4]-cluster_com[c((i*4-1):(i*4)),2]))$estimate
  fisher_cluster[4,i]=fisher.test(rbind(cluster_com[c((i*4-1):(i*4)),2],table(meta2$race_douche)[3:4]-cluster_com[c((i*4-1):(i*4)),2]))$p.value
}
fisher_cluster



# #taxa associated with proverty in 4 rd groups
# stir_norm_dm=stir_norm_d[,apply(stir_norm_d,2,function(i){length(which(i>38))})>10]
# pmat=matrix(nrow=ncol(stir_norm_dm),ncol=4)
# pmat_dir=matrix(nrow=ncol(stir_norm_dm),ncol=4)
# for (i in 1:ncol(stir_norm_dm)){
#   pmat[i,1]=try(summary(lm(stir_norm_dm[meta2$rd=="1_0",i]~meta2$C_POVERT[meta2$rd=="1_0"]))$coefficients[2,4])
#   pmat[i,2]=try(summary(lm(stir_norm_dm[meta2$rd=="1_1",i]~meta2$C_POVERT[meta2$rd=="1_1"]))$coefficients[2,4])
#   pmat[i,3]=try(summary(lm(stir_norm_dm[meta2$rd=="2_0",i]~meta2$C_POVERT[meta2$rd=="2_0"]))$coefficients[2,4])
#   pmat[i,4]=try(summary(lm(stir_norm_dm[meta2$rd=="2_1",i]~meta2$C_POVERT[meta2$rd=="2_1"]))$coefficients[2,4])
# 
#   pmat_dir[i,1]=try(summary(lm(stir_norm_dm[meta2$rd=="1_0",i]~meta2$C_POVERT[meta2$rd=="1_0"]))$coefficients[2,3])
#   pmat_dir[i,2]=try(summary(lm(stir_norm_dm[meta2$rd=="1_1",i]~meta2$C_POVERT[meta2$rd=="1_1"]))$coefficients[2,3])
#   pmat_dir[i,3]=try(summary(lm(stir_norm_dm[meta2$rd=="2_0",i]~meta2$C_POVERT[meta2$rd=="2_0"]))$coefficients[2,3])
#   pmat_dir[i,4]=try(summary(lm(stir_norm_dm[meta2$rd=="2_1",i]~meta2$C_POVERT[meta2$rd=="2_1"]))$coefficients[2,3])
# }
# rownames(pmat)=colnames(stir_norm_dm)
# fdrs=matrix(p.adjust(pmat,method="fdr"),ncol=4)
# fdrs[,1]=p.adjust(pmat[,1],method="fdr")
# fdrs[,2]=p.adjust(pmat[,2],method="fdr")
# fdrs[,3]=p.adjust(pmat[,3],method="fdr")
# fdrs[,4]=p.adjust(pmat[,4],method="fdr")
# 
# length(which(fdrs[,1]<0.1))
# length(which(fdrs[,2]<0.1))
# length(which(fdrs[,3]<0.1))
# length(which(fdrs[,4]<0.1))
# 
# pmat1=data.frame(-log10(pmat))*sign(pmat_dir)
# 
# gplots=list()
# tax_lab1=rownames(pmat)
# tax_lab1[fdrs[,1]>0.1&fdrs[,2]>0.1]=""
# gplots[[1]]=ggplot(pmat1, mapping=aes(x=X1, y=X2)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P) white women no douching" , y = "-log10(P) white women douching")+xlim(-5.5, 5.5)+ylim(-5.5, 5.5)+
#   geom_text_repel(aes(label =tax_lab1),size = 3.5)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
# 
# tax_lab2=rownames(pmat)
# tax_lab2[fdrs[,3]>0.1&fdrs[,4]>0.1]=""
# gplots[[2]]=ggplot(pmat1, mapping=aes(x=X3, y=X4)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P) black women no douching" , y = "-log10(P) black women douching")+xlim(-5.5, 5.5)+ylim(-5.5, 5.5)+
#   geom_text_repel(aes(label =tax_lab2),size = 3.5)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
# 
# tax_lab3=rownames(pmat)
# tax_lab3[fdrs[,1]>0.1&fdrs[,3]>0.1]=""
# gplots[[3]]=ggplot(pmat1, mapping=aes(x=X1, y=X3)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P) white women no douching" , y = "-log10(P) black women Non-douching")+xlim(-5.5, 5.5)+ylim(-5.5, 5.5)+
#   geom_text_repel(aes(label =tax_lab3),size = 3.5)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
# 
# tax_lab4=rownames(pmat)
# tax_lab4[fdrs[,2]>0.1&fdrs[,4]>0.1]=""
# gplots[[4]]=ggplot(pmat1, mapping=aes(x=X2, y=X4)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P) white women douching" , y = "-log10(P) black women Douching")+xlim(-5.5, 5.5)+ylim(-5.5, 5.5)+
#   geom_text_repel(aes(label =tax_lab4),size = 3.5)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
# 
# ggexport(
#   plotlist = gplots, filename = paste0("cor_poverty.pdf"),
#   ncol = 2, nrow = 2,width=10,height=10
# )


#taxa associated with edu in 4 rd groups
# stir_norm_dm=stir_norm_d[,apply(stir_norm_d,2,function(i){length(which(i>38))})>10]
# pmat=matrix(nrow=ncol(stir_norm_dm),ncol=4)
# pmat_dir=matrix(nrow=ncol(stir_norm_dm),ncol=4)
# for (i in 1:ncol(stir_norm_dm)){
#   pmat[i,1]=try(summary(lm(stir_norm_dm[meta2$rd=="1_0",i]~meta2$C_MOMEDU[meta2$rd=="1_0"]))$coefficients[2,4])
#   pmat[i,2]=try(summary(lm(stir_norm_dm[meta2$rd=="1_1",i]~meta2$C_MOMEDU[meta2$rd=="1_1"]))$coefficients[2,4])
#   pmat[i,3]=try(summary(lm(stir_norm_dm[meta2$rd=="2_0",i]~meta2$C_MOMEDU[meta2$rd=="2_0"]))$coefficients[2,4])
#   pmat[i,4]=try(summary(lm(stir_norm_dm[meta2$rd=="2_1",i]~meta2$C_MOMEDU[meta2$rd=="2_1"]))$coefficients[2,4])
#   
#   pmat_dir[i,1]=try(summary(lm(stir_norm_dm[meta2$rd=="1_0",i]~meta2$C_MOMEDU[meta2$rd=="1_0"]))$coefficients[2,3])
#   pmat_dir[i,2]=try(summary(lm(stir_norm_dm[meta2$rd=="1_1",i]~meta2$C_MOMEDU[meta2$rd=="1_1"]))$coefficients[2,3])
#   pmat_dir[i,3]=try(summary(lm(stir_norm_dm[meta2$rd=="2_0",i]~meta2$C_MOMEDU[meta2$rd=="2_0"]))$coefficients[2,3])
#   pmat_dir[i,4]=try(summary(lm(stir_norm_dm[meta2$rd=="2_1",i]~meta2$C_MOMEDU[meta2$rd=="2_1"]))$coefficients[2,3])
# }
# rownames(pmat)=colnames(stir_norm_dm)
# fdrs=matrix(p.adjust(pmat,method="fdr"),ncol=4)
# fdrs[,1]=p.adjust(pmat[,1],method="fdr")
# fdrs[,2]=p.adjust(pmat[,2],method="fdr")
# fdrs[,3]=p.adjust(pmat[,3],method="fdr")
# fdrs[,4]=p.adjust(pmat[,4],method="fdr")
# 
# length(which(fdrs[,1]<0.1))
# length(which(fdrs[,2]<0.1))
# length(which(fdrs[,3]<0.1))
# length(which(fdrs[,4]<0.1))
# 
# pmat1=data.frame(-log10(pmat))*sign(pmat_dir)
# 
# gplots=list()
# tax_lab1=rownames(pmat)
# tax_lab1[fdrs[,1]>0.1&fdrs[,2]>0.1]=""
# gplots[[1]]=ggplot(pmat1, mapping=aes(x=X1, y=X2)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P) white women no douching" , y = "-log10(P) white women douching")+xlim(-5.5, 5.5)+ylim(-5.5, 5.5)+
#   geom_text_repel(aes(label =tax_lab1),size = 3.5)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
# 
# tax_lab2=rownames(pmat)
# tax_lab2[fdrs[,3]>0.1&fdrs[,4]>0.1]=""
# gplots[[2]]=ggplot(pmat1, mapping=aes(x=X3, y=X4)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P) black women no douching" , y = "-log10(P) black women douching")+xlim(-5.5, 5.5)+ylim(-5.5, 5.5)+
#   geom_text_repel(aes(label =tax_lab2),size = 3.5)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
# 
# tax_lab3=rownames(pmat)
# tax_lab3[fdrs[,1]>0.1&fdrs[,3]>0.1]=""
# gplots[[3]]=ggplot(pmat1, mapping=aes(x=X1, y=X3)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P) white women no douching" , y = "-log10(P) black women no douching")+xlim(-5.5, 5.5)+ylim(-5.5, 5.5)+
#   geom_text_repel(aes(label =tax_lab3),size = 3.5)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
# 
# tax_lab4=rownames(pmat)
# tax_lab4[fdrs[,2]>0.1&fdrs[,4]>0.1]=""
# gplots[[4]]=ggplot(pmat1, mapping=aes(x=X2, y=X4)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P) white women douching" , y = "-log10(P) black women douching")+xlim(-5.5, 5.5)+ylim(-5.5, 5.5)+
#   geom_text_repel(aes(label =tax_lab4),size = 3.5)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
# 
# ggexport(
#   plotlist = gplots, filename = paste0("cor_education.pdf"),
#   ncol = 2, nrow = 2,width=10,height=10
# )


#metadata for race and douching
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

ggexport(
  plotlist = gplots, filename = paste0("meta_cor_race_douche.pdf"),
  ncol = 2, nrow = 3,width=15,height=25
)


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


