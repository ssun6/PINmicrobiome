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
input = file.path(pipeRoot,"input/metadata/")
output = file.path(moduleDir,"output/")

stir_count_melt=read.table(file=paste0(input, "Engel_16S_stirrups_summary_97_070819.txt"),sep="|",quote="")
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


meta_or=read.csv(file=paste0(input, "PIN0108_MicrobCov_NoSNP_20190607.csv"),row.names=1,header=T,quote="")
#new
meta_new=read.csv(file=paste0(input, "PIN0113_MbCov_AllDiet_20200706.csv"),row.names=1,header=T,quote="")
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
write.csv(stir_norm1,file=paste0(output, "stir_norm1.csv"))

# filename=read.table(file="/Users/shansun/filenames.txt",sep=" ",quote="")
# filename=data.frame(filename[1:2136,])
# filename$sample=sapply(strsplit(as.character(filename[,1]),"-"),"[[",1)
# filename$R=sapply(strsplit(as.character(filename[,1]),"_"),"[[",2)
# 
# filename1=filename[filename$R=="R1.hg19unmapped.fastq",]
# filename2=filename[filename$R=="R2.hg19unmapped.fastq",]
# filename3=merge(filename1,filename2,by="sample",all=T)
# filename4=filename3[match(rownames(stir_norm1),filename3$sample),]
# write.csv(filename4,file="filename_format.csv")

stir_shan=vegan::diversity(stir_norm1,index = "shannon", MARGIN = 1, base = exp(1))
stir_simp=vegan::diversity(stir_norm1,index = "simpson", MARGIN = 1, base = exp(1))
stir_invsimp=vegan::diversity(stir_norm1,index = "invsimpson", MARGIN = 1, base = exp(1))
meta1$shannon=stir_shan
meta1$simp=stir_simp
meta1$invsimp=stir_invsimp
meta1$lcrisp=log10(stir_norm1[,177]+1)
meta1$liners=log10(stir_norm1[,183]+1)



#Fig1 race
dir.create(paste0(output, "race/"))
raceOut <- paste0(output, "race/")


stir_norm_d1=log10(stir_norm1+1)
gen_pcoa=capscale(stir_norm_d1~1,distance="bray")
var_per=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1:6]*100,2)
pcoa_p=paste("PCoA",c(1:6)," (",var_per,"%)",sep="")

col_d=c("steelblue","tomato")
col2=adjustcolor(col_d[factor(meta1$race1)], alpha.f = 1)

## ** Figure 1a
pdf(paste0(raceOut, "pcoa_race.pdf"),width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,1))
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),display="sites",type="none",cex.lab=2,xlab=pcoa_p[1],ylab=pcoa_p[2])
points(pcoa12,"sites",col=adjustcolor(col2, alpha.f = 0.5),pch=16,cex=1.5)
for (n in 1:2){
  ordiellipse(pcoa12, meta1$race1, kind="se", conf=0.95, lwd=4, draw = "lines", col=col_d[n],show.groups=levels(factor(meta1$race1))[n],label=T,font=2,cex=1) 
}
legend("topright",c("Black","White"),col=col_d,cex=1.5,pch=16,bty = "n")
dev.off()

wilcox <- wilcox.test(meta1$shannon~meta1$race1)
#W = 97923, p-value = 2.137e-05

## ** Figure 1d
pdf(paste0(raceOut, "shannon_race_box.pdf"),width=8,height=6)
ggboxplot(meta1, x = "race1", y = "shannon", color = "race1", 
          add = "jitter", show.legend=TRUE, palette =col_d, ggtheme = theme_bw()+
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  plot.title= element_text(size=20), 
                  panel.background = element_blank(), 
                  axis.line = element_blank(),
                  axis.text = element_text(size=20), 
                  axis.title=element_text(size=20),
                  legend.text=element_text(size=20),
                  legend.key.size = unit(1,"line"), 
                  legend.title=element_text(size=20)))+
  labs( x = "Race", y = "Shannon Index")+
  ggtitle(paste0("P = ", formatC(wilcox$p.value, format = "e", digits = 2)))+
  theme(plot.title = element_text(hjust = 0.5))+
  rotate_x_text(0)
dev.off()

wilcox <- wilcox.test(meta1$lcrisp~meta1$race1)
#W = 72505, p-value = 0.001153

## ** Figure 1d
pdf(paste0(raceOut, "lcrisp_race_box.pdf"),width=6,height=6)
ggboxplot(meta1, x = "race1", y = "lcrisp", color = "race1", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race", y = "L.crispatus")+  
  ggtitle(paste0("P = ", round(wilcox$p.value, 4)))+
  theme(plot.title = element_text(hjust = 0.5))+
  rotate_x_text(0)
dev.off()

wilcox <- wilcox.test(meta1$liners~meta1$race1)
#W = 98934, p-value = 5.407e-06

## ** Figure 1d
pdf(paste0(raceOut, "liners_race_box.pdf"),width=6,height=6)
ggboxplot(meta1, x = "race1", y = "liners", color = "race1", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race", y = "L.iners")+  
  ggtitle(paste0("P = ", formatC(wilcox$p.value, format = "e", digits = 2)))+
  theme(plot.title = element_text(hjust = 0.5))+
  rotate_x_text(0)
dev.off()



#case
dir.create(paste0(output, "case/"))
caseOut <- paste0(output, "case/")

col_d=c("purple","turquoise3")
col2=adjustcolor(col_d[factor(meta1$case1)], alpha.f = 1)


## ** Figure 1b
pdf(paste0(caseOut, "pcoa_case.pdf"),width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,1))
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),display="sites",type="none",cex.lab=2,xlab=pcoa_p[1],ylab=pcoa_p[2])
points(pcoa12,"sites",col=adjustcolor(col2, alpha.f = 0.5),pch=16,cex=1.5)
for (n in 1:2){
  ordiellipse(pcoa12, meta1$case1, kind="se", conf=0.95, lwd=4, draw = "lines", col=col_d[n],show.groups=levels(factor(meta1$case1))[n],label=T,font=2,cex=1) 
}
legend("topright",c("SPTB","Term"),col=col_d,cex=1.5,pch=16,bty = "n")
dev.off()

wilcox <- wilcox.test(meta1$shannon~meta1$case1)
#W = 60691, p-value = 0.1156

## ** Figure 1e
pdf(paste0(caseOut, "shannon_case_box.pdf"),width=8,height=6)
ggboxplot(meta1, x = "case1", y = "shannon", color = "case1", add = "jitter", show.legend=TRUE,palette =rev(col_d),
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Case", y = "Shannon Index")+
  ggtitle(paste0("P = ", round(wilcox$p.value, 3)))+
  theme(plot.title = element_text(hjust = 0.5))+
  rotate_x_text(0)
dev.off()


wilcox <- wilcox.test(meta1$lcrisp~meta1$case1)
#W = 47946, p-value = 0.002645

## ** Figure 1e
pdf(paste0(caseOut, "lcrisp_case_box.pdf"),width=6,height=6)
ggboxplot(meta1, x = "case1", y = "lcrisp", color = "case1", add = "jitter", show.legend=TRUE,palette =rev(col_d),
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Case", y = "L.crispatus")+
  ggtitle(paste0("P = ", round(wilcox$p.value, 4)))+
  theme(plot.title = element_text(hjust = 0.5))+
  rotate_x_text(0)
dev.off()

wilcox <- wilcox.test(meta1$liners~meta1$case1)
#W = 60547, p-value = 0.128

## ** Figure 1e
pdf(paste0(caseOut, "liners_case_box.pdf"),width=6,height=6)
ggboxplot(meta1, x = "case1", y = "liners", color = "case1", add = "jitter", show.legend=TRUE,palette =rev(col_d),
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),  
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Case", y = "L.iners")+
  ggtitle(paste0("P = ", round(wilcox$p.value, 3)))+
  theme(plot.title = element_text(hjust = 0.5))+
  rotate_x_text(0)
dev.off()


write.csv(meta1,file=paste0(output, "metadata.csv"))
write.csv(stir_norm1,file=paste0(output, "stir_norm1.csv"))

