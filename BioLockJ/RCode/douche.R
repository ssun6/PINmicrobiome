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

# setwd("~/Documents/BioLockJ_pipelines/vaginalMicrobiome_2020Dec08/04_Douche/script/")

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/metadata/")
output = file.path(moduleDir,"output/")
clusterPath <- file.path(dir(pipeRoot, full.names = TRUE, pattern = "Cluster"), "output/")

stir_count_melt=read.table(file=paste0(input, "Engel_16S_stirrups_summary_97_070819.txt"),sep="|",quote="")
stir_count_melt1=stir_count_melt[,c(1,2,4)]
stir_count_melt1[,3]=as.numeric(as.character(stir_count_melt1[,3]))
stir_count2=acast(stir_count_melt1, V1~ V2,sum)
stir_count2[is.na(stir_count2)]=0

stir_count_melt_p1=stir_count_melt[stir_count_melt[,3]=="AT",c(1,2,5)]
stir_count_p=acast(stir_count_melt_p1, V1~ V2)
stir_count_p[is.na(stir_count_p)]=0
stir_count_p=cbind(stir_count_p,100-rowSums(stir_count_p))
colnames(stir_count_p)[337]="other"
stir_count_p=stir_count_p[which(rowSums(stir_count2)>1000),]
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

stir_norm=stir_count/rowSums(stir_count)*mean(rowSums(stir_count))
meta=meta[intersect(rownames(stir_norm),rownames(meta)),]
stir_norm1=stir_norm[intersect(rownames(stir_norm),rownames(meta)),]
stir_norm1=stir_norm1[which(meta$C_MRACE%in%c(1,2)&meta$WHYPTB2%in%c(0,1,2)),]

meta1=read.csv(file=paste0(clusterPath, "metadata1.csv"),row.names=1)

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


## ** Figure 3a
pdf(paste0(output, "pcoa_race_douche.pdf"),width=8,height=8)
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



## ** Figure 3c
pdf(paste0(output, "shannon_race_douche_box.pdf"),width=10,height=7)
ggboxplot(meta2, x = "race_douche", y = "shannon", color = "race_douche", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race_Douching", y = "Shannon Index")+rotate_x_text(45)
dev.off()



## ** Figure 3c
pdf(paste0(output, "lcrisp_race_douche_box.pdf"),width=10,height=7)
ggboxplot(meta2, x = "race_douche", y = "lcrisp", color = "race_douche", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race_Douching", y = "L.crispatus")+rotate_x_text(45)
dev.off()



## ** Figure 3c
pdf(paste0(output, "liners_race_douche_box.pdf"),width=10,height=7)
ggboxplot(meta2, x = "race_douche", y = "liners", color = "race_douche", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race_Douching", y = "L.iners")+rotate_x_text(45)
dev.off()

TukeyHSD(aov(meta2$shannon~meta2$race_douche))

TukeyHSD(aov(meta2$lcrisp~meta2$race_douche))

TukeyHSD(aov(meta2$liners~meta2$race_douche))


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



## ** Figure 3d
pdf(paste0(output, "race_cor_douche.pdf"),width=5,height=5)
p=ggplot(pmat1, mapping=aes(x=X1, y=X2)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P) by race in douching group" , y = "-log10(P) by race in no douching group")+xlim(-5.5, 5.5)+ylim(-5.5, 5.5)
tax_lab=rownames(pmat)
tax_lab[fdrs[,1]>0.1&fdrs[,2]>0.1]=""
p+geom_text_repel(aes(label =tax_lab),size = 3.5)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
dev.off()



## ** Figure 3e
pdf(paste0(output, "douche_cor_race.pdf"),width=5,height=5)
p=ggplot(pmat1, mapping=aes(x=X3, y=X4)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P) by douching in Black women" , y = "-log10(P) by douching in White women")+xlim(-7.5, 7.5)+ylim(-7.5, 7.5)
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



## ** Figure 3b
pdf(paste0(output, "cluster_douche_race_bar_stack.pdf"),width=8,height=7)
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

write.csv(meta2,file=paste0(output, "metadata2.csv"))
