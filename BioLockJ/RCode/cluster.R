#Author: Shan Sun

#BioLockJ configuration: Ali Sorgen
#Date: 12-08-20
#Description: 

#Libraries
library(reshape2)
library(vegan)
library(ggpubr)

rm(list = ls())

# setwd("~/Documents/BioLockJ_pipelines/vaginalMicrobiome_2020Dec08/03_Cluster/script/")

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/metadata/")
output = file.path(moduleDir,"output/")
racecasePath <- file.path(dir(pipeRoot, full.names = TRUE, pattern = "_RaceCase"), "output/")

meta1=read.csv(file=paste0(racecasePath, "metadata.csv"),row.names=1)


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
write.csv(sample_cluster,file=paste0(output, "JF_cluster_group.csv"),row.names = F)

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
write.csv(fisher_case,file=paste0(output, "fisher_cluster_case.csv"))

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



## ** Figure 2a
pdf(file=paste0(output, "class_bar.pdf"),height=12,width=12)
par(mfrow=c(2,1),mar=c(3,5,3,3))
barplot(t(stir_norm1pm),col=col_b,space=0, border = NA,axes=F,axisnames = FALSE,ylab="Percentage",cex.lab=2)
axis(side = 2,cex.axis=1)
plot.new()
legend("left",rev(colnames(stir_norm1pm)),fill=rev(col_b),cex=2,bty="n")
dev.off()

stir_norm_d1=log10(stir_norm1+1)
gen_pcoa=capscale(stir_norm_d1~1,distance="bray")
var_per=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1:6]*100,2)
pcoa_p=paste("PCoA",c(1:6)," (",var_per,"%)",sep="")


#PCoA plot colored by vagitype 
col12=c("blue","red","grey","lightgreen","yellow","cyan","orange")
col2=adjustcolor(col12[factor(meta1$cluster_JF_maj1)], alpha.f = 1)


## ** Figure 2b
pdf(paste0(output, "pcoa_cluster_JF_maj.pdf"),width=8,height=8)
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


## ** Figure 2c
pdf(paste0(output, "pcoa_race_case.pdf"),width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,1))
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),display="sites",type="none",cex.lab=2,xlab=pcoa_p[1],ylab=pcoa_p[2])
points(pcoa12,"sites",col=adjustcolor(col2, alpha.f = 0.5),pch=16,cex=1.5)
for (n in 1:4){
  ordiellipse(pcoa12, meta1$Race_Case, kind="se", conf=0.95, lwd=4, draw = "lines", col=col_r[n],show.groups=names(table(meta1$Race_Case))[n],label=T,font=2,cex=1) 
}
legend("topright",paste0(names(table(meta1$race_case))," (",names(table(meta1$Race_Case)),")"),col=col_r,cex=1,pch=16,bty = "n")
dev.off()

adonis(stir_norm_d1~meta1$Race,permutations=999)
adonis(stir_norm_d1[meta1$C_MRACE==1,]~meta1$C_CASE2[meta1$C_MRACE==1],permutations=999)
adonis(stir_norm_d1[meta1$C_MRACE==2,]~meta1$C_CASE2[meta1$C_MRACE==2],permutations=999)


stir_shan=vegan::diversity(stir_norm1,index = "shannon", MARGIN = 1, base = exp(1))
shan_tukey=TukeyHSD(aov(stir_shan~meta1$cluster_JF_maj))


## ** Figure 2e
pdf(paste0(output, "cluster_shannon_box.pdf"),width=8,height=5)
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



## ** Figure 2d
pdf(paste0(output, "cluster_perc_bar.pdf"),width=8,height=5)
ggplot(data=cluster_com, aes(x=cluster, y=perc_case, fill=cluster)) +geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values=c("red","blue","green","orange"))+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=20),  
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs( y = "SPTB (%)")+rotate_x_text(45)+ylim(0, 35)
dev.off()


## ** Figure 2f
pdf(paste0(output, "cluster_race_bar_stack.pdf"),width=6,height=7)
ggplot(cluster_com, aes(fill=cluster, y=perc_race, x=race)) + 
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values=c("red","blue","green","orange"))+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=20),  
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs( y = "Percentage (%)")+rotate_x_text(45)+ylim(0, 100)
dev.off()

cluster_com1=cluster_com[1:8,]


## ** Figure 2g
pdf(paste0(output, "cluster_race_perc_bar.pdf"),width=5,height=7)
ggplot(data=cluster_com1, aes(x=cluster, y=perc_race_case, fill=race)) +geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values=c("green3","orchid"))+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=20),  
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs( x="Vagitype",y = "SPTB (%)")+rotate_x_text(45)+ylim(0, 30)
dev.off()

#overall percentage of SPTB
race_case=data.frame(table(meta1$race_case))
race_case$Percentage=rep(race_case[c(1,3),2]/(race_case[c(1,3),2]+race_case[c(2,4),2])*100,each=2)
race_case$Race=c("Black","Black","White","White")
race_case$type=rep("overall",4)


## ** Figure 2g
pdf(paste0(output, "race_perc_bar.pdf"),width=4,height=7)
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


write.csv(meta1,file=paste0(output, "metadata1.csv"))
