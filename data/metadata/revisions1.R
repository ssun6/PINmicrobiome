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
library(ALDEx2)

meta1=read.csv(file="/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/meta824.csv",row.names=1,header = T)
write.csv(meta1,file="/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/meta824_public.csv")

meta1=read.csv(file="/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/meta824_public.csv",row.names=1,header = T)
meta1=meta1[,c(9,18,342,343,339,355,356,357,358,359,361,362,363)]

stir_norm=read.csv(file="/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/stir824.csv",header=T,row.names=1)
stir_norm1=stir_norm/rowSums(stir_norm)*mean(rowSums(stir_norm))

#Fig1 d,e,f
stir_shan=vegan::diversity(stir_norm,index = "shannon", MARGIN = 1, base = exp(1))
meta1$shannon=stir_shan
meta1$lcrisp=log10(stir_norm1[,266]+1)
meta1$liners=log10(stir_norm1[,272]+1)

col_d=c("steelblue","lightblue","tomato","pink")
pdf("/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/shannon_race_case_box.pdf",width=10,height=6)
ggboxplot(meta1, x = "race_case", y = "shannon", color = "race_case", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race", y = "Shannon Index")+rotate_x_text(45)
dev.off()

pdf("/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/lcrisp_race_case_box.pdf",width=10,height=6)
ggboxplot(meta1, x = "race_case", y = "lcrisp", color = "race_case", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race", y = "L.crispatus")+rotate_x_text(45)
dev.off()

pdf("/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/liners_race_case_box.pdf",width=10,height=6)
ggboxplot(meta1, x = "race_case", y = "liners", color = "race_case", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race", y = "L.iners")+rotate_x_text(45)
dev.off()

#Wilcoxon test of shannon diversity
wilcox.test(meta1$shannon[meta1$race_case=="White_Term"],meta1$shannon[meta1$race_case=="Black_Term"])
#W = 43915, p-value = 0.000954
wilcox.test(meta1$shannon[meta1$race_case=="White_SPTB"],meta1$shannon[meta1$race_case=="Black_SPTB"])
#W = 2841, p-value = 0.006471
wilcox.test(meta1$shannon[meta1$race_case=="White_Term"],meta1$shannon[meta1$race_case=="White_SPTB"])
#W = 16271, p-value = 0.7145
wilcox.test(meta1$shannon[meta1$race_case=="Black_Term"],meta1$shannon[meta1$race_case=="Black_SPTB"])
#W = 10282, p-value = 0.1169

#aldex2 test of l.iners and l.crisp
aldex_rc=list()
levels1=levels(factor(meta1$rc))
k=1
for (i in 1:3){
  for (j in c(2:4)[i:3]){
    level2=levels1[c(i,j)]
    meta3=meta1[meta1$rc%in%level2,]
    stir_norm3=stir_norm[rownames(meta3),]
    aldex_rc1 <- aldex(t(stir_norm3), meta3$rc, mc.samples=128,test="t",effect=TRUE,include.sample.summary=FALSE,verbose=TRUE)
    aldex_rc[[k]]=aldex_rc1
    k=k+1
  }
}
levels1
#"1_0" "1_1" "2_0" "2_1"

#"1_0 vs 1_1"
aldex_rc[[1]][c("Lactobacillus_iners","Lactobacillus_crispatus_cluster"),]
#we.ep    we.eBH     wi.ep    wi.eBH
# Lactobacillus_iners             0.128115000 0.8548637 0.14221458 0.8523197
# Lactobacillus_crispatus_cluster 0.003188132 0.3560813 0.01261505 0.6142383

#"1_0 vs 2_0"
aldex_rc[[2]][c("Lactobacillus_iners","Lactobacillus_crispatus_cluster"),]
#we.ep    we.eBH     wi.ep    wi.eBH
# Lactobacillus_iners             8.429624e-08 0.0000336878 1.176172e-05 0.001595333
# Lactobacillus_crispatus_cluster 1.204583e-03 0.0433724074 1.031871e-03 0.040797117

#"1_1 vs 2_1"
aldex_rc[[5]][c("Lactobacillus_iners","Lactobacillus_crispatus_cluster"),]
#we.ep    we.eBH     wi.ep    wi.eBH
# Lactobacillus_iners             0.097482 0.6876615 0.5038633 0.9007496
# Lactobacillus_crispatus_cluster 0.334709 0.8551282 0.2379326 0.7998107

#"2_0 vs 2_1"
aldex_rc[[6]][c("Lactobacillus_iners","Lactobacillus_crispatus_cluster"),]
#we.ep    we.eBH     wi.ep    wi.eBH
# Lactobacillus_iners             0.90838351 0.9885254 0.5939503 0.9620750
# Lactobacillus_crispatus_cluster 0.06821294 0.8500248 0.1110376 0.8786263

#Table S3
#merge results to a table
aldex_rc1=merge(aldex_rc[[1]],aldex_rc[[2]],by=0, all=T)
rownames(aldex_rc1)=aldex_rc1[,1]
for (k in 3:6){
  aldex_rc1=merge(aldex_rc1,aldex_rc[[k]],by=0, all=T)
  rownames(aldex_rc1)=aldex_rc1[,1]
}

aldex_rc2=-log10(aldex_rc1[,c(13,24,35,46,57,68)])*sign(aldex_rc1[,c(11,22,33,44,55,66)])
aldex_rc3=aldex_rc1[,c(14,25,36,47,58,69)]
aldex_rc4=aldex_rc3[order(apply(aldex_rc3,1,min)),]
n=length(which(apply(aldex_rc3,1,min)<0.1))
aldex_rc5=aldex_rc4[1:n,]
colnames(aldex_rc5)=c("white_term vs white_SPTB","white_term vs black_term","white_term vs black_SPTB","white_SPTB vs black_term","white_SPTB vs black_SPTB","black_term vs black_SPTB")
aldex_rc6=aldex_rc5[,c(1,5,2,6)]#FDRs

aldex_rc7=aldex_rc1[,-c(1:5)]
colnames(aldex_rc7)=gsub("1_0","W_T",colnames(aldex_rc7))
colnames(aldex_rc7)=gsub("1_1","W_S",colnames(aldex_rc7))
colnames(aldex_rc7)=gsub("2_0","B_T",colnames(aldex_rc7))
colnames(aldex_rc7)=gsub("2_1","B_S",colnames(aldex_rc7))
write.csv(aldex_rc7,"/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/aldex2_rc.csv")
aldex_rc8=aldex_rc7[,c(6,8,9,17,19,20,50,52,53,61,63,64)]
colnames(aldex_rc8)=paste(rep(c("effect","we.ep","we.eBH"),4),rep(c("white_term vs white_SPTB","white_term vs black_term","white_SPTB vs black_SPTB","black_term vs black_SPTB"),each=3))
#table S3
write.csv(aldex_rc8,"/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/aldex2_rc_format.csv")

#Fig.3 c,d,e and Table S4
#remove missing douching data
meta2=meta1[which(!is.na(meta1$douche)),]
stir_norm2=stir_norm[rownames(meta2),]
match(rownames(meta2),rownames(stir_norm2))

#Wilcoxon test of shannon diversity
wilcox.test(meta2$shannon[meta2$race_douche=="White_No douching"],meta2$shannon[meta2$race_douche=="White_Douching"])
#W = 8083, p-value = 0.003853
wilcox.test(meta2$shannon[meta2$race_douche=="White_No douching"],meta2$shannon[meta2$race_douche=="Black_Douching"])
#W = 8399, p-value = 0.0007118
wilcox.test(meta2$shannon[meta2$race_douche=="White_No douching"],meta2$shannon[meta2$race_douche=="Black_No douching"])
#W = 5702, p-value = 0.0005974

aldex_rd=list()
levels1=levels(factor(meta2$rd))
k=1
for (i in 1:3){
  for (j in c(2:4)[i:3]){
    level2=levels1[c(i,j)]
    meta3=meta2[meta2$rd%in%level2,]
    stir_norm3=stir_norm2[rownames(meta3),]
    aldex_rd1 <- aldex(t(stir_norm3), meta3$rd, mc.samples=128,test="t",effect=TRUE,include.sample.summary=FALSE,verbose=TRUE)
    aldex_rd[[k]]=aldex_rd1
    k=k+1
  }
}
levels1
#"1_0" "1_1" "2_0" "2_1"

#"1_0 vs 1_1"
aldex_rd[[1]][c("Lactobacillus_iners","Lactobacillus_crispatus_cluster"),]
# we.ep       we.eBH        wi.ep       wi.eBH
# Lactobacillus_iners             1.263442e-08 4.914790e-06 6.646421e-08 2.585458e-05
# Lactobacillus_crispatus_cluster 4.153328e-06 7.889718e-04 2.070182e-06 3.997272e-04

#"1_0 vs 2_0"
aldex_rd[[2]][c("Lactobacillus_iners","Lactobacillus_crispatus_cluster"),]

# we.ep       we.eBH        wi.ep      wi.eBH
# Lactobacillus_iners             6.764561e-07 0.0002515908 7.890312e-06 0.002616448
# Lactobacillus_crispatus_cluster 1.481512e-05 0.0026294012 2.715074e-05 0.004610929

#"1_0 vs 2_1"
aldex_rd[[3]][c("Lactobacillus_iners","Lactobacillus_crispatus_cluster"),]
# we.ep       we.eBH        wi.ep       wi.eBH
# Lactobacillus_iners             8.959264e-08 0.0000352995 1.356080e-06 0.0005339863
# Lactobacillus_crispatus_cluster 3.412727e-04 0.0311048222 1.330178e-04 0.0141148814

aldex_rd1=merge(aldex_rd[[1]],aldex_rd[[2]],by=0, all=T)
rownames(aldex_rd1)=aldex_rd1[,1]
for (k in 3:6){
  aldex_rd1=merge(aldex_rd1,aldex_rd[[k]],by=0, all=T)
  rownames(aldex_rd1)=aldex_rd1[,1]
}

aldex_rd2=-log10(aldex_rd1[,c(13,24,35,46,57,68)])*sign(aldex_rd1[,c(11,22,33,44,55,66)])
aldex_rd3=aldex_rd1[,c(14,25,58,69)]
length(which(aldex_rd3[,1]<0.1))#3
length(which(aldex_rd3[,2]<0.1))#4
length(which(aldex_rd3[,3]<0.1))#0
length(which(aldex_rd3[,4]<0.1))#0
write.csv(aldex_rd1[,-c(1:5)],"/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/aldex2_rd.csv")
aldex_rd4=aldex_rd3[order(apply(aldex_rd3,1,min)),]
n=length(which(apply(aldex_rd3,1,min)<0.1))
aldex_rd5=aldex_rd4[1:n,]
colnames(aldex_rd5)=c("white_no douching vs white douching","white_no douching vs black_no douching","white_douching vs black douching","black_no douching vs black_douching")

aldex_rd7=aldex_rd1[,-c(1:5)]
colnames(aldex_rd7)=gsub("1_0","white_no douching",colnames(aldex_rd7))
colnames(aldex_rd7)=gsub("1_1","white douching",colnames(aldex_rd7))
colnames(aldex_rd7)=gsub("2_0","black_no douching",colnames(aldex_rd7))
colnames(aldex_rd7)=gsub("2_1","black douching",colnames(aldex_rd7))
write.csv(aldex_rd7,"/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/aldex2_rd.csv")
aldex_rd8=aldex_rd7[,c(6,8,9,17,19,20,50,52,53,61,63,64)]
colnames(aldex_rd8)=paste(rep(c("effect","we.ep","we.eBH"),4),rep(c("white_no douching vs white_douching","white_no douching vs black_no douching","white_douching vs black douching","black_no douching vs black douching"),each=3))
#table S4
write.csv(aldex_rd8,"/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/aldex2_rd_format.csv")

#Fig.3d
rownames(aldex_rd2)[grep("Gard",rownames(aldex_rd2))]="Gardnerella_spp"
pdf("/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/race_douche.pdf",width=6,height=6)
p=ggplot(aldex_rd2, mapping=aes(x=we.ep.x.2, y=we.ep.y)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P)*direction(change) by race in douching group" , y = "-log10(P)*direction(change) by race in no douching group")+xlim(-8, 8)+ylim(-8, 8)
tax_lab=rownames(aldex_rd2)
tax_lab[aldex_rd3[,3]>0.5&aldex_rd3[,2]>0.5]=""
p+geom_text_repel(aes(label =tax_lab),size = 3.5,max.overlaps=20)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
dev.off()

#Fig.3e
pdf("/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/douche_race.pdf",width=6,height=6)
p=ggplot(-aldex_rd2, mapping=aes(x=we.ep.y.2, y=we.ep.x)) +geom_point(color = 'blue')+ theme_classic(base_size = 15) + labs(title="",x ="-log10(P)*direction(change) by douching in Black women" , y = "-log10(P)*direction(change) by douching in white women")+xlim(-8, 8)+ylim(-8, 8)
tax_lab=rownames(aldex_rd2)
tax_lab[aldex_rd3[,4]>0.5&aldex_rd3[,1]>0.5]=""
p+geom_text_repel(aes(label =tax_lab),size = 3.5,max.overlaps=100)+geom_hline(yintercept=0, linetype="dotted")+geom_vline(xintercept=0, linetype="dotted")+theme(axis.text=element_text(size=12),axis.title=element_text(size=12))
dev.off()


#Table 2a
#logistic regression models
meta1=read.csv(file="/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/meta824.csv",row.names=1,header = T)
stir_norm=read.csv(file="/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/stir824.csv",header=T,row.names=1)
stir_norm1=stir_norm/rowSums(stir_norm)*mean(rowSums(stir_norm))

stir_shan=vegan::diversity(stir_norm,index = "shannon", MARGIN = 1, base = exp(1))
meta1$shannon=stir_shan
meta1$lcrisp=log10(stir_norm1[,266]+1)
meta1$liners=log10(stir_norm1[,272]+1)

log_mat=matrix(nrow=3, ncol=4)
rownames(log_mat)=c('lcrisp',"liners","Shannon")
colnames(log_mat)=c('Overall',"Black","white",'interaction_P')

mylogit1 <- glm(C_CASE2 ~ lcrisp+Race+edu_cat+C_BMI+CIGS_46A, data = meta1, family = "binomial")
summary(mylogit1)
#lcrisp       -0.212839   0.075121  -2.833  0.00461 **
exp(summary(mylogit1)$coef[2,1])#0.616181
exp(confint.default(mylogit1, level = 0.95))
#lcrisp       0.6976243 0.9365016
log_mat[1,1]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")

mylogit1 <- glm(C_CASE2 ~ liners+Race+edu_cat+C_BMI+CIGS_46A, data = meta1, family = "binomial")
summary(mylogit1)
#liners        0.05814    0.09058   0.642   0.5210 
exp(summary(mylogit1)$coef[2,1])#1.059858
exp(confint.default(mylogit1, level = 0.95))
#liners       0.8874457 1.265767
log_mat[2,1]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")


mylogit1 <- glm(C_CASE2 ~ shannon+Race+edu_cat+C_BMI+CIGS_46A, data = meta1, family = "binomial")
summary(mylogit1)
#shannon       0.17987    0.15880   1.133   0.2573
exp(summary(mylogit1)$coef[2,1])#1.197058
exp(confint.default(mylogit1, level = 0.95))
#shannon      0.8768959 1.634113
log_mat[3,1]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")


meta1n=meta1[na.omit(which(meta1$race1=="Black")),]
mylogit1 <- glm(C_CASE2 ~ lcrisp+edu_cat+C_BMI+CIGS_46A, data = meta1n, family = "binomial")
summary(mylogit1)
#lcrisp       -1.660e-01  1.146e-01  -1.448    0.148
exp(summary(mylogit1)$coef[2,1])#0.8470208
exp(confint.default(mylogit1, level = 0.95))
#lcrisp       0.67656825 1.060417
log_mat[1,2]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")

mylogit1 <- glm(C_CASE2 ~ liners+edu_cat+C_BMI+CIGS_46A, data = meta1n, family = "binomial")
summary(mylogit1)
#liners        1.843e-02  1.477e-01   0.125    0.901 
exp(summary(mylogit1)$coef[2,1])#1.018597
exp(confint.default(mylogit1, level = 0.95))
#liners       0.7626445 1.360450
log_mat[2,2]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")


mylogit1 <- glm(C_CASE2 ~ shannon+edu_cat+C_BMI+CIGS_46A, data = meta1n, family = "binomial")
summary(mylogit1)
#shannon       2.638e-01  2.247e-01   1.174   0.2404 
exp(summary(mylogit1)$coef[2,1])#1.301873
exp(confint.default(mylogit1, level = 0.95))
#shannon      0.83806320 2.0223704
log_mat[3,2]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")

meta1n=meta1[na.omit(which(meta1$race1=="White")),]
mylogit1 <- glm(C_CASE2 ~ lcrisp+edu_cat+C_BMI+CIGS_46A, data = meta1n, family = "binomial")
summary(mylogit1)
#lcrisp       -0.22845    0.10048  -2.274   0.0230 *
exp(summary(mylogit1)$coef[2,1])#0.7957655
exp(confint.default(mylogit1, level = 0.95))
#lcrisp       0.6535191 0.9689736
log_mat[1,3]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")

mylogit1 <- glm(C_CASE2 ~ liners+edu_cat+C_BMI+CIGS_46A, data = meta1n, family = "binomial")
summary(mylogit1)
#liners        0.05825    0.11829   0.492   0.6224 
exp(summary(mylogit1)$coef[2,1])#1.059978
exp(confint.default(mylogit1, level = 0.95))
#liners       0.8406392 1.336547
log_mat[2,3]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")

mylogit1 <- glm(C_CASE2 ~ shannon+edu_cat+C_BMI+CIGS_46A, data = meta1n, family = "binomial")
summary(mylogit1)
#shannon       0.11032    0.22825   0.483   0.6289 
exp(summary(mylogit1)$coef[2,1])#1.116638
exp(confint.default(mylogit1, level = 0.95))
#shannon      0.7138775 1.7466296
log_mat[3,3]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")

mylogit3 <- glm(C_CASE2 ~ lcrisp*Race+edu_cat+C_BMI+CIGS_46A, data = meta1, family = "binomial")
summary(mylogit3)
#lcrisp:RaceW -0.066200   0.151388  -0.437    0.662
log_mat[1,4]=round(summary(mylogit3)$coefficients[9,4],2)

mylogit3 <- glm(C_CASE2 ~ liners*Race+edu_cat+C_BMI+CIGS_46A, data = meta1, family = "binomial")
summary(mylogit3)
#liners:RaceW  0.03655    0.18576   0.197   0.8440
log_mat[2,4]=round(summary(mylogit3)$coefficients[9,4],2)

mylogit3 <- glm(C_CASE2 ~ shannon*Race+edu_cat+C_BMI+CIGS_46A, data = meta1, family = "binomial")
summary(mylogit3)
#shannon:RaceW -0.140812   0.318199  -0.443   0.6581 
log_mat[3,4]=round(summary(mylogit3)$coefficients[9,4],2)
write.csv(log_mat,file="/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/table2a.csv")

#table 2b
log_mat=matrix(nrow=4, ncol=4)
rownames(log_mat)=c('lcrisp',"liners","Lacto_other","Others")
colnames(log_mat)=c('Overall',"Black","white",'interaction_P')

meta1$cluster_JF_maj_new=factor(meta1$cluster_JF_maj,levels=c( "L.iners","L.crispatus","Lacto_other","Others"))
mylogit1 <- glm(C_CASE2 ~ cluster_JF_maj_new+Race+edu_cat+C_BMI+CIGS_46A, data = meta1, family = "binomial")
summary(mylogit1)
# cluster_JF_maj_newL.crispatus -0.533280   0.276249  -1.930   0.0536 .
# cluster_JF_maj_newLacto_other  0.418357   0.294759   1.419   0.1558  
# cluster_JF_maj_newOthers       0.383712   0.243050   1.579   0.1144  
exp(summary(mylogit1)$coef[2:4,1])
# cluster_JF_maj_newL.crispatus cluster_JF_maj_newLacto_other      cluster_JF_maj_newOthers 
# 0.5866777                     1.5194632                     1.4677220 
exp(confint.default(mylogit1, level = 0.95))
#cluster_JF_maj_newL.crispatus 0.3413942 1.008192
# cluster_JF_maj_newLacto_other 0.8526895 2.707631
# cluster_JF_maj_newOthers      0.9115054 2.363352
log_mat[1,1]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")
log_mat[3,1]=paste0(round(exp(summary(mylogit1)$coef[3,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[3,],2),collapse = " "),")")
log_mat[4,1]=paste0(round(exp(summary(mylogit1)$coef[4,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[4,],2),collapse = " "),")")


meta1n=meta1[na.omit(which(meta1$race1=="Black")),]
mylogit1 <- glm(C_CASE2 ~ cluster_JF_maj_new+edu_cat+C_BMI+CIGS_46A, data = meta1n, family = "binomial")
summary(mylogit1)
exp(summary(mylogit1)$coef[2:4,1])
exp(confint.default(mylogit1, level = 0.95))
log_mat[1,2]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")
log_mat[3,2]=paste0(round(exp(summary(mylogit1)$coef[3,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[3,],2),collapse = " "),")")
log_mat[4,2]=paste0(round(exp(summary(mylogit1)$coef[4,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[4,],2),collapse = " "),")")


meta1n=meta1[na.omit(which(meta1$race1=="White")),]
mylogit1 <- glm(C_CASE2 ~ cluster_JF_maj_new+edu_cat+C_BMI+CIGS_46A, data = meta1n, family = "binomial")
summary(mylogit1)
exp(summary(mylogit1)$coef[2:4,1])
exp(confint.default(mylogit1, level = 0.95))
log_mat[1,3]=paste0(round(exp(summary(mylogit1)$coef[2,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[2,],2),collapse = " "),")")
log_mat[3,3]=paste0(round(exp(summary(mylogit1)$coef[3,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[3,],2),collapse = " "),")")
log_mat[4,3]=paste0(round(exp(summary(mylogit1)$coef[4,1]),2)," (",paste(round(exp(confint.default(mylogit1, level = 0.95))[4,],2),collapse = " "),")")

mylogit3 <- glm(C_CASE2 ~ cluster_JF_maj_new*Race+edu_cat+C_BMI+CIGS_46A, data = meta1, family = "binomial")
summary(mylogit3)
#cluster_JF_maj_newL.crispatus:RaceW  0.008617   0.573063   0.015   0.9880  
#cluster_JF_maj_newLacto_other:RaceW -0.884803   0.625357  -1.415   0.1571  
#cluster_JF_maj_newOthers:RaceW      -0.493630   0.489631  -1.008   0.3134  
log_mat[c(1,3,4),4]=round(summary(mylogit3)$coefficients[11:13,4],2)
log_mat[2,]="ref"
write.csv(log_mat,file="/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/table2b.csv")

#table 1
meta1=read.csv(file="/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/meta824.csv",row.names=1,header = T)

aggregate(meta1$C_MAGE24~meta1$race_case,FUN=mean)
"  meta1$race_case meta1$C_MAGE24
1      Black_SPTB       25.35714
2      Black_Term       24.08696
3      White_SPTB       26.89888
4      White_Term       27.43733"
aggregate(meta1$C_MAGE24~meta1$race_case,FUN=sd)
"  meta1$race_case meta1$C_MAGE24
1      Black_SPTB       5.830214
2      Black_Term       5.358396
3      White_SPTB       7.168500
4      White_Term       6.392734"
by(meta1$C_MAGE24,meta1$race_case,function(i){length(which(is.na(i)))})
"[1] 0
------------------------------------------------------------------------------------- 
meta1$race_case: Black_Term
[1] 0
------------------------------------------------------------------------------------- 
meta1$race_case: White_SPTB
[1] 0
------------------------------------------------------------------------------------- 
meta1$race_case: White_Term
[1] 0"
wilcox.test(meta1$C_MAGE24[meta1$Race=="W"]~meta1$case[meta1$Race=="W"])
#W = 15564, p-value = 0.3226
wilcox.test(meta1$C_MAGE24[meta1$Race=="B"]~meta1$case[meta1$Race=="B"])
#W = 13051, p-value = 0.08007


aggregate(meta1$C_MOMEDU~meta1$race_case,FUN=mean)
"  meta1$race_case meta1$C_MOMEDU
1      Black_SPTB       12.71429
2      Black_Term       12.38909
3      White_SPTB       13.12360
4      White_Term       14.07467"
aggregate(meta1$C_MOMEDU~meta1$race_case,FUN=sd)
"  meta1$race_case meta1$C_MOMEDU
1      Black_SPTB       1.668788
2      Black_Term       1.995769
3      White_SPTB       2.823683
4      White_Term       3.316588"
by(meta1$C_MOMEDU,meta1$race_case,function(i){length(which(is.na(i)))})
"meta1$race_case: Black_SPTB
[1] 0
------------------------------------------------------------------------------------- 
meta1$race_case: Black_Term
[1] 1
------------------------------------------------------------------------------------- 
meta1$race_case: White_SPTB
[1] 0
------------------------------------------------------------------------------------- 
meta1$race_case: White_Term
[1] 0"
wilcox.test(meta1$C_MOMEDU[meta1$Race=="B"]~meta1$case[meta1$Race=="B"])
#W = 12823, p-value = 0.1145
wilcox.test(meta1$C_MOMEDU[meta1$Race=="W"]~meta1$case[meta1$Race=="W"])
#W = 13996, p-value = 0.01655



aggregate(meta1$C_BMI~meta1$race_case,FUN=mean)
"  meta1$race_case meta1$C_BMI
1      Black_SPTB    27.72654
2      Black_Term    27.43577
3      White_SPTB    25.22310
4      White_Term    25.56359"
aggregate(meta1$C_BMI~meta1$race_case,FUN=sd)
"  meta1$race_case meta1$C_BMI
1      Black_SPTB    7.888277
2      Black_Term    7.989241
3      White_SPTB    6.556866
4      White_Term    6.926640"
by(meta1$C_BMI,meta1$race_case,function(i){length(which(is.na(i)))})
"meta1$race_case: Black_SPTB
[1] 6
------------------------------------------------------------------------------------- 
meta1$race_case: Black_Term
[1] 9
------------------------------------------------------------------------------------- 
meta1$race_case: White_SPTB
[1] 2
------------------------------------------------------------------------------------- 
meta1$race_case: White_Term
[1] 7"

wilcox.test(meta1$C_BMI[meta1$Race=="B"]~meta1$case[meta1$Race=="B"])
#W = 10751, p-value = 0.6632
wilcox.test(meta1$C_BMI[meta1$Race=="W"]~meta1$case[meta1$Race=="W"])
#W = 15716, p-value = 0.7919


a=data.frame(table(paste(meta1$race_case,meta1$CIGS_46A)))
b=rep(table(meta1$race_case),each=3)
a$Perc=a$Freq/b
"            Var1 Freq       Perc
1   Black_SPTB 0    8 0.09523810
2   Black_SPTB 1   66 0.78571429
3  Black_SPTB NA   10 0.11904762
4   Black_Term 0   25 0.09057971
5   Black_Term 1  215 0.77898551
6  Black_Term NA   36 0.13043478
7   White_SPTB 0   30 0.33707865
8   White_SPTB 1   51 0.57303371
9  White_SPTB NA    8 0.08988764
10  White_Term 0   82 0.21866667
11  White_Term 1  272 0.72533333
12 White_Term NA   21 0.05600000"

fisher.test(matrix(table(paste(meta1$case[meta1$Race=="B" & !is.na(meta1$CIGS_46A)],meta1$CIGS_46A[meta1$Race=="B" &!is.na(meta1$CIGS_46A)])),nrow=2))
"p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.3874121 2.5276880
sample estimates:
odds ratio
1.042288 "
fisher.test(matrix(table(paste(meta1$case[meta1$Race=="W" & !is.na(meta1$CIGS_46A)],meta1$CIGS_46A[meta1$Race=="W" &!is.na(meta1$CIGS_46A)])),nrow=2))
"p-value = 0.01603
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
1.120558 3.350594
sample estimates:
odds ratio
1.947953"


a=data.frame(table(paste(meta1$race_case,meta1$C_PARITY)))
b=c(rep(table(meta1$race_case)[1],3),rep(table(meta1$race_case)[2],3),rep(table(meta1$race_case)[3],2),rep(table(meta1$race_case)[4],3))
a$Perc=a$Freq/b
"            Var1 Freq        Perc
1   Black_SPTB 0   25 0.297619048
2   Black_SPTB 1   58 0.690476190
3  Black_SPTB NA    1 0.011904762
4   Black_Term 0  118 0.427536232
5   Black_Term 1  157 0.568840580
6  Black_Term NA    1 0.003623188
7   White_SPTB 0   35 0.393258427
8   White_SPTB 1   54 0.606741573
9   White_Term 0  168 0.448000000
10  White_Term 1  205 0.546666667
11 White_Term NA    2 0.005333333"

fisher.test(matrix(table(paste(meta1$case[meta1$Race=="B" & !is.na(meta1$C_PARITY)],meta1$C_PARITY[meta1$Race=="B" &!is.na(meta1$C_PARITY)])),nrow=2))
"p-value = 0.04106
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.3240504 0.9963076
sample estimates:
odds ratio
0.5743642 "
fisher.test(matrix(table(paste(meta1$case[meta1$Race=="W" & !is.na(meta1$C_PARITY)],meta1$C_PARITY[meta1$Race=="W" &!is.na(meta1$C_PARITY)])),nrow=2))
"p-value = 0.3441
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.4776224 1.2985747
sample estimates:
odds ratio
0.791295 "

meta1$C_MARITA1=meta1$C_MARITA
meta1$C_MARITA1[meta1$C_MARITA1!=2]=1
a=data.frame(table(paste(meta1$race_case,meta1$C_MARITA1)))
b=c(rep(table(meta1$race_case)[1],3),rep(table(meta1$race_case)[2],3),rep(table(meta1$race_case)[3],2),rep(table(meta1$race_case)[4],2))
a$Perc=a$Freq/b
"            Var1 Freq        Perc
1   Black_SPTB 1   61 0.726190476
2   Black_SPTB 2   22 0.261904762
3  Black_SPTB NA    1 0.011904762
4   Black_Term 1  215 0.778985507
5   Black_Term 2   59 0.213768116
6  Black_Term NA    2 0.007246377
7   White_SPTB 1   34 0.382022472
8   White_SPTB 2   55 0.617977528
9   White_Term 1  117 0.312000000
10  White_Term 2  258 0.688000000"

tab1=table(paste(meta1$case[meta1$Race=="B" & !is.na(meta1$C_MARITA1)],meta1$C_MARITA1[meta1$Race=="B" &!is.na(meta1$C_MARITA1)]))
mat1=matrix(tab1,ncol=2)
fisher.test(mat1)
"p-value = 0.3703
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.4197487 1.4134159
sample estimates:
odds ratio 
0.7614818 "

tab1=table(paste(meta1$case[meta1$Race=="W" & !is.na(meta1$C_MARITA1)],meta1$C_MARITA1[meta1$Race=="W" &!is.na(meta1$C_MARITA1)]))
mat1=matrix(tab1,ncol=2)
fisher.test(mat1)
"p-value = 0.2104
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.8145274 2.2566473
sample estimates:
odds ratio 
1.362284 "


aggregate(meta1$C_POVERT~meta1$race_case,FUN=mean)
"  meta1$race_case meta1$C_POVERT
1      Black_SPTB       121.8451
2      Black_Term       134.5459
3      White_SPTB       268.1467
4      White_Term       310.3499"
aggregate(meta1$C_POVERT~meta1$race_case,FUN=sd)
"  meta1$race_case meta1$C_POVERT
1      Black_SPTB       89.26584
2      Black_Term      109.76557
3      White_SPTB      227.18766
4      White_Term      251.21791"
by(meta1$C_POVERT,meta1$race_case,function(i){length(which(is.na(i)))})
"meta1$race_case: Black_SPTB
[1] 13
------------------------------------------------------------------------------------- 
meta1$race_case: Black_Term
[1] 58
------------------------------------------------------------------------------------- 
meta1$race_case: White_SPTB
[1] 14
------------------------------------------------------------------------------------- 
meta1$race_case: White_Term
[1] 32"

wilcox.test(meta1$C_POVERT[meta1$Race=="B"]~meta1$case[meta1$Race=="B"])
#W = 7491.5, p-value = 0.6862
wilcox.test(meta1$C_POVERT[meta1$Race=="W"]~meta1$case[meta1$Race=="W"])
#W = 11396, p-value = 0.1216

a=data.frame(table(paste(meta1$race_case,meta1$douche)))
b=rep(table(meta1$race_case),each=3)
a$Perc=a$Freq/b
"                     Var1 Freq      Perc
1     Black_SPTB Douching   24 0.2857143
2           Black_SPTB NA   50 0.5952381
3  Black_SPTB No douching   10 0.1190476
4     Black_Term Douching   86 0.3115942
5           Black_Term NA  122 0.4420290
6  Black_Term No douching   68 0.2463768
7     White_SPTB Douching   25 0.2808989
8           White_SPTB NA   38 0.4269663
9  White_SPTB No douching   26 0.2921348
10    White_Term Douching   77 0.2053333
11          White_Term NA  125 0.3333333
12 White_Term No douching  173 0.4613333"

fisher.test(matrix(table(paste(meta1$case[meta1$Race=="B" & !is.na(meta1$douch)],meta1$douch[meta1$Race=="B" &!is.na(meta1$douch)])),nrow=2))
"p-value = 0.1276
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.805367 4.748677
sample estimates:
odds ratio
1.891462 "
fisher.test(matrix(table(paste(meta1$case[meta1$Race=="W" & !is.na(meta1$douch)],meta1$douch[meta1$Race=="W" &!is.na(meta1$douch)])),nrow=2))
"p-value = 0.01499
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
1.115176 4.162106
sample estimates:
odds ratio
2.154406 "

aggregate(meta1$CES_SCOR~meta1$race_case,FUN=mean)
"  meta1$race_case meta1$CES_SCOR
1      Black_SPTB       19.57872
2      Black_Term       18.82011
3      White_SPTB       15.57377
4      White_Term       15.52824"
aggregate(meta1$CES_SCOR~meta1$race_case,FUN=sd)
"meta1$race_case meta1$CES_SCOR
1      Black_SPTB       11.31563
2      Black_Term       10.39400
3      White_SPTB       10.09366
4      White_Term       11.31969"
by(meta1$CES_SCOR,meta1$race_case,function(i){length(which(is.na(i)))})
"meta1$race_case: Black_SPTB
[1] 45
------------------------------------------------------------------------------------- 
meta1$race_case: Black_Term
[1] 99
------------------------------------------------------------------------------------- 
meta1$race_case: White_SPTB
[1] 28
------------------------------------------------------------------------------------- 
meta1$race_case: White_Term
[1] 68"
wilcox.test(meta1$CES_SCOR[meta1$Race=="B"]~meta1$case[meta1$Race=="B"])
#W = 3542.5, p-value = 0.7977
wilcox.test(meta1$CES_SCOR[meta1$Race=="W"]~meta1$case[meta1$Race=="W"])
#W = 9690.5, p-value = 0.6668


aggregate(meta1$LENEG_C~meta1$race_case,FUN=mean)
"  meta1$race_case meta1$LENEG_C
1      Black_SPTB      5.189189
2      Black_Term      3.816568
3      White_SPTB      4.213115
4      White_Term      3.592233"
aggregate(meta1$LENEG_C~meta1$race_case,FUN=sd)
"  meta1$race_case meta1$LENEG_C
1      Black_SPTB      4.551177
2      Black_Term      2.981404
3      White_SPTB      3.596826
4      White_Term      3.181022"
by(meta1$LENEG_C,meta1$race_case,function(i){length(which(is.na(i)))})
"meta1$race_case: Black_SPTB
[1] 47
------------------------------------------------------------------------------------- 
meta1$race_case: Black_Term
[1] 107
------------------------------------------------------------------------------------- 
meta1$race_case: White_SPTB
[1] 28
------------------------------------------------------------------------------------- 
meta1$race_case: White_Term
[1] 66"


wilcox.test(meta1$LENEG_C[meta1$Race=="B"]~meta1$case[meta1$Race=="B"])
#W = 3615.5, p-value = 0.1344
wilcox.test(meta1$LENEG_C[meta1$Race=="W"]~meta1$case[meta1$Race=="W"])
#W = 10292, p-value = 0.2529


meta1$ges_age=meta1$C_DELGS2
meta1$ges_age[meta1$C_DELGS2<32]="<32"
meta1$ges_age[meta1$C_DELGS2>=32 & meta1$C_DELGS2<35 ]="32-34"
meta1$ges_age[meta1$C_DELGS2>=35 & meta1$C_DELGS2<37 ]="35-36"
meta1$ges_age[meta1$C_DELGS2>=37 ]=">=37"

wilcox.test(meta1$C_DELGS2[meta1$race1=="Black"]~meta1$case1[meta1$race1=="Black"])
#W = 0, p-value < 2.2e-16
wilcox.test(meta1$C_DELGS2[meta1$race1=="White"]~meta1$case1[meta1$race1=="White"])
#W = 0, p-value < 2.2e-16

mean(meta1$C_DELGS2[meta1$race_case=="White_Term"])
#39.392
sd(meta1$C_DELGS2[meta1$race_case=="White_Term"])
#1.307502

mean(meta1$C_DELGS2[meta1$race_case=="White_SPTB"])
#34.57303
sd(meta1$C_DELGS2[meta1$race_case=="White_SPTB"])
#1.731313

table(meta1$ges_age[meta1$race_case=="White_SPTB"])
"  <32 32-34 35-36 
    6    23    60 "
table(meta1$ges_age[meta1$race_case=="White_SPTB"])/sum(table(meta1$ges_age[meta1$race_case=="White_SPTB"]))
"       <32      32-34      35-36 
0.06741573 0.25842697 0.67415730 "


wilcox.test(meta1$C_DELGS2[meta1$race_case=="White_Term"],meta1$C_DELGS2[meta1$race_case=="White_SPTB"])
#W = 33375, p-value < 2.2e-16

mean(meta1$C_DELGS2[meta1$race_case=="Black_Term"])
#39.39493
sd(meta1$C_DELGS2[meta1$race_case=="Black_Term"])
#1.375041

mean(meta1$C_DELGS2[meta1$race_case=="Black_SPTB"])
#33.78571
sd(meta1$C_DELGS2[meta1$race_case=="Black_SPTB"])
#2.634669

table(meta1$ges_age[meta1$race_case=="Black_SPTB"])
"<32 32-34 35-36 
   14    24    46"
table(meta1$ges_age[meta1$race_case=="Black_SPTB"])/sum(table(meta1$ges_age[meta1$race_case=="Black_SPTB"]))
"       <32     32-34     35-36 
0.1666667 0.2857143 0.5476190  "


table(meta1$WHYPTB2[meta1$race_case=="White_SPTB"])
"1  2
58 31"
table(meta1$WHYPTB2[meta1$race_case=="White_SPTB"])/sum(table(meta1$WHYPTB2[meta1$race_case=="White_SPTB"]))
"       1         2 
0.6516854 0.3483146 "

table(meta1$WHYPTB2[meta1$race_case=="Black_SPTB"])
"1  2
53 31"
c(53,31)/84
"0.6309524 0.3690476"



#Fig. S2
#read counts analysis
stir_count_melt=read.table(file="/Users/shansun/Google\ Drive/engel/test/test/Engel_16S_stirrups_summary_97_070819.txt",sep="|",quote="")
stir_count_melt1=stir_count_melt[,c(1,2,4)]
stir_count_melt1[,3]=as.numeric(as.character(stir_count_melt1[,3]))
stir_count2=acast(stir_count_melt1, V1~ V2,sum)
stir_count2[is.na(stir_count2)]=0
stir_count2=stir_count2[1:1068,]#remove controls
stir_count3=stir_count2
rownames(stir_count3)=sapply(strsplit(rownames(stir_count3),"\\-"), "[[", 1)

meta1=read.csv(file="/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/meta824_public.csv",row.names=1,header = T)
stir_count5=stir_count3[intersect(rownames(stir_count3),rownames(meta1)),]
count_sum=log10(rowSums(stir_count5))

meta1$reads=count_sum
wilcox.test(count_sum~meta1$race1)
#W = 78697, p-value = 0.1547
wilcox.test(count_sum~meta1$case1)
#W = 52224, p-value = 0.1419
wilcox.test(count_sum~meta1$Douche)
#W = 27226, p-value = 0.168

pdf("read_counts_meta.pdf",width=6,height=7)
par(mfrow=c(2,2))
box1=boxplot(reads~race1,data=meta1,border=c("red","blue"),col="white",main="Wilcoxon test P = 0.155",xlab="Race",ylab="Read counts (log10)",cex.main=0.8,las=1)
stripchart(reads~race1,data=meta1[!meta1$reads %in% box1$out,],vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = c("red","blue"))

box1=boxplot(reads~case1,data=meta1,border=c("red","blue"),col="white",main="Wilcoxon test P = 0.142",xlab="Case",ylab="Read counts (log10)",cex.main=0.8,las=1)
stripchart(reads~case1,data=meta1[!meta1$reads %in% box1$out,],vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = c("red","blue"))

box1=boxplot(reads~douche,data=meta1,border=c("red","blue"),col="white",main="Wilcoxon test P = 0.168",xlab="Douching",ylab="Read counts (log10)",cex.main=0.8,las=1)
stripchart(reads~douche,data=meta1[!meta1$reads %in% box1$out,],vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = c("red","blue"))
dev.off()

#Fig. S3
#rarefy
stir_count_melt=read.table(file="/Users/shansun/Google\ Drive/engel/test/test/Engel_16S_stirrups_summary_97_070819.txt",sep="|",quote="")
stir_count_melt1=stir_count_melt[,c(1,2,4)]
stir_count_melt1[,3]=as.numeric(as.character(stir_count_melt1[,3]))
stir_count2=acast(stir_count_melt1, V1~ V2,sum)
stir_count2[is.na(stir_count2)]=0
stir_count2=stir_count2[1:1068,]

stir_count2=stir_count2[which(rowSums(stir_count2)>1000),]
rownames(stir_count2)=sapply(strsplit(rownames(stir_count2),"\\-"), "[[", 1)

meta1=read.csv(file="/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/meta824.csv",row.names=1,header = T)
stir_norm=stir_count2[intersect(rownames(stir_count2),rownames(meta1)),]
stir_norm1=rrarefy(stir_norm, min(rowSums(stir_norm)))

stir_shan=vegan::diversity(stir_norm1,index = "shannon", MARGIN = 1, base = exp(1))
stir_simp=vegan::diversity(stir_norm1,index = "simpson", MARGIN = 1, base = exp(1))
stir_invsimp=vegan::diversity(stir_norm1,index = "invsimpson", MARGIN = 1, base = exp(1))
meta1$shannon=stir_shan
meta1$simp=stir_simp
meta1$invsimp=stir_invsimp
meta1$lcrisp=log10(stir_norm1[,266]+1)
meta1$liners=log10(stir_norm1[,272]+1)

stir_norm_d1=log10(stir_norm1+1)
gen_pcoa=capscale(stir_norm_d1~1,distance="bray")
var_per=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1:6]*100,2)
pcoa_p=paste("PCoA",c(1:6)," (",var_per,"%)",sep="")

meta1$race_case=factor(meta1$race_case,levels=c("Black_Term","Black_SPTB","White_Term","White_SPTB"))
col_d=c("steelblue","lightblue","tomato","pink")
pdf("new_rarefy_shannon_race_case_box.pdf",width=10,height=6)
ggboxplot(meta1, x = "race_case", y = "shannon", color = "race_case", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race", y = "Shannon Index")+rotate_x_text(45)
dev.off()

pdf("new_rarefy_lcrisp_race_case_box.pdf",width=10,height=6)
ggboxplot(meta1, x = "race_case", y = "lcrisp", color = "race_case", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race", y = "L.crispatus")+rotate_x_text(45)
dev.off()

pdf("new_rarefy_liners_race_case_box.pdf",width=10,height=6)
ggboxplot(meta1, x = "race_case", y = "liners", color = "race_case", add = "jitter", show.legend=TRUE,palette =col_d,
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.title= element_text(size=20),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.key.size = unit(1,"line"), legend.title=element_text(size=20)))+labs( x = "Race", y = "L.iners")+rotate_x_text(45)
dev.off()

wilcox.test(meta1$shannon[meta1$race_case=="White_Term"],meta1$shannon[meta1$race_case=="Black_Term"])
#W = 43956, p-value = 0.001015
wilcox.test(meta1$shannon[meta1$race_case=="White_SPTB"],meta1$shannon[meta1$race_case=="Black_SPTB"])
#W = 2897, p-value = 0.01069
wilcox.test(meta1$shannon[meta1$race_case=="White_Term"],meta1$shannon[meta1$race_case=="White_SPTB"])
#W = 16079, p-value = 0.5929
wilcox.test(meta1$shannon[meta1$race_case=="Black_Term"],meta1$shannon[meta1$race_case=="Black_SPTB"])
#W = 10299, p-value = 0.1217

aldex_rc=list()
levels1=levels(factor(meta1$rc))
k=1
for (i in 1:3){
  for (j in c(2:4)[i:3]){
    level2=levels1[c(i,j)]
    meta3=meta1[meta1$rc%in%level2,]
    stir_norm3=stir_norm1[rownames(meta3),]
    aldex_rc1 <- aldex(t(stir_norm3), meta3$rc, mc.samples=128,test="t",effect=TRUE,include.sample.summary=FALSE,verbose=TRUE)
    aldex_rc[[k]]=aldex_rc1
    k=k+1
  }
}
levels1
#"1_0" "1_1" "2_0" "2_1"
aldex_rc[[1]][c("Lactobacillus_iners","Lactobacillus_crispatus_cluster"),]
# we.ep    we.eBH      wi.ep    wi.eBH
# Lactobacillus_iners             0.15598235 0.8517007 0.15219908 0.8569447
# Lactobacillus_crispatus_cluster 0.02457244 0.6660337 0.04485716 0.7530022
aldex_rc[[2]][c("Lactobacillus_iners","Lactobacillus_crispatus_cluster"),]
# we.ep       we.eBH        wi.ep     wi.eBH
# Lactobacillus_iners             1.812591e-07 4.967335e-05 1.070631e-05 0.00242048
# Lactobacillus_crispatus_cluster 6.810038e-03 1.738327e-01 7.834914e-03 0.18570674

aldex_rc[[5]][c("Lactobacillus_iners","Lactobacillus_crispatus_cluster"),]
# we.ep    we.eBH     wi.ep    wi.eBH
# Lactobacillus_iners             0.1182927 0.7264532 0.4222185 0.8840877
# Lactobacillus_crispatus_cluster 0.2498124 0.7963321 0.1989688 0.7789969

aldex_rc[[6]][c("Lactobacillus_iners","Lactobacillus_crispatus_cluster"),]
# we.ep    we.eBH      wi.ep    wi.eBH
# Lactobacillus_iners             0.88558739 0.9871359 0.77350552 0.9744896
# Lactobacillus_crispatus_cluster 0.08374493 0.8187446 0.08198509 0.8277226

#race
col_d=c("steelblue","tomato")
col2=adjustcolor(col_d[factor(meta1$race1)], alpha.f = 1)
pdf("rarefy_pcoa_race.pdf",width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,1))
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),display="sites",type="none",cex.lab=2,xlab=pcoa_p[1],ylab=pcoa_p[2])
points(pcoa12,"sites",col=adjustcolor(col2, alpha.f = 0.5),pch=16,cex=1.5)
for (n in 1:2){
  ordiellipse(pcoa12, meta1$race1, kind="se", conf=0.95, lwd=4, draw = "lines", col=col_d[n],show.groups=levels(factor(meta1$race1))[n],label=T,font=2,cex=1)
}
legend("topright",c("Black","White"),col=col_d,cex=1.5,pch=16,bty = "n")
dev.off()

#case
col_d=c("purple","turquoise3")
col2=adjustcolor(col_d[factor(meta1$case1)], alpha.f = 1)
pdf("rarefy_pcoa_case.pdf",width=8,height=8)
par(mar=c(5,5,5,5),mfrow=c(1,1))
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),display="sites",type="none",cex.lab=2,xlab=pcoa_p[1],ylab=pcoa_p[2])
points(pcoa12,"sites",col=adjustcolor(col2, alpha.f = 0.5),pch=16,cex=1.5)
for (n in 1:2){
  ordiellipse(pcoa12, meta1$case1, kind="se", conf=0.95, lwd=4, draw = "lines", col=col_d[n],show.groups=levels(factor(meta1$case1))[n],label=T,font=2,cex=1)
}
legend("topright",c("SPTB","Term"),col=col_d,cex=1.5,pch=16,bty = "n")
dev.off()

#Fig. S7
#format table for VALENCIA
stir_count=read.csv(file="/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/stir_norm1.csv",header=T,row.names=1)
crispatus=rowSums(stir_count[,grep("Lactobacillus_crispatus",colnames(stir_count))])
stir_count=stir_count[,-grep("Lactobacillus_crispatus",colnames(stir_count))]
stir_count$Lactobacillus_crispatus=crispatus

colnames(stir_count)[grep("Lactobacillus_gasseri",colnames(stir_count))]="Lactobacillus_gasseri"

jensenii=rowSums(stir_count[,grep("Lactobacillus_jensenii",colnames(stir_count))])
stir_count=stir_count[,-grep("Lactobacillus_jensenii",colnames(stir_count))]
stir_count$Lactobacillus_jensenii=jensenii

colnames(stir_count)[grep("Lachnospiraceae_BVAB1",colnames(stir_count))]="BVAB1"
colnames(stir_count)[grep("Clostridiales_BVAB2",colnames(stir_count))]="BVAB2"
colnames(stir_count)[grep("Clostridiales_BVAB3",colnames(stir_count))]="BVAB3"

Prevotella=rowSums(stir_count[,grep("Prevotella",colnames(stir_count))])
stir_count=stir_count[,-grep("Prevotella",colnames(stir_count))]
stir_count$g_Prevotella=Prevotella

Streptococcus=rowSums(stir_count[,grep("Streptococcus",colnames(stir_count))])
stir_count=stir_count[,-grep("Streptococcus",colnames(stir_count))]
stir_count$g_Streptococcus=Streptococcus

Enterococcus=rowSums(stir_count[,grep("Enterococcus",colnames(stir_count))])
stir_count=stir_count[,-grep("Enterococcus",colnames(stir_count))]
stir_count$g_Enterococcus=Enterococcus

Bifidobacterium=rowSums(stir_count[,grep("Bifidobacterium",colnames(stir_count))])
stir_count=stir_count[,-grep("Bifidobacterium",colnames(stir_count))]
stir_count$g_Bifidobacterium=Bifidobacterium

Staphylococcus=rowSums(stir_count[,grep("Staphylococcus",colnames(stir_count))])
stir_count=stir_count[,-grep("Staphylococcus",colnames(stir_count))]
stir_count$g_Staphylococcus=Staphylococcus

Mycoplasma=rowSums(stir_count[,grep("Mycoplasma",colnames(stir_count))])
stir_count=stir_count[,-grep("Mycoplasma",colnames(stir_count))]
stir_count$g_Mycoplasma=Mycoplasma

Megasphaera=rowSums(stir_count[,grep("Megasphaera",colnames(stir_count))])
stir_count=stir_count[,-grep("Megasphaera",colnames(stir_count))]
stir_count$g_Megasphaera=Megasphaera

colnames(stir_count)[grep("TM7_OTU",colnames(stir_count))]="BVAB_TM7"

ncol1=ncol(stir_count)
sum1=rowSums(stir_count)
stir_count$sampleID=rownames(stir_count)
stir_count$read_count=sum1

stir_count=cbind(stir_count[,(ncol1+1):(ncol1+2)],stir_count[,1:ncol1])
write.csv(stir_count,file="/Users/shansun/Google Drive/engel/new_paper/stir_count_valencia.csv",row.names=F)

valencia=read.csv(file="/Users/shansun/Google Drive/engel/new_paper/valencia_results.csv",header=T,row.names=1)
dim(valencia)
table(valencia$CST)
"I   II   III  IV-A IV-B IV-C    V
190   55  388   47   75   12   57"

table(valencia$CST)/sum(table(valencia$CST))
# I         II        III       IV-A       IV-B       IV-C          V 
# 0.23058252 0.06674757 0.47087379 0.05703883 0.09101942 0.01456311 0.06917476 

table(meta1$cluster_JF_maj1)
# Gardnerella_vaginalis                      Lachnospiraceae_BVAB1 
# 38                                         46 
# Lactobacillus_crispatus_cluster              Lactobacillus_gasseri_cluster 
# 192                                         51 
# Lactobacillus_iners Lactobacillus_jensenii.fornicalis.psittaci 
# 385                                         42 
# Others 
# 70 

valencia1=valencia[,2:268]
valencia$ab_taxa=colnames(valencia1)[apply(valencia1,1,function(i){which(i==max(i))})]
table(valencia$CST)
"   I   II  III IV-A IV-B IV-C    V 
 190   55  388   47   75   12   57"

table(valencia$ab_taxa[valencia$CST=="I"])
"                  BVAB1   Gardnerella_vaginalis Lactobacillus_crispatus
1                       1                     188  "
sort(valencia[which(valencia$ab_taxa!="Lactobacillus_crispatus" & valencia$CST=="I"),"Lactobacillus_crispatus"]/41124.62)

valencia[which(valencia$ab_taxa=="BVAB1" & valencia$CST=="I"),"BVAB1"]/41124.62
#0.3727406
valencia[which(valencia$ab_taxa=="Gardnerella_vaginalis" & valencia$CST=="I"),"Gardnerella_vaginalis"]/41124.62
#0.4847515

table(valencia$ab_taxa[valencia$CST=="II"])
"    Gardnerella_vaginalis Lactobacillus_delbrueckii     Lactobacillus_gasseri       Lactobacillus_iners
2                         1                        51                         1 "
valencia[which(valencia$ab_taxa=="Gardnerella_vaginalis" & valencia$CST=="II"),"Gardnerella_vaginalis"]/41124.62
#0.5988770 0.5082356
valencia[which(valencia$ab_taxa=="Lactobacillus_delbrueckii" & valencia$CST=="II"),"Lactobacillus_delbrueckii"]/41124.62
#0.9306325
valencia[which(valencia$ab_taxa=="Lactobacillus_iners" & valencia$CST=="II"),"Lactobacillus_iners"]/41124.62
#0.5102676

table(valencia$ab_taxa[valencia$CST=="III"])
"        Atopobium_vaginae   Fusobacterium_cluster48             g_Megasphaera              g_Mycoplasma
2                         1                         1                         2
g_Prevotella     Gardnerella_vaginalis   Lactobacillus_crispatus Lactobacillus_delbrueckii
1                         1                         4                         2
Lactobacillus_iners      other     Sneathia_sanguinegens
372                        1                         1 "


sort(valencia[which(valencia$ab_taxa!="Lactobacillus_iners" & valencia$CST=="III"),"Lactobacillus_iners"]/41124.62)
#0.007293772 0.158302996 0.202704818 0.209432444 0.225550676 0.241765567
#0.243421624 0.250111818 0.308714856 0.403896751 0.420060019 0.420204154
#0.428298892 0.467080406 0.473534327 0.490963675

table(valencia$ab_taxa[valencia$CST=="IV-A"])
" BVAB1        g_Mycoplasma Lactobacillus_iners
45                   1                   1 "
sort(valencia[which(valencia$ab_taxa!="BVAB1" & valencia$CST=="IV-A"),"BVAB1"]/41124.62)
#6.673117e-05 1.457231e-01

table(valencia$ab_taxa[valencia$CST=="IV-B"])
" Atopobium_vaginae                BVAB_TM7           g_Megasphaera            g_Mycoplasma            g_Prevotella   Gardnerella_vaginalis
17                       1                       8                       2                       6                      37
Lactobacillus_cluster17     Lactobacillus_iners          Sneathia_amnii
1                       2                       1"


table(valencia$ab_taxa[valencia$CST=="IV-C"])
"Alloscardovia_omnicolens        g_Bifidobacterium          g_Streptococcus
1                        6                        5"

table(valencia$ab_taxa[valencia$CST=="V"])
"Gardnerella_vaginalis    Lactobacillus_iners Lactobacillus_jensenii
2                     11                     44"

sort(valencia[which(valencia$ab_taxa!="Lactobacillus_jensenii" & valencia$CST=="V"),"Lactobacillus_jensenii"]/41124.62)
#0.3755327 0.3953264 0.4181805 0.4317404 0.4339155 0.4465976 0.4525937 0.4580497
#0.4592283 0.4738008 0.4756799 0.4769017 0.4844307

valencia2=valencia[,283:285]
sample_cluster=read.csv(file="/Users/shansun/Google\ Drive/engel/new_paper/0914/JF_cluster_group.csv",row.names = 1,header=T)
cluster2=merge(valencia2,sample_cluster,by=0,all=T)
rownames(cluster2)=cluster2[,1]
table(cluster2$cluster[valencia$CST=="I"])
"         Gardnerella_vaginalis           Lachnospiraceae_BVAB1 Lactobacillus_crispatus_cluster
1                               1                             188 "
table(cluster2$cluster[valencia$CST=="II"])
"Gardnerella_vaginalis     Lactobacillus_delbrueckii Lactobacillus_gasseri_cluster           Lactobacillus_iners
2                             1                            51                             1 "
table(cluster2$cluster[valencia$CST=="III"])
"Atopobium_vaginae         Fusobacterium_cluster48           Gardnerella_vaginalis Lactobacillus_crispatus_cluster       Lactobacillus_delbrueckii
2                               1                               1                               4                               2
Lactobacillus_iners         Megasphaera_OTU70_type1             Mycoplasma_girerdii                          notype             Prevotella_cluster2
369                               1                               2                                    4                               1
Sneathia_sanguinegens
1 "
table(cluster2$cluster[valencia$CST=="V"])
" Gardnerella_vaginalis                        Lactobacillus_iners Lactobacillus_jensenii/fornicalis/psittaci
2                                         13                                         42 "
table(cluster2$cluster[valencia$CST=="IV-A"])
"Lachnospiraceae_BVAB1    Mycoplasma_hominis                notype
45                     1                     1 "
table(cluster2$cluster[valencia$CST=="IV-B"])
"Atopobium_vaginae   Gardnerella_vaginalis Lactobacillus_cluster17     Lactobacillus_iners Megasphaera_OTU70_type1 Megasphaera_OTU71_type2
14                      32                       1                       2                       4                       1
Mycoplasma_girerdii      Mycoplasma_hominis                  notype        Prevotella_amnii          Sneathia_amnii
1                              1                             16                       2                       1 "
table(cluster2$cluster[valencia$CST=="IV-C"])
"Alloscardovia_omnicolens       Bifidobacterium_bifidum Bifidobacterium_breve_cluster      Streptococcus_agalactiae
1                             1                             5                             5 "
cluster3=merge(valencia,sample_cluster,by=0,all=T)
rownames(cluster3)=cluster3[,1]
colnames(cluster3)[286:287]=c("dominate_tax","vagitype")
write.csv(cluster3,file="/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/vagitype_CST_comparison.csv")

stir_norm1=read.csv(file="/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/stir824.csv",header=T,row.names=1)
meta1$Case=c("T","S")[factor(meta1$C_CASE2)]
meta1$Race=c("W","B")[factor(meta1$C_MRACE)]
meta1$Douche=factor(c("D","N")[factor(meta1$douche)],levels=c("D","N"))
meta1$race_case=factor(paste(meta1$race1,meta1$case1,sep="_"))
meta1$Race_Case=factor(paste(meta1$Race,meta1$Case,sep="_"))

stir_shan=vegan::diversity(stir_norm1,index = "shannon", MARGIN = 1, base = exp(1))
stir_simp=vegan::diversity(stir_norm1,index = "simpson", MARGIN = 1, base = exp(1))
stir_invsimp=vegan::diversity(stir_norm1,index = "invsimpson", MARGIN = 1, base = exp(1))
meta1$shannon=stir_shan
meta1$simp=stir_simp
meta1$invsimp=stir_invsimp
meta1$lcrisp=log10(stir_norm1[,177]+1)
meta1$liners=log10(stir_norm1[,183]+1)

meta1$cluster_JF=cluster2$cluster
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


meta1$CST=cluster2$CST
meta1$CST_case=paste(meta1$CST,meta1$C_CASE2)
meta1$CST_race=paste(meta1$CST,meta1$C_MRACE)
meta1$CST_race_case=paste(meta1$CST,meta1$C_CASE2)
meta1$CST=factor(meta1$CST)

# Fig. S7b
pdf("/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/pcoa_CST.pdf",width=8,height=8)
stir_norm_d1=log10(stir_norm1+1)
gen_pcoa=capscale(stir_norm_d1~1,distance="bray")
var_per=round((gen_pcoa$CA$eig/sum(gen_pcoa$CA$eig))[1:6]*100,2)
pcoa_p=paste("PCoA",c(1:6)," (",var_per,"%)",sep="")

col12=c("blue","red","darkgrey","lightgreen","yellow","cyan","orange")
col2=adjustcolor(col12[factor(meta1$CST)], alpha.f = 1)
par(mar=c(5,5,5,5),mfrow=c(1,1))
pcoa12=ordiplot(gen_pcoa,choices=c(1,2),display="sites",type="none",cex.lab=2,xlab=pcoa_p[1],ylab=pcoa_p[2])
points(pcoa12,"sites",col=adjustcolor(col2, alpha.f = 0.5),pch=16,cex=1.5)
# for (i in 1:7){
#   ordiellipse(pcoa12, meta1$CST, kind="se", conf=0.95, lwd=4, draw = "lines", col=col12[i],show.groups=levels(factor(meta1$CST))[i],label=T,font=2,cex=1)
# }
legend("topright",levels(factor(meta1$CST)),col=col12,cex=0.8,pch=16,bty = "n")
dev.off()

# Fig. S7d
cl_case=table(paste(meta1$CST,meta1$case1))
fisher_case=matrix(nrow=7,ncol=14)
for (m in 1:7){
  for (n in 1:7){
    mat1=rbind(cl_case[(m*2-1):(m*2)],cl_case[(n*2-1):(n*2)])
    fisher_case[m,n]=fisher.test(mat1)$estimate
    fisher_case[m,n+7]=fisher.test(mat1)$p.value
  }
}
rownames(fisher_case)=levels(meta1$CST)
colnames(fisher_case)=paste(rep(levels(meta1$CST),2),c(rep("r2",7),rep("P",7)))
write.csv(fisher_case,file="/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/fisher_cluster_case.csv")
#I is significantly lower than III, IV-A, IV-B, V while the rest are not significantly different.
pdf("/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/cluster_perc_bar.pdf",width=8,height=5)
ggplot(data=cluster_com, aes(x=cluster, y=perc_case, fill=cluster)) +geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values=col12)+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=20),
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs( y = "SPTB (%)")+rotate_x_text(0)+ylim(0, 35)
dev.off()

#Fig. S7 panel a
stir_norm1p=stir_norm1/rowSums(stir_norm1)
stir_norm1pm1=stir_norm1p[,apply(stir_norm1p,2,function(i){length(which(i>0.3))})>=5]
Other=1-rowSums(stir_norm1pm1)
stir_norm1pm=cbind(stir_norm1pm1,Other)
stir_norm1pm=stir_norm1pm[order(meta1$CST),]
meta1pm=meta1[order(meta1$CST),]
sortnames=c(rownames(stir_norm1pm[meta1pm$CST=="I",])[order(stir_norm1pm[meta1pm$CST=="I","Lactobacillus_crispatus_cluster"],decreasing = T)],
            rownames(stir_norm1pm[meta1pm$CST=="II",])[order(stir_norm1pm[meta1pm$CST=="II","Lactobacillus_gasseri_cluster"],decreasing = T)],
            rownames(stir_norm1pm[meta1pm$CST=="III",])[order(stir_norm1pm[meta1pm$CST=="III","Lactobacillus_iners"],decreasing = T)],
            rownames(stir_norm1pm[meta1pm$CST=="IV-A",])[order(stir_norm1pm[meta1pm$CST=="IV-A","Lachnospiraceae_BVAB1"],decreasing = T)],
            rownames(stir_norm1pm[meta1pm$CST=="IV-B",])[order(stir_norm1pm[meta1pm$CST=="IV-B","Gardnerella_vaginalis"],decreasing = T)],
            rownames(stir_norm1pm[meta1pm$CST=="IV-C",])[order(stir_norm1pm[meta1pm$CST=="IV-C","Streptococcus_agalactiae"],decreasing = T)],
            rownames(stir_norm1pm[meta1pm$CST=="V",])[order(stir_norm1pm[meta1pm$CST=="V","Lactobacillus_jensenii.fornicalis.psittaci"],decreasing = T)])
stir_norm1pm =stir_norm1pm[sortnames,]
stir_norm1pm=stir_norm1pm[,order(colSums(stir_norm1pm))]

colnames(stir_norm1pm)
col_b=c("hotpink","steelblue","purple","pink","grey","yellow","orange","lightgreen","cyan","red","blue")
pdf(file="/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/class_bar.pdf",height=12,width=12)
par(mfrow=c(2,1),mar=c(3,5,3,3))
barplot(t(stir_norm1pm),col=col_b,space=0, border = NA,axes=F,axisnames = FALSE,ylab="Percentage",cex.lab=2)
axis(side = 2,cex.axis=1)
plot.new()
legend("left",rev(colnames(stir_norm1pm)),fill=rev(col_b),cex=2,bty="n")
dev.off()

#Fig. S7 panel e
shan_tukey=TukeyHSD(aov(stir_shan~meta1$CST))
# diff          lwr        upr     p adj
# II-I       0.25690017  0.027291953  0.4865084 0.0170152
# III-I      0.15428320  0.021502862  0.2870635 0.0111057
# IV-A-I     0.76412472  0.519831755  1.0084177 0.0000000
# IV-B-I     1.03897270  0.834479898  1.2434655 0.0000000
# IV-C-I     0.46241909  0.016073972  0.9087642 0.0366429
# V-I        0.29514219  0.068679442  0.5216049 0.0024129
# III-II    -0.10261696 -0.318673387  0.1134395 0.7998660
# IV-A-II    0.50722456  0.209350735  0.8050984 0.0000123
# IV-B-II    0.78207253  0.515863847  1.0482812 0.0000000
# IV-C-II    0.20551892 -0.272261072  0.6832989 0.8649150
# V-II       0.03824202 -0.245192511  0.3216766 0.9996879
# IV-A-III   0.60984152  0.378239491  0.8414436 0.0000000
# IV-B-III   0.88468950  0.695539421  1.0738396 0.0000000
# IV-C-III   0.30813588 -0.131391598  0.7476634 0.3703431
# V-III      0.14085899 -0.071851666  0.3535696 0.4429764
# IV-B-IV-A  0.27484798 -0.004125476  0.5538214 0.0565750
# IV-C-IV-A -0.30170564 -0.786713722  0.1833024 0.5221701
# V-IV-A    -0.46898253 -0.764438550 -0.1735265 0.0000649
# IV-C-IV-B -0.57655361 -1.042784346 -0.1103229 0.0050757
# V-IV-B    -0.74383051 -1.007331005 -0.4803300 0.0000000
# V-IV-C    -0.16727690 -0.643553251  0.3089995 0.9450599
#alpha diversity of I is significantly lower than the others
pdf("/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/CST_shannon_box.pdf",width=8,height=5)
ggboxplot(meta1, x = "CST", y = "shannon", color = "CST", palette = col12, add = "jitter", lwd = 0.8, xlab = "CST",ylab = "Shannon Index",
          ggtheme = theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
                                     panel.background = element_blank(), axis.line = element_blank(),axis.text = element_text(size=20),title=element_text(size=5),
                                     axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20)))+rotate_x_text(0)
dev.off()

#Fig. S7 panel f
#stacked cluster composition for race
names1=c(names(table(paste(meta1$CST,meta1$Race,meta1$case1)))[1:14],"IV-A W SPTB",names(table(paste(meta1$CST,meta1$Race,meta1$case1)))[15:27])
meta1$CST_race_case=factor(paste(meta1$CST,meta1$Race,meta1$case1),levels=names1)
cluster_com=data.frame(table(meta1$CST_race_case))
cluster_com$perc=cluster_com$Freq/829*100
cluster_com$cluster=sapply(strsplit(as.character(cluster_com$Var1)," "),"[[",1)
cluster_com$race=sapply(strsplit(as.character(cluster_com$Var1)," "),"[[",2)
cluster_com$case=sapply(strsplit(as.character(cluster_com$Var1)," "),"[[",3)
cluster_com$perc_race=cluster_com$Freq/rep(c(360,360,464,464),7)*100
perc_race_case=rep(table(meta1$CST_race_case)[seq(1,28,2)]/(table(meta1$CST_race_case)[seq(1,28,2)]+table(meta1$CST_race_case)[seq(2,28,2)]),each=2)
cluster_com$perc_race_case=perc_race_case*100
perc_case=(cluster_com$Freq[seq(1,28,4)]+cluster_com$Freq[seq(3,28,4)])/by(cluster_com$Freq,cluster_com$cluster,sum)
cluster_com$perc_case=rep(perc_case,each=4)*100

pdf("/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/cluster_race_bar_stack.pdf",width=6,height=7)
ggplot(cluster_com, aes(fill=cluster, y=perc_race, x=race)) +
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values=col12)+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=20),
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs( y = "Percentage (%)")+rotate_x_text(0)+ylim(0, 100)
dev.off()
fisher_cluster_race=matrix(nrow=2,ncol=7)
for (i in 1:7){
  a=cluster_com[seq(1,28,2),2]+cluster_com[seq(2,28,2),2]
  b=c(sum(a[seq(1,14,2)]),sum(a[seq(2,14,2)]))
  fisher_cluster_race[1,i]=fisher.test(matrix(c(a[c(2*i-1,2*i)],b-a[c(2*i-1,2*i)]),nrow=2))$estimate
  fisher_cluster_race[2,i]=fisher.test(matrix(c(a[c(2*i-1,2*i)],b-a[c(2*i-1,2*i)]),nrow=2))$p.value
}
colnames(fisher_cluster_race)=names(table(meta1$CST))
rownames(fisher_cluster_race)=c("OR","P")
write.csv(fisher_cluster_race,file="/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/race_difference_in_CST_composition.csv")
"              I          II          III         IV-A      IV-B       IV-C            V
OR 0.5632643309 0.301265609 2.029074e+00 8.140761e+00 0.6596199 0.25398065 0.3580336800
P  0.0008467345 0.000204307 7.918298e-07 2.562837e-09 0.1125720 0.07761342 0.0008346839
"
#I, II, III, IV-A and V are different between races.


#Fig. S7 panel g
cluster_com1=cluster_com[c(1:4,9:12),]
pdf("/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/cluster_race_perc_bar.pdf",width=5,height=7)
ggplot(data=cluster_com1, aes(x=cluster, y=perc_race_case, fill=race)) +geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values=c("green3","orchid"))+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=20),
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs( x="Vagitype",y = "SPTB (%)")+rotate_x_text(0)+ylim(0, 30)
dev.off()

#overall percentage of SPTB
race_case=data.frame(table(meta1$race_case))
race_case$Percentage=rep(race_case[c(1,3),2]/(race_case[c(1,3),2]+race_case[c(2,4),2])*100,each=2)
race_case$Race=c("Black","Black","White","White")
race_case$type=rep("overall",4)
pdf("/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/race_perc_bar.pdf",width=4,height=7)
ggplot(data=race_case, aes(x=type, y=Percentage, fill=Race)) +geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values=c("green3","orchid"))+
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = 'white', colour = 'black'),axis.text = element_text(size=20),
        axis.title=element_text(size=20),legend.text=element_text(size=20),legend.title=element_text(size=20))+labs( x="",y = "SPTB (%)")+rotate_x_text(0)+ylim(0, 30)
dev.off()

fisher_cluster=matrix(nrow=2,ncol=7)
for (i in 1:7){
  fisher_cluster[1,i]=fisher.test(matrix(cluster_com[c((i*4-3):(i*4)),2],nrow=2))$estimate
  fisher_cluster[2,i]=fisher.test(matrix(cluster_com[c((i*4-3):(i*4)),2],nrow=2))$p.value
}
colnames(fisher_cluster)=names(table(meta1$CST))
rownames(fisher_cluster)=c("OR","P")
write.csv(fisher_cluster,file="/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/race_difference_for_preterm_ratio.csv")
# I        II       III     IV-A      IV-B      IV-C         V
# OR 1.2429470 1.4477166 1.0587509      Inf 0.8327972 6.7081726 2.0942563
# P  0.6471883 0.6889449 0.9037945 0.166456 0.7884191 0.3181818 0.2943632
#Race doesn't modify the associations between CST and sPTB


#Fig. S5
#verify the taxa in Elovitz et al
#Sneathia sanguinegens
#Mobiluncus curtisii/mulieris
#Mageeibacillus indolicus
#g Megasphaera
#Porphyromonas asaccharolytica
#Prevotella buccalis
#g Atopobium

stir_val=read.csv(file="/Users/shansun/Google Drive/engel/new_paper/stir_count_valencia.csv",row.names=1,header=T)
Atopobium=rowSums(stir_val[,grep("Atopobium",colnames(stir_val))])
stir_val=stir_val[,-grep("Atopobium",colnames(stir_val))]
stir_val$g_Atopobium=Atopobium

list_e=c(grep("Sneathia_sanguinegens",colnames(stir_val)),grep("Mobiluncus",colnames(stir_val)),grep("Megasphaera",colnames(stir_val)),
         grep("Porphyromonas_uenonis_asaccharolytica",colnames(stir_val)),grep("Atopobium",colnames(stir_val)))

pdf("/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/Elovitz_taxa_verify.pdf",width=10,height=8,onefile = T)
par(mfrow=c(2,3))
j=1
pmat_t=matrix(nrow=6,ncol=3)
for (i in list_e){
  pmat_t[j,1]=wilcox.test(log10(stir_val[,i]+1)~meta1$Case)$p.value
  pmat_t[j,2]=wilcox.test(log10(stir_val[,i]+1)[meta1$Race=="B"]~meta1$Case[meta1$Race=="B"])$p.value
  pmat_t[j,3]=wilcox.test(log10(stir_val[,i]+1)[meta1$Race=="W"]~meta1$Case[meta1$Race=="W"])$p.value
  main1=paste(colnames(stir_val)[i],"P =",round(pmat_t[j,1],3),"\nW P = ",round(pmat_t[j,2],3),"B P = ",round(pmat_t[j,3],3))
  boxplot(log10(stir_val[,i]+1)~meta1$Race_Case,xlab="race case",ylab="log normalized abundance",main=main1,border=c("red","blue","orange","green"),col="white",cex.lab=1.2,cex.axis=1.2)
  j=j+1
}
dev.off()

#Fig. S6
stir_norm1=read.csv(file="/Users/shansun/Google\ Drive/engel/new_paper/0208NM/revision1/stir824.csv",header=T,row.names=1)
stir_norm1p=stir_norm1/rowSums(stir_norm1)
stir_norm1pm1=stir_norm1p[,apply(stir_norm1p,2,function(i){length(which(i>0.3))})/824>0.01]
Other=1-rowSums(stir_norm1pm1)
stir_norm1pm=cbind(stir_norm1pm1,Other)
stir_norm1pm=stir_norm1pm[order(meta1$cluster_JF_maj1),]
meta1pm=meta1[order(meta1$cluster_JF_maj1),]

type_l=c("Lactobacillus_crispatus_cluster","Lactobacillus_iners","Lactobacillus_gasseri_cluster","Lactobacillus_jensenii.fornicalis.psittaci","Lachnospiraceae_BVAB1","Gardnerella_vaginalis")
sortnames=vector()
for (i in type_l){
  sortnames1=rownames(stir_norm1pm)[which(meta1pm$cluster_JF_maj1==i)][order(stir_norm1pm[which(meta1pm$cluster_JF_maj1==i),i],decreasing = T)]
  print(length(sortnames1))
  sortnames=c(sortnames,sortnames1)
}
sortnames1=rownames(stir_norm1pm)[which(meta1pm$cluster_JF_maj1=="Others")]
sortnames=c(sortnames,sortnames1)

stir_norm1pm =stir_norm1pm[sortnames,]
name_l1=c("Lactobacillus_iners","Lactobacillus_crispatus_cluster",
          "Lactobacillus_jensenii.fornicalis.psittaci","Lactobacillus_gasseri_cluster","Gardnerella_vaginalis",
          "Lachnospiraceae_BVAB1","Other","Megasphaera_OTU70_type1","Atopobium_vaginae")
stir_norm1pm=stir_norm1pm[,rev(match(name_l1,colnames(stir_norm1pm)))]

meta1$douche[is.na(meta1$douche)]="Unknown"
meta1$rcd=paste(meta1$race_case,meta1$douche,sep="_")
colnames(stir_norm1pm)
col_b=c("purple","pink","grey","yellow","orange","lightgreen","cyan","red","blue")
name_l2=c("Lactobacillus_iners","Lactobacillus_crispatus_cluster",
          "Lactobacillus_jensenii\\fornicalis\\psittaci","Lactobacillus_gasseri_cluster","Gardnerella_spp",
          "Lachnospiraceae_BVAB1","Other","Megasphaera_OTU70_type1","Atopobium_vaginae")

pdf(file="/Users/shansun/Google Drive/engel/new_paper/0208NM/revision1/class_bar_stratify.pdf",height=50,width=12,onefile=T)
par(mfrow=c(13,1),mar=c(3,5,3,3))
for (n in names(table(meta1$rcd))){
  stir_norm1pm2=t(stir_norm1pm[which(sortnames%in%rownames(meta1)[meta1$rcd==n]),])
  barplot(stir_norm1pm2,col=col_b,space=0, border = NA,axes=F,axisnames = FALSE,ylab="Percentage",cex.lab=2,main=n)
  axis(side = 2,cex.axis=1)
}
plot.new()
legend("left",name_l2,fill=rev(col_b),cex=2,bty="n")
dev.off()


