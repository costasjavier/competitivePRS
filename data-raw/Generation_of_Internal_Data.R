###############################
# Generation of Internal Data #
###############################

#1. Genes hg19, data taken from MAGMA
######################################

hg19<-read.table("~/OldcompetitivePRS/data-raw/NCBI37.3.gene.loc",header = F)
head(hg19)
levels(hg19$V2)
hg19$V2<-gsub(pattern = "X",replacement = "23", x = hg19$V2)
hg19$V2<-gsub(pattern = "Y",replacement = "24", x = hg19$V2)
table(hg19$V2)


#2. Categorizing SNPs by deciles using SNPsnap properties
##########################################################

SNPsnap<-read.table("~/SNPsetPRS/SNPsnap_min.gz",sep="\t",header=T)
SNPsnap$snpid<-SNPsnap$snpID
SNPsnap$snpID<-NULL
SNPsnap<-SNPsnap[,c(5,1,2,3,4)]

summary(SNPsnap)
#Remove the unique SNP with "friends_ld05" = NA
SNPsnap.woNAs<-subset(SNPsnap,!is.na(SNPsnap$friends_ld05))

summary(SNPsnap.woNAs$dist_nearest_gene_snpsnap_protein_coding==Inf)
SNPsnap.limit<-subset(SNPsnap.woNAs,SNPsnap.woNAs$dist_nearest_gene_snpsnap_protein_coding!=Inf)
summary(SNPsnap.limit$dist_nearest_gene_snpsnap_protein_coding)

#Replace Inf by the largest value + 1
SNPsnap.inf<-subset(SNPsnap.woNAs,SNPsnap.woNAs$dist_nearest_gene_snpsnap_protein_coding==Inf)
SNPsnap.inf$dist_nearest_gene_snpsnap_protein_coding<-rep(2762001,dim(SNPsnap.inf)[1])

SNPsnap.corrected<-rbind(SNPsnap.limit,SNPsnap.inf)
summary(SNPsnap.corrected)

#Remove xMHC. We used the limits of de Bakker et al. Nat Genet. 2006 Oct; 38(10): 1166–1172:
#"extended MHC region (defined by the SLC17A2 gene at the telomeric end to the DAXX gene at the centromeric end of chromosome 6)."
#This corresponds to positions 25912984 to 33290793 in hg19.
library(stringr)
splitted<-str_split_fixed(SNPsnap.corrected$snpid,pattern = ":",n = 2)
SNPsnap.corrected.pos<-cbind(splitted,SNPsnap.corrected)
head(SNPsnap.corrected.pos)
SNPsnap.corrected.pos$pos<-as.numeric(levels(SNPsnap.corrected.pos$`2`)[SNPsnap.corrected.pos$`2`])
tail(SNPsnap.corrected.pos)
lox<-SNPsnap.corrected.pos$`1`==6 & SNPsnap.corrected.pos$pos>25912984 & SNPsnap.corrected.pos$pos<33290793
summary(lox)

summary(SNPsnap.corrected.pos$`1`==6)
summary(as.integer(SNPsnap.corrected.pos$`2`)>25912984)
summary(as.integer(SNPsnap.corrected.pos$`2`)<33290793)
xMHC<-subset(SNPsnap.corrected.pos,lox)
head(xMHC)
tail(xMHC)
summary(xMHC$pos)
min(xMHC$pos)
max(xMHC$pos)

SNPsnap.corrected.woMHC<-subset(SNPsnap.corrected.pos,!lox)


#Classification by percentile
##############################

summary(SNPsnap.corrected.woMHC$snp_maf)
SNPsnap.corrected.woMHC$snp_maf.dec<-cut(x=SNPsnap.corrected.woMHC$snp_maf,breaks=quantile(SNPsnap.corrected.woMHC$snp_maf,probs=seq(0,1,by=0.1)),include.lowest=T,labels=F)
table(SNPsnap.corrected.woMHC$snp_maf.dec)


#80.3% of SNPs present low (0-4) gene counts. The classes will be from 0 to 4, and two additional classes, evenly distributed
#Detemination of quartiles
table(SNPsnap.corrected.woMHC$gene_count)
1-((dim(SNPsnap.corrected.woMHC)[1]-1901654-2915025-1462691-819503-516141)/dim(SNPsnap.corrected.woMHC)[1])
gc.more.4<-subset(SNPsnap.corrected.woMHC,SNPsnap.corrected.woMHC$gene_count>4)
quantile(x = gc.more.4$gene_count,probs=c(0,0.5,1))
SNPsnap.corrected.woMHC$gene_count.dec<-cut(x=SNPsnap.corrected.woMHC$gene_count,breaks=c(0,0.1,1,2,3,4,8,140),include.lowest=T,labels=F)
table(SNPsnap.corrected.woMHC$gene_count.dec)

summary(SNPsnap.corrected.woMHC$dist_nearest_gene_snpsnap_protein_coding)
SNPsnap.corrected.woMHC$dist_nearest_gene_snpsnap_protein_coding.dec<-cut(x=SNPsnap.corrected.woMHC$dist_nearest_gene_snpsnap_protein_coding,breaks=quantile(SNPsnap.corrected.woMHC$dist_nearest_gene_snpsnap_protein_coding,probs=seq(0,1,by=0.1)),include.lowest=T,labels=F)
table(SNPsnap.corrected.woMHC$dist_nearest_gene_snpsnap_protein_coding.dec)

summary(SNPsnap.corrected.woMHC$friends_ld05)
SNPsnap.corrected.woMHC$friends_ld05.dec<-cut(x=SNPsnap.corrected.woMHC$friends_ld05,breaks=quantile(SNPsnap.corrected.woMHC$friends_ld05,probs=seq(0,1,by=0.1)),include.lowest=T,labels=F)
table(SNPsnap.corrected.woMHC$friends_ld05.dec)

head(SNPsnap.corrected.woMHC)

#Make subset column according to classification in deciles
colnames(SNPsnap.corrected.woMHC)
SNPsnap.corrected.woMHC$subset<-paste(SNPsnap.corrected.woMHC$snp_maf.dec,paste(SNPsnap.corrected.woMHC$gene_count.dec,paste(SNPsnap.corrected.woMHC$dist_nearest_gene_snpsnap_protein_coding.dec,SNPsnap.corrected.woMHC$friends_ld05.dec,sep=":"),sep=":"),sep=":")
head(SNPsnap.corrected.woMHC)
tail(SNPsnap.corrected.woMHC)
subsets<-as.data.frame(table(SNPsnap.corrected.woMHC$subset))
summary(subsets$Freq)
#There are many classes with very low number of SNPs
#Calculation of percent of SNPs with at least 100 matched SNPs
subsets.less101<-subset(subsets,subsets$Freq<101)
(sum(subsets$Freq)-sum(subsets.less101$Freq))/sum(subsets$Freq)#99.6%
#Calculation of percent of SNPs with at least 1000 matched SNPs
subsets.less1001<-subset(subsets,subsets$Freq<1001)
(sum(subsets$Freq)-sum(subsets.less1001$Freq))/sum(subsets$Freq)#84.86%

#Create a data.frame with the two relevant data
SNPsnap.subsets<-SNPsnap.corrected.woMHC[,c("snpid","subset")]

#3. Create Internal Data
library(devtools)
setwd("R/Rprojects/competitivePRS/")
use_data(hg19,SNPsnap.subsets,internal = T,overwrite = T)
