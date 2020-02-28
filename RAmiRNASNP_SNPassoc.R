# Prepare Tables and Figures with R script
BiocManager::install("SNPassoc")
library("SNPassoc")
library("scales")


calcOddsRatio <- function(mymatrix,alpha=0.05,referencerow=3){
  mymatrix<-mymatrix+1
  numrow <- nrow(mymatrix)
  myrownames <- rownames(mymatrix)
  OR<-c()
  for (ii in 1:numrow){
    rowname <- myrownames[ii]
    DiseaseUnexposed <- mymatrix[referencerow,1]
    ControlUnexposed <- mymatrix[referencerow,2]
    if (ii != referencerow){ 
      DiseaseExposed <- mymatrix[ii,1]
      ControlExposed <- mymatrix[ii,2] 
      totExposed <- DiseaseExposed + ControlExposed
      totUnexposed <- DiseaseUnexposed + ControlUnexposed
      probDiseaseGivenExposed <- DiseaseExposed/totExposed
      probDiseaseGivenUnexposed <- DiseaseUnexposed/totUnexposed
      probControlGivenExposed <- ControlExposed/totExposed
      probControlGivenUnexposed <- ControlUnexposed/totUnexposed
      pvalue=fisher.test(matrix(c(DiseaseExposed,ControlExposed,DiseaseUnexposed,ControlUnexposed),2,2,byrow=T))$p.value
      # calculate the odds ratio
      oddsRatio <- (probDiseaseGivenExposed*probControlGivenUnexposed)/(probControlGivenExposed*probDiseaseGivenUnexposed)
      # calculate a confidence interval
      confidenceLevel <- (1 - alpha)*100
      sigma <- sqrt((1/DiseaseExposed)+(1/ControlExposed)+(1/DiseaseUnexposed)+(1/ControlUnexposed))
      # sigma is the standard error of our estimate of the log of the odds ratio
      z <- qnorm(1-(alpha/2))
      lowervalue <- oddsRatio * exp(-z * sigma)
      uppervalue <- oddsRatio * exp( z * sigma)
      temp=paste(rownames(mymatrix)[ii]," ",round(oddsRatio,2)," (",round(lowervalue,2),",",round(uppervalue,2),")"," P=",pvalue,sep="")
      print(temp)
      or<-c(round(oddsRatio,2),round(lowervalue,2),round(uppervalue,2),pvalue)
      OR=rbind(OR,or)
    }
  }
  return(OR)
}

epitabe1<-function(data){
  xx<-data.frame(SNP1=data[,grep(snp1,colnames(data))],SNP2=data[,grep(snp2,colnames(data))])
  x1<-subset(xx,SNP1==0 & SNP2==0)
  x2<-subset(xx,SNP1>0 & SNP2==0)
  x3<-subset(xx,SNP1==0 & SNP2>0)
  x4<-subset(xx,SNP1>0 & SNP2>0)
  xxx<-c(nrow(x1),nrow(x2),nrow(x3),nrow(x4))
  return(xxx)
}

epitabe2<-function(data){
  xx<-data.frame(SNP1=data[,grep(snp1,colnames(data))],SNP2=data[,grep(snp2,colnames(data))])
  x1<-subset(xx,SNP1<2 & SNP2 <2)
  x2<-subset(xx,SNP1==2 & SNP2 <2)
  x3<-subset(xx,SNP1<2 & SNP2==2)
  x4<-subset(xx,SNP1==2 & SNP2==2)
  xxx<-c(nrow(x1),nrow(x2),nrow(x3),nrow(x4))
  return(xxx)
}

epitabe3<-function(data){
  xx<-data.frame(SNP1=data[,grep(snp1,colnames(data))],SNP2=data[,grep(snp2,colnames(data))])
  x1<-subset(xx,SNP1==0 & SNP2 ==0)
  x2<-subset(xx,SNP1==1 & SNP2==1)
  x3<-subset(xx,SNP1==2 & SNP2==2)
  xxx<-c(nrow(x1),nrow(x2),nrow(x3))
  return(xxx)
}

percent <- function(x, digits = 1, format = "f", ...){
  result<-paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  return(result)
}

data1<-read.table("ROI.RA3000.counts.assoc",head=T,as.is=T)
data2<-read.table("ROI.RA3000.freq.assoc",head=T,as.is=T)
data1<-data1[order(data1$P,decreasing=F),]
data2<-data2[order(data2$P,decreasing=F),]
data2$F_A=percent(data2$F_A)
data2$F_U=percent(data2$F_U)
data2$A1=paste(data2$A1,"/",data2$A2,sep="")
data2$F_A=paste(data2$F_A,"(",data1$C_A,")",sep="")
data2$F_U=paste(data2$F_U,"(",data1$C_U,")",sep="")
data2$OR=paste(round(data2$OR,2),"(",round(data2$L95,2),"-",round(data2$U95,2),")",sep="")
data2<-data2[,-c(7,8,12,13)]
data2<-data2[,c(1:6,8:9,7)]
write.csv(data2,file="table1.allelic.chisq.table.csv",quote=F,row.names=F)

data<-read.table("ROI.RA3000.model",head=T,as.is=T)
GENO=na.omit(subset(data,TEST=="GENO"))
GENO=GENO[order(GENO$P,decreasing=F),]
write.csv(GENO,file="table2.GENO.chisq.table.csv",quote=F,row.names=F)

data<-read.table("ROI.RA3000.model",head=T,as.is=T)
input=na.omit(subset(data,TEST=="DOM"))
input=input[order(input$P,decreasing=F),]
write.csv(input,file="table3.DOM.chisq.table.csv",quote=F,row.names=F)

data<-read.table("ROI.RA3000.model",head=T,as.is=T)
input=na.omit(subset(data,TEST=="REC"))
input=input[order(input$P,decreasing=F),]
write.csv(input,file="table4.REC.chisq.table.csv",quote=F,row.names=F)

data<-read.table("ROI.RA3000.model",head=T,as.is=T)
input=na.omit(subset(data,TEST=="TREND"))
input=input[order(input$P,decreasing=F),]
write.csv(input,file="table6.TREND.chisq.table.csv",quote=F,row.names=F)

############################################################################
data1<-read.table("ROI.RA3000.counts.assoc",head=T,as.is=T)

x<-head(data1)


frq<-read.table("ROI.RA3000.frq",head=T)
data1<-data1[-na.omit(match(subset(frq,MAF<0.1)$SNP,data1$SNP)),]
p<-p.adjust(data1$P,method = "fdr")
data1$fdr<-p
data1<-data1[order(data1$fdr,decreasing = F),]
head(data1,n=10)

data<-read.table("RA_GWASmeta_Asian_v2.txt",head=T,sep="\t",as.is=T)
head(data)
colnames(data)<-c("SNP","CHR","POS","A1","A2","OR","L95","U95","P")
data$SE=(data$U95-data$L95)/(2*1.96)
write.table(data,file="RA_GWASmeta_Asian_Guo2020.txt",sep="\t",quote=F,col.names = T,row.names = F)



data<-read.table("ROI.RA3000.assoc.logistic",head=T,as.is=T)
x<-head(data)




plink --meta-analysis ROI.RA3000.assoc.logistic RA_GWASmeta_Asian_Guo2020.txt --out meta

data1<-read.table("ROI.RA3000.assoc.logistic",head=T,as.is=T)
data2<-read.table("RA_GWASmeta_Asian_Guo2020.txt",head=T,as.is=T)
data2<-data2[,-5]
SNP<-data1[(data1$SNP %in% data2$SNP),]$SNP
data1<-data1[match(SNP,data1$SNP),]
data2<-data2[match(SNP,data2$SNP),]
head(data2)
write.table(data2,file="RA_GWASmeta_Asian_Guo2020_v3.txt",sep="\t",quote=F,col.names = T,row.names = F)
write.table(data.frame(data2$SNP,data2$A1),file="RA_GWASmeta_Asian_Guo2020.A1.txt",sep="\t",quote=F,col.names = F,row.names = F)

meta<-read.table("plink.meta",head=T)

match(data1$SNP,meta$SNP)

output<-merge(data1,meta,by="SNP")
write.table(int,file="RA_GWASmeta_Asian_Guo2020.miRNA.A1.txt",sep="\t",quote=F,col.names = T,row.names = F)





MyReferenceAllele<-subset(data1,OR<1)[,c(2,7)]
head(data1)
MySigSNPs<-subset(data1,P<0.05)[,c(2)]
write.table(MyReferenceAllele,file="MyReferenceAllele.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(MySigSNPs,file="MySigSNPs.txt",sep="\t",quote=F,col.names=F,row.names=F)

############################################################################
system("bedtools intersect -wao -a ~/hpc/db/hg19/refGene.hg19.V2.bed.txt -b hsa.gff3.hg19.bed | grep -v '\-1' > hsa.gff3.dbSNP153.bed")
data<-read.table("ROI.RA3000.epi.cc",head=T,as.is=T)
miR<-read.table("hsa.gff3.dbSNP153.bed",head=F,as.is=T)
miR1<-miR[match(data$SNP1,miR$V9),5]
miR2<-miR[match(data$SNP2,miR$V9),5]
data<-data.frame(data,miR1,miR2)
write.csv(unique(data),file="table3.epistasis.csv",quote=F,row.names=F)

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/RA/meta3000/MIR")
system("plink --bfile ROI.RA3000.dbsnp --allow-no-sex --reference-allele MyReferenceAllele.txt --recodeA --out ROI.RA3000")

###########################################################################
epi<-read.table("ROI.RA3000.epi.cc",head=T,as.is=T)
head(epi)
chisq<-read.table("ROI.RA3000.counts.assoc",head=T,as.is=T)
head(chisq)
input<-data.frame(epi,chisq[match(epi$SNP1,chisq$SNP),c(10,9)],chisq[match(epi$SNP2,chisq$SNP),c(10,9)])
sub<-input
data<-read.table("ROI.RA3000.raw",head=T,as.is=T)

casen<-c()
contn<-c()
OR<-c()
Name<-c()
for(i in 1:nrow(sub)){
snp1<-sub[i,2]
snp2<-sub[i,4]
name1<-colnames(data)[grep(c(snp1),colnames(data))]
name2<-colnames(data)[grep(c(snp2),colnames(data))]
ele1<-unlist(strsplit(name1,"_"))
ele2<-unlist(strsplit(name2,"_"))
yy1<-"-/-"
yy2<-paste(paste(ele1[2],sep=""),"&",paste("-",sep=""))
yy3<-paste(paste("-",sep=""),"&",paste(ele2[2],sep=""))
yy4<-paste(paste(ele1[2],ele1[2],sep=""),"&",paste(ele2[2],ele2[2],sep=""))
Name<-c(Name,rep(paste(name1,"&",name2),4))
NC<-subset(data,PHENOTYPE=="1")
RA<-subset(data,PHENOTYPE=="2")
mx<-data.frame(RA=epitabe(RA),N=epitabe(NC))
rownames(mx)<-c("11","21","12","22")
#mx<-data.frame(RA=epitabe2(RA),N=epitabe2(NC))
#rownames(mx)<-c("11","12","22")
#mx<-data.frame(RA=epitabe3(RA),N=epitabe3(NC))
#rownames(mx)<-c("11","12","22")
or=calcOddsRatio(mx,referencerow=1)
colnames(or)<-c("OR","L95","U95","P")
or<-rbind(c(1,1,1,1),or)
rownames(or)<-c(yy1,yy2,yy3,yy4)
OR<-rbind(OR,or)
casen<-c(casen,epitabe(RA))
contn<-c(contn,epitabe(NC))
}
matrix<-data.frame(Name,casen,contn,OR)
matrix
write.csv(matrix,file="table4.epistasis.csv",quote=F,row.names=F)

#########################################################################################
#########################################################################################
## 

plink --bfile ROI.RA3000.dbsnp --allow-no-sex --extract MySigSNPs.txt --reference-allele MyReferenceAllele.txt --recodeA --out ROI.RA3000
data1<-read.table("ROI.RA3000.counts.assoc",head=T,as.is=T)
MyReferenceAllele<-subset(data1,OR<1)[,c(2,7)]
MyReferenceAllele
MySigSNPs<-subset(data1,P<0.005)[,c(2)]
write.table(MyReferenceAllele,file="MyReferenceAllele.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(MySigSNPs,file="MySigSNPs.txt",sep="\t",quote=F,col.names=F,row.names=F)


data1<-read.table("ROI.RA3000.counts.assoc",head=T,as.is=T)
MySigSNPs<-subset(data1,P<0.05)[,c(2)]
MySigSNPs
write.table(MySigSNPs,file="MySigSNPs.txt",sep="\t",quote=F,col.names=F,row.names=F)

epi<-read.table("ROI.RA3000.epi.cc",head=T,as.is=T)
chisq<-read.table("ROI.RA3000.counts.assoc",head=T,as.is=T)
data<-read.table("ROI.RA3000.raw",head=T,as.is=T)
data$rowSum<-rowSums(data[,7:ncol(data)])
head(data)

CN<-table(subset(data,PHENOTYPE==1)$rowSum)
RA<-table(subset(data,PHENOTYPE==2)$rowSum)
CN
RA
input<-rbind(RA[1:11],CN[1:11])
rownames(input)<-c("RA","CN")
input<-t(input)
input
OR=calcOddsRatio(input,referencerow=1)
CI95<-c(1,paste(OR[,1],"(",OR[,2],"-",OR[,3],")",sep=""))
P=c(1,OR[,4])
cumulative<-data.frame(input,CI95,P)
write.csv(cumulative,file="table_5.cumulative.csv",quote=F)


#########################################################################################
#########################################################################################
## Table 1. 
table1<-read.csv("table1.allelic.chisq.table.csv")
kg1000<-read.table("1000PopulationFrequency.txt",head=T,sep="\t")
gnomad<-read.table("GnomadFrequency.txt",head=T,sep="\t")
head(gnomad)
EAS<-kg1000[match(table1$SNP,kg1000$dbSNP),]$EAS.Frequency
AFR<-kg1000[match(table1$SNP,kg1000$dbSNP),]$AFR.Frequency
EUR<-kg1000[match(table1$SNP,kg1000$dbSNP),]$EUR.Frequency
out<-data.frame(table1,EAS,EUR,AFR)
write.csv(out,file="table_1.allelic.csv",row.names = F,quote=F)



