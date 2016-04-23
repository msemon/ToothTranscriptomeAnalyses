# http://www.r-bloggers.com/rstudio-and-github/

## R Code to reproduce Figure 2b et c

setwd("/home/msemon/Documents/PapierSouris/RSCRIPT_FOR_GITHUB/")

library(DESeq2)
library(plyr)

# Reading the data, using DESEq2 to get normalised (rld) values
counts <- read.table("TabCountsMouse",sep="\t",h=T,row.names="id")
counts <- counts[1:(nrow(counts)-5),]

timepoints <- rep(seq(14.5,18,by=0.5),each=2)
tooth <- rep(c("lo","up"),8)

colData <- data.frame(std=timepoints,tooth=tooth)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ 1)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
vsd <- varianceStabilizingTransformation(dds)
rsd <- rlogTransformation(dds)

datarsd <- assay(rsd)

# This function simply produces all 2187 possibilities obtained by increasing (1), decreasing (-1) or keeping stable (0) the expression level 
# from one step to the other

f=function(x){
    b <- sapply(1:7,function(d){
    a <- x%%3
    x <<- x%/%3
    return(a-1)
  })
}

code=data.frame(t(sapply(1:3^7,f)))

# The resulting expression theoretical expression profiles can be obtained by cumulative sum. Expression levels is set at 0 for
# the first time point. Flat profile is removed
codet=data.frame(t(apply(code,1,function(x){c(0,cumsum(x))})))
codet=codet[-1093,]

## This function computes the correlation between normalized expression levels and the theoretical profiles.
# Then these correlations are clustered using K-mean function, to obtain the main tendencies in the data.

getProfiles=function(toothtype="up",dat=datarsd,nb=10)
{
x <- row.names(colData)[colData$tooth==toothtype]
dat1 <- data.frame(t(apply(dat[,x],1,function(s){s-s[1]})))

# to remove one quartile of genes that vary the least during the timecourse
maxmin <- apply(dat[,x],1,function(x){max(x)-min(x)})
dat1 <- dat1[maxmin>quantile(maxmin)[2],]

## correlate each true gene profile with all theoretical profiles
corGene <- apply(dat1,1,function(x){apply(codet,1,function(y){cor(x,y,method="spearman")})})
corGene <- t(corGene)

# get the main profiles using the kmeans methods
fitClus <- kmeans(na.omit(corGene), nb) 

# plot the main profiles
profils <- aggregate(dat1,by=list(fitClus$cluster),FUN=median)
colo=rainbow(10)
plot(1:8,profils[1,-1],type="l",ylim=c(min(profils[,-1]),max(profils[,-1])),col=colo[1],xlab="time point",ylab="expr. profile")
sapply(2:10,function(i){lines(1:8,profils[i,-1],type="l",col=colo[i])})
legend("topleft",col=colo,pch=20,legend=paste(1:10," (",fitClus$size,")"),ncol=2,cex=0.8)


return(fitClus)
}

## this step may take a while: Hence, we run it for demo on only a sample of 1000 genes

s=sample(1:nrow(datarsd),1000)
fitClusUp <- getProfiles(toothtype="up",dat=datarsd[s,],nb=10)
fitClusLo <- getProfiles(toothtype="lo",dat=datarsd[s,],nb=10)


## This function computes the correspondance between expression profiles in both tooth types

getConservationProfiles=function(toothtype1="up",toothtype2="lo",dat=datarsd,fitClus1=fitClusUp,nb=10)
{
  # first, get the genes for which the expression profiles in toothtype1 is well assigned to its cluster (cor > 0.7)
  x <- row.names(colData)[colData$tooth==toothtype1]
  dat1 <- data.frame(t(apply(dat[,x],1,function(s){s-s[1]})))
  dat1=dat1[names(fitClus1$cluster),]
  profils1 <- aggregate(dat1,by=list(fitClus1$cluster),FUN=median)
  n=1:nb
  assignToothtype1OwnCluster=apply(dat1,1,function(x){apply(profils1[,-1],1,function(j){cor.test(x,j,method="spearman")$estimate})})
  assignToothtype1OwnCluster1=apply(assignToothtype1OwnCluster,2,function(x){z=NA;if(max(x)>0.7){z=n[x==max(x)][1]};return(z)})
  
  # Then, get the genes for which the expression profile in toothtype2 is well assigned to one cluster as defined in toothtype1
  y <- row.names(colData)[colData$tooth==toothtype2]
  dat2 <- data.frame(t(apply(dat[,y],1,function(s){s-s[1]})))

  clusterProjToothtype2OnType1=apply(dat2,1,function(x){apply(profils1[,-1],1,function(j){cor(x,j,method="spearman")})})
  clusterProjToothtype2OnType1b=unlist(apply(clusterProjToothtype2OnType1,2,function(x){z=NA;
  if(!is.na(x[1])&max(x)>0.7){z=n[x==max(x)][1]};return(z)}))
  
  mergeclusters=merge(fitClus1$cluster,assignToothtype1OwnCluster1,by.x=0,by.y=0)
  mergeclusters1=merge(mergeclusters,clusterProjToothtype2OnType1b,by.x="Row.names",by.y=0)
  names(mergeclusters1)=c("id","cluster.tooth1","reassigned.tooth1","reassigned.tooth2")

  mergeclusters2 <- mergeclusters1[mergeclusters1$cluster.tooth1==mergeclusters1$reassigned.tooth1&!is.na(mergeclusters1$reassigned.tooth1),]
  
bw=t(table(mergeclusters2$cluster.tooth1,mergeclusters2$reassigned.tooth2,useNA="always"))
bwt=sapply(1:10,function(x){bw[x,x]})
bwt=rbind(bwt,apply(bw[1:10,1:10],2,sum)-bwt)
bwt=rbind(bwt,bw[11,1:10])
barplot(bwt,col=c("black","dark grey","light grey"))
legend(x="topright",col=c("black","dark grey","light grey"),legend=c("Same cluster","Other cluster","No best match"),pch=20)
}


getConservationProfiles(toothtype1="up",toothtype2="lo",dat=datarsd[s,],fitClus1=fitClusUp,nb=10)
getConservationProfiles(toothtype1="lo",toothtype2="up",dat=datarsd[s,],fitClus1=fitClusLo,nb=10)


