# Estimating the proportion of cusp tissue in the whole tooth germe during morphogenesis, based on the timing of apparition of each cusp

library(DESeq2)
library(ade4)

###############################################
##
## Setting the data for PCA analysis
##
###############################################


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
rsd <- rlogTransformation(dds)
datarsd <- assay(rsd)

# PCA on the whole dataset. 

p <- dudi.pca(t(data.frame(datarsd)),scann=F,nf=10,scale=F,center=T)

pcaplot <- function(p, colData, decal=0.1, nm1="Comp1",  nm2="Comp2",  axis1=1, axis2=2)
{
  perc=round(100* p$eig/sum(p$eig),1)
  plot(p$li[,c(axis1,axis2)],pch=20, col=ifelse(colData$tooth=="up","black","grey"),
       xlab=paste(nm1,perc[axis1],"%"),ylab=paste(nm1,perc[axis2],"%"))
  
  abline(h=0)
  abline(v=0)
  text(p$li[,axis1]-decal, p$li[,axis2]-decal, lab=paste(colData[,1],colData[,2]),offset=1,cex=1)
}

pcaplot(p,colData)

pcaplot(p,colData,nm1="Comp1",  nm2="Comp3",  axis1=1, axis2=3)


###############################################
##
## Setting a list of matrices with the proportion of cusps, based on in situ data. Once patterned, each cusp grows linearly
##
###############################################


# This function is simply making the cusp growing at each step of time

addcusp=function(nx,ny,nr)
{
  sapply(1:3,function(x){sapply(1:nr,function(y){
    if(contrasts[x,y]>0)
    {
      contrasts[x,y]<<-contrasts[x,y]+1
    }})})
  sapply(1:length(unlist(nx)),function(z){
    contrasts[unlist(nx)[z],unlist(ny)[z]]<<-1
  })
  return(contrasts)
}


# Each tooth germ is a grid. The lower tooth has 2 rows of cusps, each containing 3 cusps
# For each of the 8 timepoints, register the coordinates on the grid of the cusp(s) that appear at this timepoint

lnx <- list(2,0,2,0,3,c(1,3),1)
lny <- list(1,0,2,0,2,c(2,1),1)
contrasts <- matrix(0,nrow=3,ncol=2)
# Make the cusps grow of one unit per time point
llow <- lapply(1:8,function(i){addcusp(lnx[i],lny[i],nr=2)})
# llow contains the size and position of cusps at each time point
slow <- lapply(llow,sum)
slowp <- unlist(slow)/(8*6)

# Each tooth germ is a grid. Upper tooth has 3 rows of cusps, each containing 3 cusps (actually just 8 cusps are patterned ultimately)
lnx <- list(2,0,0,c(2,3),0,3,c(3,1),c(2,1))
lny <- list(1,0,0,c(2,1),0,2,c(3,1),c(3,2))
contrasts <- matrix(0,nrow=3,ncol=3)
lup <- lapply(1:8,function(i){addcusp(lnx[i],lny[i],nr=3)})
slup <- lapply(lup,sum)
slupp <- unlist(slup)/(8*8)


# Compare the heterochrony observed with PCA on the first principal component, and the difference of cusp proportion from the model
plot(seq(14.5,18,by=0.5),slowp-slupp,pch=20,col="grey",type="h",ylab="heterochrony cusp prop model")

delta=abs(diff(p$li[,1])[seq(1,15,by=2)])
plot(seq(14.5,18,by=0.5),delta,pch=20,col="grey",type="h",ylab="heterochrony PCA Comp1")


plot(slowp-slupp,delta,pch=20,xlab="deltamodel",ylab="deltaPCA")
text(slowp-slupp,delta,seq(14.5,18,by=0.5))


# Simple model with the number of patterned cusps
clow=c(1,1,2,2,3,5,6,6)
cup=c(1,1,1,3,3,4,6,8)


pdf("maturationsUPLO.pdf")
plot(slowp,slupp,pch=20,xlab="maturation LO",ylab="maturation UP",ylim=c(0,0.6),xlim=c(0,0.6))
abline(b=1,a=0)
text(slowp-slupp,delta,seq(14.5,18,by=0.5))
abline(lm(slupp~slowp),col="red")
dev.off()



nbup=cbind(pos=unlist(slup),neg=c(64-unlist(slup)))
glmup=glm(nbup~p$li[colData$tooth=="up",1],family="binomial")

nblo=cbind(pos=unlist(slow),neg=c(48-unlist(slow)))
glmlo=glm(nblo~p$li[colData$tooth=="lo",1],family="binomial")

nbtot=data.frame(rbind(nbup,nblo))
nbtot$type=c(rep("up",8),rep("lo",8))
nbtot$pca=c(p$li[colData$tooth=="up",1],p$li[colData$tooth=="lo",1])
glmtot=glm(cbind(nbtot$pos,nbtot$neg)~nbtot$pca+nbtot$type,family="binomial")


glmpca=glm(cbind(nbtot$pos,nbtot$neg)~nbtot$pca,family="binomial")
glmnull=glm(cbind(nbtot$pos,nbtot$neg)~1,family="binomial")
coef(summary(glmpca))[,4]
1-logLik(glmpca)/logLik(glmnull)
# 'log Lik.' 0.6967879 (df=2)

pchisq(137.5032-3.0693,df=1)

glmlo=glm(cbind(pos,neg)~pca,family="binomial",data=nbtot[nbtot$type=="lo",])
glmnull=glm(cbind(pos,neg)~1,family="binomial",data=nbtot[nbtot$type=="lo",])
coef(summary(glmpca))[,4] --> 3.460409e-22 
1-logLik(glmpca)/logLik(glmnull)
# 'log Lik.' 0.4154333 (df=2)

1-glmlo$deviance/glmlo$null.deviance
0.9774495

pdf("PCA_modeleGLMlo.pdf")
plot(nbtot$pca[nbtot$type=="lo"],nbtot$pos[nbtot$type=="lo"]/(nbtot$neg[nbtot$type=="lo"]+nbtot$pos[nbtot$type=="lo"]),col=ifelse(nbtot$type[nbtot$type=="lo"]=="up","black","grey"),pch=20,ylim=c(0,0.6),ylab="prop model",xlab="pca axis 1")
x = seq(-40,40)
lines(x,predict(glmlo,newdata=data.frame(pca=x),type="response"),lwd=2,col="yellow")
dev.off()


