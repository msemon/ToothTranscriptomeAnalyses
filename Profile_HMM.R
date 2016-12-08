#This software is governed by the CeCILL license under French law and
#abiding by the rules of distribution of free software.  You can  use, 
#modify and/ or redistribute the software under the terms of the CeCILL
#license as circulated by CEA, CNRS and INRIA at the following URL
#http://www.cecill.info. 

#As a counterpart to the access to the source code and  rights to copy,
#modify and redistribute granted by the license, users are provided only
#with a limited warranty  and the software's author,  the holder of the
#economic rights,  and the successive licensors  have only  limited
#liability. 

#In this respect, the user's attention is drawn to the risks associated
#with loading,  using,  modifying and/or developing or reproducing the
#software by the user in light of its specific status of free software,
#that may mean  that it is complicated to manipulate,  and  that  also
#therefore means  that it is reserved for developers  and  experienced
#professionals having in-depth computer knowledge. Users are therefore
#encouraged to load and test the software's suitability as regards their
#requirements in conditions enabling the security of their systems and/or 
#data to be ensured and,  more generally, to use and operate it in the 
#same conditions as regards security. 

#The fact that you are presently reading this means that you have had
#knowledge of the CeCILL license and that you accept its terms.

# Gaussian Hidden Markov Modeling of absolute differences in gene
# expression between upper and lower teeth, for differentially
# expressed genes.


library(RHmm)
library(DESeq2)
library(plyr)


###############################################
##
## Setting the data for HMM analysis
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


# Getting genes that are differentially expressed between upper and
# lower over the whole time series

design(dds) = ~ tooth
dds <- DESeq(dds)
res <- results(dds, contrast=c("tooth","up","lo"))

thresholdpadj=0.05
res01 <- res[res$padj<thresholdpadj & !is.na(res$padj),]


# Getting genes that are differentially expressed between upper and
# lower at each time point (taking consecutive timepoints as
# replicates)

getDEperPairOfTimePoints <-function(i)
{
colData <- data.frame(std=rep(c(i,i+0.5),each=2),tooth=rep(c("lo","up"),2))
dds <- DESeqDataSetFromMatrix(countData = counts[,timepoints%in%c(i,i+0.5)],
                              colData = colData,
                              design = ~ tooth)
dds <- DESeq(dds)
res <- results(dds, contrast=c("tooth","up","lo"))
res01 <- res[res$padj<thresholdpadj & !is.na(res$padj),]
return(row.names(res01))
}

listeDEPairOfTimePoints=sapply(seq(14.5,17.5,by=0.5),function(time){getDEperPairOfTimePoints(time)})

# Merging the lists of DE genes, either over the whole time series, or
# per time point

listeAllDE=unique(c(unlist(listeDEPairOfTimePoints),row.names(res01)))

# Taking only genes that are in these list of DE genes

datarsdDE=datarsd[listeAllDE,]


######################################################
##
## Making the HMM Analysis on absolute values all DE genes 
##
######################################################

# get the columns for up and low data

x <- row.names(colData)[colData$tooth=="up"]
y <- row.names(colData)[colData$tooth=="lo"]

# We are working on the absolute value of the differential between
# upper and lower tooth, per timepoint.

diffexpr <- data.frame(datarsdDE[,x]-datarsdDE[,y])
names(diffexpr) <- c(paste("stage",seq(14.5,18,by=0.5),sep=""))

# removing genes that have too much difference between upper and lower + too weakly expressed
nb <- apply(counts[listeAllDE,],1,function(x){sum(x>50)})
dif <- apply(diffexpr,1,function(x){sum(x>0.5)})
diffexpr1 <- diffexpr[nb==16&dif==0,]

# Converting to a list for HMMFit
diffexpr2Split <- dlply(diffexpr1,1,function(x){abs(as.numeric(x))})
names(diffexpr2Split) <- row.names(diffexpr1)


# Over all time points, we model these differences through two normal
# distributions, one for high differences ("HI") and one for medium or
# small differences ("ME"), and we model the transitions between these
# states HI and ME in an HMM.


# HMMFit optimizes this modeling on all the data, with an EM
# optimization (Baum-Welch algorithm).

hmm2statesNormal <- HMMFit(obs=diffexpr2Split, nStates=2, dis="NORMAL")

# We check that both states are separate (ie the mean and var of
# Distribution parameters), and state 1 stands for "HI" differences.

hmm2statesNormal$HMM

# Given this optimized HMM, we compute, for each gene, the a
# posteriori probability of each state on each time point. Then, for
# each gene, we set as "HI" all time points for which the a posteriori
# probability of state "HI" is above a given threshold.

threshold <- 0.8

ls <- lapply(diffexpr2Split,function(d){
  fb<- data.frame(forwardBackward(hmm2statesNormal, d)$Gamma)
  names(fb) <- c("HI","ME")
  r <- data.frame(apply(fb,1,rank))
  state <- rep("ME",8)
  state[r[1,]==2&fb[,1]>threshold]="HI"
  return(state)
})
lls <- t(data.frame(ls))
reshmm2statesNormal <- data.frame(lls)
names(reshmm2statesNormal) <- seq(14.5,18,by=0.5)

# We count all types of transitions between and within states "HI" and
# "ME"

ls <- sapply(1:7,function(s){
  table(paste(reshmm2statesNormal[,s],reshmm2statesNormal[,(s+1)],sep="->"))
})
shifts2StatesNormal <- data.frame(ls)
names(shifts2StatesNormal) <- c("14.5->15","15->15.5","15.5->16","16->16.5","16.5->17","17->17.5","17.5->18")


pdf("dotchart_HMMvalabs_DEgenes.pdf", useDingbats=FALSE)
dotchart(as.matrix(shifts2StatesNormal),col=rainbow(4))
dev.off()


pdf("corrplot_HMMvalabs_DEgenes_col.pdf",  useDingbats=FALSE)
corrplot(as.matrix(shifts2StatesNormal),tl.cex=3,title="HMM _ val abs",
         method="circle",is.corr=FALSE,type="full") 
dev.off()





