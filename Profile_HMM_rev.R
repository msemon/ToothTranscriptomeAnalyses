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
library(corrplot)


###############################################
##
## Setting the data for HMM analysis
##
###############################################

#setwd('~/Documents/PapierSouris/ResoumissionGB/ReplyRev1/')
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
res01 <- res[res$padj<0.1 & !is.na(res$padj),]


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
res01 <- res[res$padj<0.1 & !is.na(res$padj),]
return(row.names(res01))
}

listeDEPairOfTimePoints=sapply(seq(14.5,17.5,by=0.5),function(time){getDEperPairOfTimePoints(time)})

# Merging the lists of DE genes, either over the whole time series, or
# per time point

listeAllDE=unique(c(unlist(listeDEPairOfTimePoints),row.names(res01)))

# Taking only genes that are in these list of DE genes

datarsdDE=datarsd[listeAllDE,]
dim(datarsdDE)
# 2477 genes

######################################################
##
## Making the HMM Analysis on allgenes 
##
######################################################

# get the columns for up and low data

colup <- row.names(colData)[colData$tooth=="up"]
collo <- row.names(colData)[colData$tooth=="lo"]

exprupDE=data.frame(datarsdDE[,colup])
exprloDE=data.frame(datarsdDE[,collo])

nb=apply(counts[listeAllDE,],1,function(x){sum(x>50)})
exprupokDE = exprupDE[nb==16,]
exprlookDE = exprloDE[nb==16,]


#######################################################
####
####  MODELING DIFFERENTIAL BETWEEN UPPPER AND LOWER
####

# We are working on the differential between upper and lower tooth,
# per timepoint.

diffexpr=data.frame(exprupokDE-exprlookDE)
names(diffexpr)=c(paste("stage",seq(14.5,18,by=0.5),sep=""))


##########################################
####
#### 
#
# If we want to filter out genes with a two large difference at any time point
# dif=apply(diffexpr,1,function(x){sum(x>0.5)})
# diffexpr1=diffexpr[dif==0,]
####
####
#######################################



diffexpr1=diffexpr


# Converting to a list for HMMFit
# In this model we take the absolute value of UP/LOW difference 

diffexpr2Split=dlply(diffexpr1,1,function(x){abs(as.numeric(x))})
names(diffexpr2Split)=row.names(diffexpr1)
svg("distriDiffexp.svg")
hist(unlist(diffexpr2Split),nclass=200,main="abs(log(up/low))")
dev.off()
# Over all time points, we model these differences through two normal
# distributions, one for high differences ("HI") and one for medium or
# small differences ("ME"), and we model the transitions between these
# states HI and ME in an HMM.


#######################################
### 1 - Gaussian modelling, 2 states 

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

threshold=0.8

ls=lapply(diffexpr2Split,function(d){
  fb=data.frame(forwardBackward(hmm2statesNormal, d)$Gamma)
  names(fb)=c("HI","ME")
  r=data.frame(apply(fb,1,rank))
  state=rep("ME",8)
  state[r[1,]==2&fb[,1]>threshold]="HI"
  return(state)
})
lls=t(data.frame(ls))
reshmm2statesNormal=data.frame(lls)
names(reshmm2statesNormal)=seq(14.5,18,by=0.5)

# We count all types of transitions between and within states "HI" and
# "ME"

ls=sapply(1:7,function(s){
  table(paste(reshmm2statesNormal[,s],reshmm2statesNormal[,(s+1)],sep="->"))
})

shifts2StatesNormal=data.frame(ls)

names(shifts2StatesNormal)=c("14.5->15","15->15.5","15.5->16","16->16.5","16.5->17","17->17.5","17.5->18")


svg("dotchart_HMMvalabs_DEgenes.svg")
dotchart(as.matrix(shifts2StatesNormal),col=rainbow(4))
dev.off()

# Corrplot representation
svg("corrplot_HMMvalabs_DEgenes.svg")
corrplot(as.matrix(shifts2StatesNormal)[-4,],is.cor=F)
dev.off()


#######################################
### 2 - Mixture modeling

# We perform similar analyses, using an HMM with two hidden states,
# where each hidden state is modeled by a mixture of gaussian
# distributions.

# We choose as the best mixture the one with the lowest AIC.

listhmm2statesM=lapply(2:15,
    function(nbM)
        {
            cat("nb models",nbM,"\r")
            HMMFit(obs=diffexpr2Split, nStates=2, dis="MIXTURE", nMixt=nbM)
        }
                       )

listhmm2states=append(list(hmm2statesNormal),listhmm2statesM)

### the vector of the successive AICs
lAIC=sapply(listhmm2states,function(x){return(x$AIC)})

### The best number of models in the list is
which.min(lAIC)

### Note that this number depends on the random starting points of the
### optimizations. However, over all our observations, the following
### results, which are the same as abode, were always reproduced.

#################
### We look at the model with the best number of gaussian
### distributions.


hmm2statesMB=listhmm2states[[which.min(lAIC)]]
hmm2statesMB$HMM

# We see that state 1 is still the "HI" model, and state 2 the "ME"
# model, since the means of the most probable models of state 1 are
# all above the means of the most probable models of state 2.


threshold=0.8

ls=lapply(diffexpr2Split,function(d){
  fb=data.frame(forwardBackward(hmm2statesMB, d)$Gamma)
  names(fb)=c("HI","ME")
  r=data.frame(apply(fb,1,rank))
  state=rep("ME",8)
  state[r[1,]==2&fb[,1]>threshold]="HI"
  return(state)
})
lls=t(data.frame(ls))
reshmm2statesMB=data.frame(lls)
names(reshmm2statesMB)=seq(14.5,18,by=0.5)

# We count all types of transitions between and within states "HI" and
# "ME"

ls=sapply(1:7,function(s){
  table(paste(reshmm2statesMB[,s],reshmm2statesMB[,(s+1)],sep="->"))
})
shifts2StatesMB=data.frame(ls)
names(shifts2StatesMB)=c("14.5->15","15->15.5","15.5->16","16->16.5","16.5->17","17->17.5","17.5->18")


svg("dotchart_HMMvalabs_DEgenes_MB.svg")
dotchart(as.matrix(shifts2StatesMB),col=rainbow(4))
dev.off()

svg("corrplot_HMMvalabs_DEgenes_MB.svg")
corrplot(as.matrix(shifts2StatesMB),is.cor=F)
dev.off()

pdf("corrplot_HMMvalabs_DEgenes_MB.pdf")
corrplot(as.matrix(shifts2StatesMB),is.cor=F,col="blue")
dev.off()


col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
                           "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))  
col3 <- colorRampPalette(c("red", "white", "blue")) 
col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
                           "cyan", "#007FFF", "blue","#00007F"))   
wb <- c("white","black")
## using these color spectrums
corrplot(as.matrix(shifts2StatesMB),is.cor=F,col="blue")


#######################################
### 3 - HMM with more hidden states

# We perform similar analyses, using an HMM with three hidden states,
# where each hidden state is modeled by a gaussian distribution.

hmm3statesNormal=HMMFit(obs=diffexpr2Split, nStates=3, dis="NORMAL")

#################
### We look at the model

hmm3statesNormal$HMM

# Since the gaussian models are quite separate, we can describe the
# states as one for high differences ("HI"), one for medium
# differences ("ME"), and one for small differences ("SM").

# We compute, for each gene, the a posteriori probability of each
# state on each time point. Then, for each gene, we set as "HI" all
# time points for which the a posteriori probability of state "HI" is
# above a given threshold.

threshold=0.8

ls=lapply(diffexpr2Split,function(d){
  fb=data.frame(forwardBackward(hmm3statesNormal, d)$Gamma)
  names(fb)=c("HI","ME","SM")
  r=data.frame(apply(fb,1,rank))
  state=rep("?",8)
  state[r[1,]==3&fb[,1]>threshold]="HI"
  state[r[2,]==3&fb[,2]>threshold]="ME"
  state[r[3,]==3&fb[,3]>threshold]="SM"
  return(state)
})
lls=t(data.frame(ls))
reshmm3statesNormal=data.frame(lls)
names(reshmm3statesNormal)=seq(14.5,18,by=0.5)


# We count all types of transitions between and within confident states 


ls=sapply(1:7,function(s){
    p=paste(reshmm3statesNormal[,s],reshmm3statesNormal[,(s+1)],sep="->")
    p=c(p,"?->?","?->HI","?->ME","?->SM","HI->?","HI->HI","HI->ME","HI->SM","ME->?","ME->ME","ME->SM","ME->HI","SM->?","SM->ME","SM->HI","SM->SM")
    return(table(p)-1)
})

shifts3StatesNormal=data.frame(ls)
names(shifts3StatesNormal)=c("14.5->15","15->15.5","15.5->16","16->16.5","16.5->17","17->17.5","17.5->18")

shifts3StatesNormal2=shifts3StatesNormal[-grep("\\?",rownames(shifts3StatesNormal)),]

svg("dotchart_3HMMvalabs_DEgenes.svg")
dotchart(as.matrix(shifts3StatesNormal2),col=rainbow(9))
dev.off()




###################################################
###################################################
######## MODELING ON NORMALIZED COUNTS

## We try to retrieve the same kind of signal, but directly modeling
## the levels of expression of both teeth together, instead of the
## differences.

## We keep all genes (even not DE), with a minimum expression level.

exprup=data.frame(datarsd[,colup])
exprlo=data.frame(datarsd[,collo])

nb=apply(counts,1,function(x){sum(x>50)})

exprupok = exprup[nb==16,]
exprlook = exprlo[nb==16,]

####################
## Since we homogenize the gene expressions, so that they can be
## considered together.

expruphom=t(apply(exprupok, 1, function(x){return(x-mean(x))}))
exprlohom=t(apply(exprlook, 1, function(x){return(x-mean(x))}))

vexpruphom=as.vector(expruphom)
vexprlohom=as.vector(exprlohom)


## Even after removing the average expression value per tooth, upper
## and lower expression levels are strongly correlated
## (p-value < 2.2e-16) :

regloup=lm(vexprlohom ~ vexpruphom)

pdf("homogen_expressions_upper_lower.pdf")
plot(vexpruphom, vexprlohom, cex=0.1)
abline(regloup,col="red")
dev.off()


summary(regloup)


##########################
### We compute HMM with 4 hidden states, each with 2D gaussian
### distribution. We hope seeing in the states several ranges of
### differences between upper and lower teeth.


exprloup=lapply(1:dim(exprlohom)[1],
    function(x){
        return(as.data.frame(t(rbind(as.numeric(expruphom[x,]),as.numeric(exprlohom[x,])))))
    })

names(exprloup)=row.names(expruphom)

hmm4loupNormal <- HMMFit(obs=exprloup, nStates=4, dis="NORMAL")

hmm4loupNormal$HMM

### Variances are quite large compared to the differences between
### means, which means that all distributions overlap a lot. However,
### from the differences between up and lo, we could say that 2 models
### can be denoted "HI" (differences around 10e-4) and other models as
### "ME" (differences around 10e-6).

abs(sapply(hmm4loupNormal$HMM$distribution$mean,diff))

rst=rank(abs(sapply(hmm4loupNormal$HMM$distribution$mean,diff)))

hiR=which(rst>=3)

###
threshold=0.8

ls=lapply(exprloup,function(d){
  fb=data.frame(forwardBackward(hmm4loupNormal, d)$Gamma)
  r=data.frame(apply(fb,1,rank))
  state=rep("ME",8)
  sapply(hiR,function(st)
      {
          state[r[st,]>=3&fb[,st]>threshold]<<-"HI"
      }
         )
  return(state)
})

lls=t(data.frame(ls))

reshmm4states=data.frame(lls)

names(reshmm4states)=seq(14.5,18,by=0.5)

# We count all types of transitions between and within states "HI" and
# "ME"

ls=sapply(1:7,function(s){
  table(paste(reshmm4states[,s],reshmm4states[,(s+1)],sep="->"))
})

shifts4States=data.frame(ls)
names(shifts4States)=c("14.5->15","15->15.5","15.5->16","16->16.5","16.5->17","17->17.5","17.5->18")


pdf("dotchart_HMM_4s_genes.pdf")
dotchart(as.matrix(shifts4States),col=rainbow(4))
dev.off()

#################
### Transitions between HI and ME are very rare, and do not show any
### pattern. With this method, the correlation between upper and lower
### totally blurs the signal of their differences.
