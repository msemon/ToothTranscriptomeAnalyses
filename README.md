ToothTranscriptomeAnalyses
=========================


# Initial install and configuration
Once you have downloaded the repo to your local environment, the first thing that you want to do is to install all of the packages that this project will require.  This way you will avoid errors when running the scripts.

```r
install.packages(c("reshape", "ggplot2"), dependencies=TRUE)

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

```

=========================

This GitHub repository includes R scripts and data files for conducting the analyses in R as described on the manuscript "Transcriptomic signatures shaped by cell proportions shed light on differences in serial organ morphogenesis".

These scripts were used to construct the Figures 2, 4 and 6 of the manuscript

