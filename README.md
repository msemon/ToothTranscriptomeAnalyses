ToothTranscriptomeAnalyses
=========================


# Initial install and configuration
Once you have downloaded the repository to your local environment, the first thing that you may do is to install all of the packages that this project will require. This way you will avoid errors when running the scripts.

The compressed version of the library RHmm can be found here:

https://cran.r-project.org/src/contrib/Archive/RHmm/RHmm_2.0.3.tar.gz

```r
install.packages(pkgs="RHmm_2.0.3.tar.gz",repos = NULL, type="source")
```

Other necessary packages are :

```r
install.packages(c("reshape", "ggplot2", "ade4", "corrplot"), dependencies=TRUE)

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

```

=========================

This GitHub repository includes R scripts and data files for conducting analyses in R as described on the manuscript "Transcriptomic signatures shaped by cell proportions shed light on differences in serial organ morphogenesis".

These scripts correspond to 

##### A model of cusp maturation (script Model_Cusp_Prop.R)

To estimate the number of cusp patterned for a given RNA-seq sample, we used an embryo harvested from the same litter and having the same weight. Following epithelium dissociation, we performed Fgf4 in situ hybridization to reveal the number and position of Fgf4 expressing enamel knots. We built a simple model for quantifying cusp patterning and maturation in time, fueled with the data on cusp number and position mentioned above. For simplicity, we chose to model cusp proportions with a simple grid model of 6 (lower molar) or 8 (upper molar) equally contributing units, corresponding to the cusps. Once patterned, the territory of each cusp expands at the same speed within a unit. The expansion of the territory is simply coded by increments of one unit (from 0: no cusp patterned to 8: fully expanded cusp) per half day. The proportion of cusp tissue per tooth germ is then computed as the sum of the values assigned to each cusp, divided by the total volume for a fully patterned tooth (total number of cusps x maximal value = 8). In this model, the cusps are still being incremented at the last time point. The fit of these models were computed using generalized linear models (glm in R), which is modeling the proportion of cusp tissue over the total molar tooth germ.

##### A method based on hidden Markov Models to model upper/lower differences (script Profile_HMM.R). 

First, we looked for a description in 2 classes of the overall distribution of upper/lower ratios. To do this, we modeled this distribution as a mixture of 2 normal distributions, and we performed a maximum likelihood estimation of this mixture, from the whole data. Next, we built an hidden Markov model with 2 probabilistic states, one for each of the previous normal distributions : “Small” (resp. “Large”) for the distribution with the smaller (resp. larger) mean. Considering all time series independently with classical forward algorithm, it is possible to compute the likelihood of the whole data, and we optimized this likelihood to fit at best the transition probabilities between the states of the HMM (using R library RHmm version 2.0.3). Finally, from this optimized HMM, we computed the a posteriori probabilities of transition between states in each interval of time points in each time series. In Fig 4b, for each type of transition and each interval of time points, we counted how many genes followed this transition with a probability above 0.8.

NB: Here we model the absolute value of the difference, but of course this is working as well when the sign is kept (in which case more classes are needed)

##### Fitting main trends in upper and lower time-series, by comparison to theoretical profiles. Comparing the trends in upper and lower datasets (sript Profile_Clusters_Comparisons.R).

To classify the genes we first generated theoretical profiles representing all possible combinations of transitions between stages (expression level being either increasing, decreasing or flat between the 8 consecutive stages), making up a total of 2,187 profiles. We then correlated each real expression profile to each theoretical profile using Spearman correlations. We clustered the obtained correlation matrix using K-means  (choosing 10 clusters, and checking that the obtained profiles are repeatable over several iterations of the method). The median of the expression profiles of these 10 clusters were drawn. Similarity of the expression dynamics between upper and lower tooth were estimated by a conservation test at the gene level. For each gene which timeprofile in a given molar best matches a given cluster, we tested whether the timeprofile in the opposite molar would also best match this given cluster. 

