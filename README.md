ToothTranscriptomeAnalyses
=========================

This GitHub repository includes R scripts and data files for conducting analyses in R as described on the manuscript "Transcriptomic signatures shaped by cell proportions shed light on differences in serial organ morphogenesis".

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

# These scripts correspond to 

##### A model of cusp maturation (script Model_Cusp_Prop.R)

To estimate the number of cusp patterned for a given RNA-seq sample, we used an embryo harvested from the same litter and having the same weight. Following epithelium dissociation, we performed Fgf4 in situ hybridization to reveal the number and position of Fgf4 expressing enamel knots. We built a simple model for quantifying cusp patterning and maturation in time, fueled with the data on cusp number and position mentioned above. For simplicity, we chose to model cusp proportions with a simple grid model of 6 (lower molar) or 8 (upper molar) equally contributing units, corresponding to the cusps. Once patterned, the territory of each cusp expands at the same speed within a unit. The expansion of the territory is simply coded by increments of one unit (from 0: no cusp patterned to 8: fully expanded cusp) per half day. The proportion of cusp tissue per tooth germ is then computed as the sum of the values assigned to each cusp, divided by the total volume for a fully patterned tooth (total number of cusps x maximal value = 8). In this model, the cusps are still being incremented at the last time point. The fit of these models were computed using generalized linear models (glm in R).

##### A method based on hidden Markov Models to model upper/lower differences (script Profile_HMM.R and Profile_HMM_rev.R). 

We designed a method, based on hidden Markov Models, to track differences in time profile, in a way which explicitly models the temporal series. First, we looked for a description in 2 classes of the overall distribution of the absolute value of upper/lower ratios, taking into account the temporal contraint. To do this, we modeled this distribution as a mixture of 2 normal distributions, organized in an hidden Markov model. This makes 7 free parameters to estimate, using 8 time points from 2477 genes. These 2 distributions (hidden states) correspond to medium (ME) and high (HI) upper/lower difference in expression levels. Considering all time series independently with the classical forward algorithm, it is possible to compute the likelihood of the whole data, and we optimized this likelihood to fit at best the transition probabilities between the states of the HMM the whole data (2477 genes with a differential expression between upper and lower), using R library RHmm version 2.0.3. Finally, from this optimized HMM, we computed the a posteriori probabilities of transition between states for each gene, in each interval of time points in each time series. The scripts produces a corrplot, in which, for each type of transition and each interval of time points, we counted how many genes followed this transition with a probability above 0.8.  in an attempt to better fit the distribution of the ratios, we performed similar analyses, using an HMM with two hidden states, where each hidden state is now modeled by a mixture of gaussian distributions. 

NB: We also modeled the temporal profiles of expression in each tooth (the expression being normalised by the average value for each gene) since we are not interested in the average expression level. We computed a HMM with 4 hidden states, each with a 2 dimension gaussian distribution corresponding to both upper and lower expression levels. We saw that the hidden states correspond to joined distributions that are very similar in upper and lower profiles (either HI in upper and in lower, or LO in upper and in lower). We then continued to estimate the transitions between HI and LO, and found they are very rare, and do not show any pattern. With this method, the correlation between upper and lower profiles totally blurs the signal of their differences. In other words, we can fear that the signal of the upper/lower difference has been hidden by modeling both distributions together, instead of the difference directly. So in our case, it is better to model the upper/lower ratio than the joint value of expression levels in upper and lower. But it is in the script as we thought it might be useful to others.

##### Finding main trends in upper and lower time-series, by comparison to theoretical profiles. Comparing the trends in upper and lower profiles (sript Profile_Clusters_Comparisons.R).

To classify the genes we first generated theoretical profiles representing all possible combinations of transitions between stages (expression level being either increasing, decreasing or flat between the 8 consecutive stages), making up a total of 2,187 profiles. We then correlated each real expression profile to each theoretical profile using Spearman correlations. We clustered the obtained correlation matrix using K-means  (choosing 10 clusters, and checking that the obtained profiles are repeatable over several iterations of the method). The median of the expression profiles of these 10 clusters were drawn. Similarity of the expression dynamics between upper and lower tooth were estimated by a conservation test at the gene level. For each gene which timeprofile in a given molar best matches a given cluster, we tested whether the timeprofile in the opposite molar would also best match this given cluster. 

