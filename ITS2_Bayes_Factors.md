---
title: "Symbiont_and_ITS2_Bayes_Factors "
author: "Matias Gomez"
date: "18/1/2020"
output:md_document


## R Markdown

## Bayes Factor Analysis
## https://richarddmorey.github.io/BayesFactor/#ctables

## Usefull links, code adapted from these links:

## https://www.rdocumentation.org/packages/fmsb/versions/0.7.0/topics/pairwise.fisher.test
## https://www.r-bloggers.com/some-options-for-testing-tables/
## https://rdrr.io/cran/BayesFactor/man/plot.BFBayesFactor.html

## Install package
## BAYES FACTOR


```r
#install.packages("BayesFactor")
library(BayesFactor)
```

## Contingency tables 

## The BayesFactor package implements versions of Gunel and Dickey's (1974) contingency table Bayes factor tests. Bayes factors for contingency tests are computed using the contingencyTableBF function. The necessary arguments are a matrix of cell frequencies and details about the sampling plan that produced the data.


  
## We can perform a Bayes factor analysis using the contingencyTableBF function:
## We used sampleType="indepMulti" and fixedMargin="cols" to specify that the columns are assumed to be sampled as independent multinomials with their total fixed. See the help at ?contingencyTableBF for more details about possible sampling plans and the priors.



## Bayes factor contingency table for Lineage and Site

## 4 sites 2 lineages


```r
library(BayesFactor)
Contingency_Table <- matrix(c(11,1,9,8,4,1,0,0                ),nrow=4,ncol=2,byrow=FALSE,dimnames=list("Site"=c("Colon","Solarte", "Cristobal", "Bastimentos"), "Clades"=c("PAN_1","PAN_2")))

Contingency_Table
```


```r
bf = contingencyTableBF(Contingency_Table, sampleType = "indepMulti", fixedMargin = "cols")
bf
bf_poisson = contingencyTableBF(Contingency_Table, sampleType = "poisson", fixedMargin = "cols")
bf_poisson

bf_jointMulti = contingencyTableBF(Contingency_Table, sampleType = "jointMulti", fixedMargin = "cols")
bf_jointMulti
```


```r
plot(bf, addDenom=T)
plot(bf_poisson, addDenom=T)
plot(bf_jointMulti, addDenom=T)
```

## 4 sites 3 lineages


```r
library(BayesFactor)

Contingency_Table <- matrix(c(6,1,9,8,4,1,0,0,5,0,0,0               ),nrow=4,ncol=3,byrow=FALSE,dimnames=list("Site"=c("Colon","Solarte", "Cristobal", "Bastimentos"), "Clades"=c("PAN_1","PAN_2", "PAN_3")))



Contingency_Table
```


```r
bf = contingencyTableBF(Contingency_Table, sampleType = "indepMulti", fixedMargin = "cols")
bf
bf_poisson = contingencyTableBF(Contingency_Table, sampleType = "poisson", fixedMargin = "cols")
bf_poisson

bf_jointMulti = contingencyTableBF(Contingency_Table, sampleType = "jointMulti", fixedMargin = "cols")
bf_jointMulti
```

```r
plot(bf, addDenom=T)
plot(bf_poisson, addDenom=T)
plot(bf_jointMulti, addDenom=T)
```

## We can also use the posterior function to estimate the difference in probabilities of been tolerant or susceptible between clades, assuming the non-independence alternative:
## For the independent multinomial sampling plan, the chains will contain the individual cell probabilities and the marginal column probabilities. 

## 10000 iterations. Make sure you choose the right distribution for your your sampling, in this case "indepMulti"


```r
bf = contingencyTableBF(Contingency_Table, sampleType = "indepMulti", fixedMargin = "cols")
bf

chains = posterior(bf, iterations = 100000)
```

## We first need to compute the conditional probabilities from the results:
## Here we quantify the probability of tolerance for each clade.



```r
ColonGivenPAN1 = chains[,"pi[1,1]"] / chains[,"pi[*,1]"]
ColonGivenPAN2 = chains[,"pi[1,2]"] / chains[,"pi[*,2]"]
ColonGivenPAN3 = chains[,"pi[1,3]"] / chains[,"pi[*,3]"]
```

## …and then plot the MCMC estimate of the difference between them:


```r
plot(mcmc(ColonGivenPAN1 - ColonGivenPAN2), main = "Increase in probability of being found on Colon between lineages (PAN_1-PAN_2)")
plot(mcmc(ColonGivenPAN1 - ColonGivenPAN3), main = "Increase in probability of being found on Colon between lineages(PAN_1-PAN_3)")

plot(mcmc(ColonGivenPAN2 - ColonGivenPAN1), main = "Increase in probability of being found on Colon between lineages(PAN_2-PAN_1)")
plot(mcmc(ColonGivenPAN2 - ColonGivenPAN3), main = "Increase in probability of being found on Colon between lineages(PAN_2-PAN_3)")


plot(mcmc(ColonGivenPAN3 - ColonGivenPAN1), main = "Increase in probability of being found on Colon between lineages(PAN_3-PAN_1)")
plot(mcmc(ColonGivenPAN3 - ColonGivenPAN2), main = "Increase in probability of being found on Colon between lineages(PAN_3-PAN_2)")
```



```r
SolarteGivenPAN1 = chains[,"pi[2,1]"] / chains[,"pi[*,1]"]
SolarteGivenPAN2 = chains[,"pi[2,2]"] / chains[,"pi[*,2]"]
SolarteGivenPAN3 = chains[,"pi[2,3]"] / chains[,"pi[*,3]"]
```



```r
plot(mcmc(SolarteGivenPAN1 - SolarteGivenPAN2), main = "Increase in probability of being found on Solarte between lineages (PAN_1-PAN_2)")
plot(mcmc(SolarteGivenPAN1 - SolarteGivenPAN3), main = "Increase in probability of being found on Solarte between lineages(PAN_1-PAN_3)")

plot(mcmc(SolarteGivenPAN2 - SolarteGivenPAN1), main = "Increase in probability of being found on Solarte between lineages(PAN_2-PAN_1)")
plot(mcmc(SolarteGivenPAN2 - SolarteGivenPAN3), main = "Increase in probability of being found on Solarte between lineages(PAN_2-PAN_3)")


plot(mcmc(SolarteGivenPAN3 - SolarteGivenPAN1), main = "Increase in probability of being found on Solarte between lineages(PAN_3-PAN_1)")
plot(mcmc(SolarteGivenPAN3 - SolarteGivenPAN2), main = "Increase in probability of being found on Solarte between lineages(PAN_3-PAN_2)")
```


```r
CristobalGivenPAN1 = chains[,"pi[3,1]"] / chains[,"pi[*,1]"]
CristobalGivenPAN2 = chains[,"pi[3,2]"] / chains[,"pi[*,2]"]
CristobalGivenPAN3 = chains[,"pi[3,3]"] / chains[,"pi[*,3]"]
```



```r
plot(mcmc(CristobalGivenPAN1 - CristobalGivenPAN2), main = "Increase in probability of being found on Solarte between lineages (PAN_1-PAN_2)")
plot(mcmc(CristobalGivenPAN1 - CristobalGivenPAN3), main = "Increase in probability of being found on Solarte between lineages(PAN_1-PAN_3)")

plot(mcmc(CristobalGivenPAN2 - CristobalGivenPAN1), main = "Increase in probability of being found on Solarte between lineages(PAN_2-PAN_1)")
plot(mcmc(CristobalGivenPAN2 - CristobalGivenPAN3), main = "Increase in probability of being found on Solarte between lineages(PAN_2-PAN_3)")


plot(mcmc(CristobalGivenPAN3 - CristobalGivenPAN1), main = "Increase in probability of being found on Solarte between lineages(PAN_3-PAN_1)")
plot(mcmc(CristobalGivenPAN3 - CristobalGivenPAN2), main = "Increase in probability of being found on Solarte between lineages(PAN_3-PAN_2)")
```


```r
BastimentosGivenPAN1 = chains[,"pi[4,1]"] / chains[,"pi[*,1]"]
BastimentosGivenPAN2 = chains[,"pi[4,2]"] / chains[,"pi[*,2]"]
BastimentosGivenPAN3 = chains[,"pi[4,3]"] / chains[,"pi[*,3]"]
```



```r
plot(mcmc(BastimentosGivenPAN1 - BastimentosGivenPAN2), main = "Increase in probability of being found on Solarte between lineages (PAN_1-PAN_2)")
plot(mcmc(BastimentosGivenPAN1 - BastimentosGivenPAN3), main = "Increase in probability of being found on Solarte between lineages(PAN_1-PAN_3)")

plot(mcmc(BastimentosGivenPAN2 - BastimentosGivenPAN1), main = "Increase in probability of being found on Solarte between lineages(PAN_2-PAN_1)")
plot(mcmc(BastimentosGivenPAN2 - BastimentosGivenPAN3), main = "Increase in probability of being found on Solarte between lineages(PAN_2-PAN_3)")


plot(mcmc(BastimentosGivenPAN3 - BastimentosGivenPAN1), main = "Increase in probability of being found on Solarte between lineages(PAN_3-PAN_1)")
plot(mcmc(BastimentosGivenPAN3 - BastimentosGivenPAN2), main = "Increase in probability of being found on Solarte between lineages(PAN_3-PAN_2)")
```



## Bayes factor contingency table for Lineage and Symbiont


## 4 genera 3 lineages


```r
library(BayesFactor)
Contingency_Table <- matrix(c(8,8,5,2,3,1,0,1,4,0,0,0),nrow=4,ncol=3,byrow=FALSE,dimnames=list("Site"=c("Symbiodinium","Brevolium", "Cladocopium", "Durusdinium"), "Clades"=c("PAN_1","PAN_2", "PAN_3")))

Contingency_Table
```



```{ r eval=FALSE}
fisher<- fisher.test(Contingency_Table, simulate.p.value = TRUE, B=1000000)
fisher

bf = contingencyTableBF(Contingency_Table, sampleType = "indepMulti", fixedMargin = "cols")
bf
bf_poisson = contingencyTableBF(Contingency_Table, sampleType = "poisson", fixedMargin = "cols")
bf_poisson

bf_jointMulti = contingencyTableBF(Contingency_Table, sampleType = "jointMulti", fixedMargin = "cols")
bf_jointMulti
```


```r
plot(bf, addDenom=T)
plot(bf_poisson, addDenom=T)
plot(bf_jointMulti, addDenom=T)
```

