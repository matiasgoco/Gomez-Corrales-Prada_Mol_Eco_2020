---
title: "SNP filtering, Clone correction, DAPC, Structure, NJ tree"
author: "Matias Gomez"
date: "5/22/2020"
output: md_document
---

## Reanalysis of Dziedzic et al. (2019) O. faveolata data (https://www.ncbi.nlm.nih.gov//bioproject/PRJNA413258, https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=2&WebEnv=NCID_1_14569396_130.14.18.48_5555_1574118114_1163998098_0MetA0_S_HStore&o=libraryselection_s%3Aa&s=SRR8844779,SRR8844780,SRR8844781,SRR8844782,SRR8844783,SRR8844784,SRR8844785,SRR8844786,SRR8844787,SRR8844788,SRR8844789,SRR8844790,SRR8844792,SRR8844793,SRR8844794,SRR8844795,SRR8844796,SRR8844797,SRR8844798,SRR8844799,SRR8844800,SRR8844801,SRR8844802,SRR8844803,SRR8844804,SRR8844805,SRR8844808,SRR8844809,SRR8844810,SRR8844812,SRR8844807,SRR8844806,SRR8844791,SRR8844776,SRR8844811,SRR8844778,SRR8844777,SRR8844775,SRR8844774)

## Activate previously created environment Orbicella

```bash
conda activate Orbicella
```

## SNP filtering (http://www.ddocent.com/filtering/)
## First, let´s create a working dir and change to it, then copy the raw VCF file generetated by dDocent to start with the SNP filtering


```bash
mkdir SNP_Filtering_Dziedzic_2019
cd SNP_Filtering_Dziedzic_2019
cp /RAID_STORAGE2/mgomez/Orbicella_Raw/Dziedzic_2019_Data/TotalRawSNPs.vcf .
```

## Adapted from dDdocent pipeline (http://www.ddocent.com/filtering/)

## Typically we are only interested in polymorphic loci (SNPs). We apply this filter first, since most loci are monomorphic and excluding them will greatly reduce file sizes. This is accomplished by running the following script: PolyFilter.pl -i combined.tab -n 2 -p y -o snps.tab .In this example, we selected all loci at which 2 or more genotypes were observed, writing them to a file called "snps.tab". 

## I am going to keep only biallelic SNPs first


```bash
vcftools --vcf TotalRawSNPs.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out TotalRawSNPs_BiAllelic

After filtering, kept 39 out of 39 Individuals
Outputting VCF file...
After filtering, kept 370231 out of a possible 383160 Sites
Run Time = 48.00 seconds
```

## I need to keep loci with at least 2 o or more genotypes. I can also keep loci only genotyped in 90% of individuals

## we are going to use the program VCFtools (http://vcftools.sourceforge.net) to filter our vcf file. To make this file more manageable, let’s start by applying three step filter. We are going to only keep variants that have been successfully genotyped in 90% of individuals, a minimum quality score of 30, and a minor allele count of 3
## In this code, we call vcftools, feed it a vcf file after the --vcf flag, --max-missing 0.9 tells it to filter genotypes called below 90% (across all individuals) the --mac 2 flag tells it to filter SNPs that have a minor allele count less than 3. This is relative to genotypes, so it has to be called in at least 1 homozygote and 1 heterozygote or 3 heterozygotes. The --recode flag tells the program to write a new vcf file with the filters, --recode-INFO-all keeps all the INFO flags from the old vcf file in the new one. Lastly, --out designates the name of the output 


```bash
vcftools --gzvcf TotalRawSNPs_BiAllelic.recode.vcf --max-missing 0.9 --mac 3 --minQ 30 --recode --recode-INFO-all --out TotalRawSNPs_BiAllelic_mac3_miss90

After filtering, kept 39 out of 39 Individuals
Outputting VCF file...
After filtering, kept 7037 out of a possible 370231 Sites
Run Time = 4.00 seconds
```

## Now keep loci with a Minor Allele Freq of 2.6 or present in at least 1 individual (1*100/39= 2.57) and a mininmum depth of 5x


```bash
vcftools --vcf TotalRawSNPs_BiAllelic_mac3_miss90.recode.vcf  --maf 0.026 --min-meanDP 20 --recode --recode-INFO-all --out TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20 After filtering, kept 39 out of 39 Individuals
Outputting VCF file...
After filtering, kept 6007 out of a possible 7037 Sites
Run Time = 1.00 seconds 
```

## The next step is to get rid of individuals that did not sequence well. We can do this by assessing individual levels of missing data.


```bash
vcftools --vcf TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.recode.vcf --missing-indv

```

## This will create an output called out.imiss. Let’s examine it.

```bash
cat out.imiss
```

## We can see that one  have as high as 30% missing data. Let’s take a look at a histogram


```bash
mawk '!/IN/' out.imiss | cut -f5 > totalmissing
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
```

## There is only one individuals with 30 % missing data. Now we need to create a list of individuals with more than 30 missing data. 

```bash
mawk '$5 > 0.30581' out.imiss | cut -f1 > lowDP.indv
cat lowDP.indv
```

## Since it is only one individual with 30% missing information I decided to keep it for the time being.

## First, we need to create file that has two tab separated columns. First with the sample name, second with population assignment.


```bash
vcf-query -l TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.recode.vcf > TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.recode.vcf.indv
cut -d "_" -f1 TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.recode.vcf.indv > TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.recode.vcf.indv.pop
paste TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.recode.vcf.indv TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.recode.vcf.indv.pop > popmap.txt
```

## You need to add headers the table before you uploaded to R: edit the table in R by clicking on the table and addding the headers separated by tab


```r
SNP_6007.pop <-read.table("/RAID_STORAGE2/mgomez/Orbicella_Raw/Dziedzic_2019_Data/SNP_Filtering_Dziedzic_2019/popmap.txt", header=TRUE)
SNP_6007.pop
```

## Convert the VCF file to BayeScan format using pgdspider


```bash
java -jar /usr/local/bin/PGDSpider2-cli.jar -inputfile TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.recode.vcf  -outputfile TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.BayesScan -spid BSsnp.spid
```

## Now scan for outlier loci with BayeScan (in case Bayescan does not run, try loading the vcf file to R and write a hierfstat file an then a bayescan one)


```bash
BayeScan2.1_linux64bits TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.BayesScan -o TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.BayesScan.Output -n 5000 -thin 10 -nbp 20    -pilot 5000   -burn 50000 -pr_odds 100 -threads 20
```


```bash
cp /RAID_STORAGE2/mgomez/Orbicella_Raw/Libraries_compiled/SNP_Filtering/plot_R.r .
```


```r
source("plot_R.r")
BayeScan_Results<-plot_bayescan("TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.BayesScan.Output_fst.txt",FDR=0.05)
BayeScan_Results$outliers
BayeScan_Results$nb_outliers
```


```r
BayeScan_Results_sel=read.table("TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.BayesScan.Output.sel",colClasses="numeric")
parameter="Fst1"
plot(density(BayeScan_Results_sel[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))

parameter="logL"
plot(density(BayeScan_Results_sel[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))

BayeScan_Results_Fst=read.table("TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.BayesScan.Output_fst.txt",colClasses="numeric")
parameter="alpha"

plot(density(BayeScan_Results_Fst[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))
```

## if you have the package "boa" installed, you can very easily obtain Highest Probability 
## Density Interval (HPDI) for your parameter of interest (example for the 95% interval):
## > boa.hpd(mydata[[parameter]],0.05)


```r
#install.packages("boa")
library(boa)
parameter="Fst1"
boa.hpd(BayeScan_Results_sel[[parameter]],0.05)
```

## Remove Indels and leave SNPs only 


```bash
vcftools --vcf TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20.recode.vcf --remove-indels --recode --recode-INFO-all --out TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels

After filtering, kept 39 out of 39 Individuals
Outputting VCF file...
After filtering, kept 5354 out of a possible 6007 Sites
```

## We now need to keep SNPs that are not in LD, we want to keep SNPs that are at least 1000 pb apart 


```bash
vcftools --vcf TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels.recode.vcf --thin 1000 --recode --recode-INFO-all --out TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000 

Outputting VCF file...
After filtering, kept 3560 out of a possible 5354 Sites
```


## Clone correction


```r
library(vcfR)
VCF_3560_SNPs<-read.vcfR("TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000.recode.vcf")
VCF_3560_SNPs

Genind_3560_SNPs <- vcfR2genind(VCF_3560_SNPs)
Genind_3560_SNPs
```



```r
library(poppr)
filtered_Genind_3560_SNPs<-filter_stats(Genind_3560_SNPs, distance = bitwise.dist, threshold = 1e+06 +.Machine$double.eps^0.5, 
             stats = "All", missing = "ignore",plot = TRUE,cols = NULL, hist = "scott")

print(farthest_thresh <- cutoff_predictor(filtered_Genind_3560_SNPs$farthest$THRESHOLDS))
print(nearest_thresh <- cutoff_predictor(filtered_Genind_3560_SNPs$nearest$THRESHOLDS))
print(average_thresh <- cutoff_predictor(filtered_Genind_3560_SNPs$average$THRESHOLDS))
```



```r
Genclone_3560_SNPs<- as.genclone(Genind_3560_SNPs)

mlg.filter(Genclone_3560_SNPs, distance = "bitwise.dist", algorithm = "nearest_neighbor")<-0.015
mlg.table(Genclone_3560_SNPs)
mll(Genclone_3560_SNPs)
mlg.id(Genclone_3560_SNPs)
```


## Remove putative clones from VCF and rerun all analyses


```bash
nano clones2remove.txt

PAN_Ofav29
PAN_Ofav20
PAN_Ofav15
PAN_Ofav13

```


```bash
vcftools --vcf TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000.recode.vcf --remove clones2remove.txt --recode --recode-INFO-all --out TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones

After filtering, kept 35 out of 39 Individuals
Outputting VCF file...
After filtering, kept 3560 out of a possible 3560 Sites
Run Time = 0.00 seconds

```


## Make a pop file for the next steps 


```bash
vcf-query -l TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones.recode.vcf > TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones.indv
cut -d "_" -f1 TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones.indv > TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones.indv.pop
paste TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones.indv TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones.indv.pop > TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones_pop.txt

less TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones_pop.txt
```



## You need to add headers the table before you uploaded to R: edit the table in R by clicking on the table and addding the headers separated by tab. I made all individuals PAN 1


```r
SNP_3560.pop <-read.table("TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones_pop.txt", header=TRUE)
SNP_3560.pop 
```

## Edit file so pop is a number and not a string. This needs to be done in order to run structure since it will not take string in the population column. 
## I added a the column Pop_No and assigned each individual a 1


```r
SNP_3560.pop <-read.table("TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones_pop.txt", header=TRUE)
SNP_3560.pop
```

## Load VCF on R and convert it to Structure format


```r
library(vcfR)
VCF_3560_SNPs<-read.vcfR("TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones.recode.vcf")
vcf_3560_SNPs 
```

## Now that we have a VCF thorougly filtered with unlinked, without indels and no outlier loci we can run Bayesian clustering as implemented in STRUCTURE software

## Running structure

## https://radcamp.github.io/NYC2019/05_STRUCTURE_API.html

## Although there are many newer and faster implementations of STRUCTURE, such as faststructure or admixture, the original STRUCTURE works much better with missing data, which is of course a common feature of RAD-seq data sets.

## We convert a VCF file to the structure format with the following steps. In this way a make the needed file with the pop labels but without the individuals labels. Even though this file does not have the individual`s labels, these are optinal for Structure so we can try running structure with this file to seee how it goes. 

## We start by loading the following R packages


```r
library(adegenet)
library(hierfstat)
SNP_3560.genind <- vcfR2genind(VCF_3560_SNPs)
pop(SNP_3560.genind) <- SNP_3560.pop$Pop_No
indNames(SNP_3560.genind)<-SNP_3560.pop$Sample
indNames(SNP_3560.genind)
SNP_3560_ID<-SNP_3560.pop$Sample
SNP_3560.genind
SNP_3560.hierfstat <- genind2hierfstat(SNP_3560.genind)
SNP_3560.hierfstat

write.struct(SNP_3560.hierfstat, fname = "SNP_3560_STRUCTURE.str") 
```




```bash
which python 
which structure
conda install -c bioconda structure
wget https://www.crypticlineage.net/assets/strauto_1.tar.gz

conda install -c conda-forge parallel

cp ../strauto_1.tar.gz .

tar xvzf strauto_1.tar.gz

cp ../structure_3rd_run/extraparams .
cp ../structure_3rd_run/mainparams .
cp ../structure_3rd_run/input.py .
cp ../structure_3rd_run/structureCommands .

cp ../structure_3rd_run/structure_3rd_run.sh .
mv structure_3rd_run.sh structure_revision_run.sh


sbatch structure_revision.sh

```

## You need to edit the mainparams file to specify 35 individulas 

## You need to edit the extraparams file to ALPHA=0.1


## Since the previous run was successful, I can now run my own data.
## First we modify the file input.py according to my own data.
## We use  k = 1-8, 
## Analyses were performed with five runs of 500,000 iterations each (250,000 burn‐in), with correlated allele frequencies and under the admixture model (Falush, Stephens, & Pritchard, 2003). (Kornillos et al, 2019)

## 1. How many populations are you assuming? [Integers]
## maxpops = 8

## 2. How many burnin you wish to do before collecting data [Integers]
## burnin = 250000

## 3. How long do you wish to collect the data after burnin [Integers]
## mcmc = 500000

## 4. Name of your dataset.  Don't remove quotes. No spaces allowed. Exclude the '.str' extension.  
##    e.g. dataset = "sim" if your datafile is called 'sim.str'
## dataset = "SNP_3560_STRUCTURE"

## 5. How many runs of Structure do you wish to do for every assumed cluster K? [Integers]
## kruns = 10

## 6. Number of individuals [Integers]
## numind = 35

## 7. Number of loci [Integers]
## numloci = 3560

## 8. What is the ploidy [Integers 1 through n]
## ploidy = 2

## 9. How is the missing data coded? Write inside quotes. e.g. missing = "-9"
## missing = "-9"

## 10. Does the data file contain every individual on 2 lines (0) or 1 line (1). [Boolean]
## onerowperind = 0 

## 11. Do the individuals have labels in the data file?  [Boolean]
## label = False

## 12. Are populations identified in the data file? [Boolean]
## popdata =  True

## 13. Do you wish to set the popflag parameter? [Boolean]
## popflag = False

## 14. Does the data file contain location identifiers (Not the same as population identifier) [Boolean]
## locdata = False

## 15. Does the data file contain phenotypic information? [Boolean]
## pheno = False

## 16. Does the data file contain any extra columns before the genotype data begins? [Boolean]
## extracols = False

## 17. Does the data file contain a row of marker names at the top? [Boolean]
## markers = True

## 18. Are you using dominant markers such as AFLP? [Boolean]
## dominant = False

## 19. Does the data file contain information on map locations for individual markers? [Boolean]
## mapdist = False

## 20. Is the data in correct phase? [Boolean]
## phase = True

## 21. Is the phase information provided in the data? [Boolean]
## phaseinfo = False

## 22. Does the phase follow markov chain? [Boolean]
## markov = False

## 23. Do you want to make use of parallel processing [Boolean]
##     Note that you must have GNU parallel installed and available
##     www.gnu.org/s/parallel
##     You can check availability of the program by running 'which parallel' at the 
##     command prompt. If a destination of the file is returned, then it is available.
##     If not available, it must be installed locally or through your system administrator.

## parallel = True

## 24. How would you like to define the number of cores for parallel processing ['number' or 'percent']
##     Use 'percent' if you would like to define the percentage of available cores to use.  For instance,
##     on a quad-core machine you might use 'percent' here and '75' for cores to use 3 of the 4 processors.
##     Use 'number' if you want to explicitely define the number of cores used.  

## core_def = 'number'

## 25. How many cores would you like to use [integer]
##     This represents either a pecentage or the physical number of cores as defined in core_def (#24).

## cores = 20

## 26. Would you like to automatically run structureHarvester?  [boolean]
##     Note that you must have program installed and available.
##     https://users.soe.ucsc.edu/~dearl/software/structureHarvester/

## harvest = True

##  End of questions 
##  Please do not write any other information in this file ###########


nano mainparams
Basic program parameters
#define MAXPOPS          8
#define BURNIN           250000
#define NUMREPS          500000

Input file
#define INFILE   SNP_3560_STRUCTURE.str

Data file format
#define NUMINDS          35
#define NUMLOCI          3560
#define PLOIDY           2
#define MISSING          -9
#define ONEROWPERIND     0

#define LABEL            0
#define POPDATA          1
#define POPFLAG          0
#define LOCDATA          0
#define PHENOTYPE        0
#define EXTRACOLS        0
#define MARKERNAMES	 1
#define RECESSIVEALLELES         0
#define MAPDISTANCES     0

Advanced data file options
#define PHASED           1
#define MARKOVPHASE	 0
#define NOTAMBIGUOUS     -999


## Once you have the results from Automation and Parallelization of STRUCTURE Analysis (StrAuto), go to http://taylor0.biology.ucla.edu/structureHarvester/ and upload the file ./harvester/SNP_3833_STRUCTURE_Harvester_Upload.zip to this server. Download the files results from the server and (remane the file) upload them to a directory called Harvester_results working directory and unzip it.


```bash
tar xvzf Harvester_results.gz 
```



## Most code adapted from: https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html

## PCA 

## A principal components analysis (PCA) converts the observed SNP data into a set of values of linearly
## uncorrelated variables called principal components that summarize the variation between samples.
## We can perform a PCA on our genlight object by using the glPCA function.
## To view the results of the PCA we can use the package ggplot2. We need to convert the data frame
## that contains the principal components (Orbicella_genlight_pca$scores) into the new object Orbicella_genlight_pca_scores. In addition, we will add the population values as a new column in our Orbicella_genlight_pca_scores object, in order to be able to color samples by population.

## ggplot2 will plot the PCA, color the samples by population, and create ellipses that include 95% of the data for each the population:



```r
library(ggplot2)
library(RColorBrewer)
library(colorspace)
#Install below package if necessary
#install.packages("poppr")
library(poppr)
library(colorspace)
library(vcfR)

VCF_3560_SNPs<-read.vcfR("TotalRawSNPs_BiAllelic_mac3_miss90_maf0.026_dp20_NoIndels_LD1000_No_clones.recode.vcf")
VCF_3560_SNPs

genlight_3560_SNPs <- vcfR2genlight(VCF_3560_SNPs)
genlight_3560_SNPs

genlight_3560_SNPs_pca <- glPca(genlight_3560_SNPs, nf=3)

genlight_3560_SNPs_pca_scores <- as.data.frame(genlight_3560_SNPs_pca$scores)
genlight_3560_SNPs_pca_scores

set.seed(9)
p <- ggplot(genlight_3560_SNPs_pca_scores, aes(x=PC1, y=PC2))
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major =
               element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
p <- p +theme
p
```



```r
barplot(100*genlight_3560_SNPs_pca$eig/sum(genlight_3560_SNPs_pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
head (genlight_3560_SNPs_pca$eig)
```




## DAPC
## We can further explore population assignments using a discriminant analysis of principal components (DAPC).
## Since we do not have any previous information about genetic structure, and all colonies come from the same geographical population, we need to find natural clusters in our data. We first use find.cluster (cluster identification using successive K-means): These functions implement the clustering procedure used in Discriminant Analysis of Principal Components (DAPC, Jombart et al. 2010). This procedure consists in running successive K-means with an increasing number of clusters (k), after transforming data using a principal component anal- ysis (PCA). For each model, a statistical measure of goodness of fit (by default, BIC) is computed, which allows to choose the optimal k.


```r
clusters_genlight_3560_SNPs <- find.clusters(genlight_3560_SNPs, max.n.clust = 8, n.pca =35, n.clust = 3  )
clusters_genlight_3560_SNPs
clusters_genlight_3560_SNPs$grp
#Choose the number PCs to retain (>=1): 
#35
#Choose the number of clusters (>=2: 
#3
```

## We choose to retain 13 PC that explained 80 % of the observed variance. If too few PCs (with respect to the number of individuals) are retained, useful information will be excluded from the analysis, and the resultant model will not be informative enough to accurately discriminate between groups. By contrast, if too many PCs are retained, this will have a destabilising effect on the coefficients extimated, leading to problems of overfit (adegenet package)

## Based on the lowest BIC value we choose  clusters (k=3)


```r
pop(genlight_3560_SNPs)<- clusters_genlight_3560_SNPs$grp
dapc_3560_SNPs <- dapc(genlight_3560_SNPs, n.pca = 35, n.da = 2)
dapc_3560_SNPs
scatter.dapc(dapc_3560_SNPs, scree.pca = TRUE)
```

## This DAPC scatter plot shows 3 groups because we overstimated by chossing all the PC. This is way we know want to evaluate the optimal number of PC to retain through the alfa score

## DAPC requires enough PCs to secure a space with sufficient power of discrimination but must also avoid retaining too many dimensions that lead to over-fitting.
## Using the a-score to measure the trade-off between power of discrimination and over-fitting. This score is simply the difference between the proportion of successful reassignment of the analysis (observed discrimination) and values obtained using random groups (random discrimination). It can be seen as the proportion of successful reassignment corrected for the number of retained PCs. It is implemented by a.score, which relies on repeating the DAPC analysis using randomized groups, and computing a-scores for each group, as well as the average a-score: 


```r
a_score <- optim.a.score(dapc_3560_SNPs)  
a_score
```



## Now we run dapc again but knowing that 1 PC is the optimal number and we keep the cluster assignation (k=3) previouly found with find.clusters. And wee keep 2 discriminant functions


```r
library(RColorBrewer)
pop(genlight_3560_SNPs) <- clusters_genlight_3560_SNPs$grp
colors <- brewer.pal(n = nPop(genlight_3560_SNPs), name = "Dark2")
dapc_1_PC <- dapc(genlight_3560_SNPs,clusters_genlight_3560_SNPs$grp, n.pca = 1, n.da = 2)
dapc_1_PC
dapc_1_PC$var
dapc_1_PC$grp

scatter(dapc_1_PC, col = colors, scree.pca = TRUE, scree.da = TRUE, posi.da = "topleft", posi.pca = "topright")
```

## Improving the DAPC chart: make the fond bigger and define the color palette 
##https://rdrr.io/cran/adegenet/man/dapcGraphics.html

```r
#install.packages("palettetown")
library(palettetown)
```
 

```r
#colors <- (values=c("#E84838","#5088B8","#78b850","#385860"))
colors <- (values=c("#78b850", "#E84838","#5088B8"))
scatter(dapc_1_PC,scree.da = FALSE, col =colors, cex=4, cstar=0, solid=.4,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3), posi.leg = "topright" )
clusters_genlight_3560_SNPs$grp
par(xpd=TRUE)

scatter(dapc_1_PC,scree.da = TRUE, posi.da = "bottomleft", col =colors, cex=4, cstar=0, solid=.4,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3), posi.leg = "topright" )
clusters_genlight_3560_SNPs$grp
par(xpd=TRUE)

scatter(dapc_1_PC,scree.da = FALSE,scree.pca = TRUE, posi.pca= "bottomleft", col =colors, cex=4, cstar=0, solid=.4,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3), posi.leg = "topright" )
clusters_genlight_3560_SNPs$grp
par(xpd=TRUE)

scatter(dapc_1_PC,scree.da = FALSE,scree.pca = TRUE, mstree = TRUE, posi.pca= "topright", col =colors, cex=4, cstar=0, solid=.4,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3), posi.leg = "topleft" )
clusters_genlight_3560_SNPs$grp
par(xpd=TRUE)

scatter(dapc_1_PC,scree.da = FALSE,scree.pca = TRUE, mstree = FALSE, posi.pca= "bottomleft",bg="grey" , grid = TRUE, col =colors,addaxes = TRUE, cex=4, cstar=0.5, solid=.4,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:3), posi.leg = "topright" )
```


## Final figure


```r
dapc_2_PC_2da <- dapc(genlight_3560_SNPs,clusters_genlight_3560_SNPs$grp, n.pca = 2, n.da = 2)
dapc_2_PC_2da
	
scatter (dapc_2_PC_2da, col = colors, posi.da = "bottomright", cex=4, leg=TRUE, txt.leg=paste("Cluster",1:3), posi.leg = "topright")
```


## https://www.rdocumentation.org/packages/adegenet/versions/2.0.1/topics/dapc%20graphics


```r
colors <- (values=c("#78b850","#E84838", "#5088B8"))

scatter (dapc_2_PC_2da, col = colors, scree.da = F, posi.da = "bottomright", cex=4, bg="#F5F5F5", leg=F, txt.leg=paste("Cluster",1:3), posi.leg = "bottomleft", cellipse = F, clabel=F, addaxes = T, grid=F,mstree = F, cstar = F, cgrid=T)
```



```r
dapc_1_PC_1da <- dapc(genlight_3560_SNPs,clusters_genlight_3560_SNPs$grp, n.pca = 1, n.da = 2)
dapc_1_PC_1da

scatter (dapc_1_PC_1da, col = colors, scree.da = F, posi.da = "bottomright", cex=2, bg="white", leg=F, txt.leg=paste("Cluster",1:3), posi.leg = "bottomleft", cellipse = F, clabel=F, addaxes = F, grid=F,mstree = F, cstar = F, cgrid=T, onedim.filled=T)
```


# Attempt to grapgh kernel density pot for LD2

# http://lists.r-forge.r-project.org/pipermail/adegenet-forum/2014-December/001033.html


```r
scatter(dapc_2_PC_2da,xax = 2, scree.da = F, col =colors, cex=2,, cstar=0,clab=0, grid=T, bg="white", addaxes=T, cellipse = F, mstree = F, posi.da = "topright")
clusters_genlight_3560_SNPs$grp
par(xpd=TRUE)
```






## Distance Tree

```r
new_names = values= c("10","11","12","14","16","17","18","19","1","22","23","24","25","26","28","30","31","32","33","34","35","36","37","38","39","3","40","41","44","45","5","6","7","8","9")
indNames(genlight_3560_SNPs) <- new_names
indNames(genlight_3560_SNPs)

Orbicella_tree <- aboot(genlight_3560_SNPs, tree = "upgma", distance = bitwise.dist, sample = 1000, showtree = T, cutoff = 50, quiet = T)
```



## Cutoff = 0

```r
new_names = values= c("10","11","12","14","16","17","18","19","1","22","23","24","25","26","28","30","31","32","33","34","35","36","37","38","39","3","40","41","44","45","5","6","7","8","9")
indNames(genlight_3560_SNPs) <- new_names
indNames(genlight_3560_SNPs)

Orbicella_tree <- aboot(genlight_3560_SNPs, tree = "upgma", distance = bitwise.dist, sample = 1000, showtree = T, cutoff = 0, quiet = T )
```




```r
library(poppr)
library(ape)


pop(genlight_3560_SNPs)
pop(genlight_3560_SNPs)<- dapc_1_PC$grp
# colors <- (values=c("#E84838","#5088B8","#78b850","#385860"))
colors <- (values=c("#E84838","#5088B8","#78b850"))
nPop(genlight_3560_SNPs)


tree<-plot.phylo(Orbicella_tree, cex = 0.8, font = 2, adj = 0, tip.color =  colors[pop(genlight_3560_SNPs)], use.edge.length = T)
tree
nodelabels(Orbicella_tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend("left", legend = c("PAN_1","PAN_2","PAN_3"), fill = colors, border = FALSE, bty = "o", cex = 0.8, inset = 0.1, 0)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")

#write.tree(Orbicella_tree, file = "Ofav_NJ_Tree_v2.txt", append = T, digits = 10)
```




```r
plot.phylo(Orbicella_tree, cex = 0.8, font = 2, adj = 0, tip.color =  colors[pop(genlight_3560_SNPs)], type = "phylo")
nodelabels(Orbicella_tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend("left", legend = c("PAN_1","PAN_2","PAN_3"), fill = colors, border = FALSE, bty = "o", cex = 1, inset = 0.1, 0)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different")
```



