# JOBS
JOint model viewing Bk-eQTLs as a weighted sum of Sc-eQTLs (JOBS)![image](https://github.com/LidaWangPSU/JOBS/assets/70541830/df17d42c-4a3a-49f7-a418-ee0718be9c63)


## Table of contents
* [Introduction](#Introduction)
* [Installation](#Installation)
* [Quick tutorial](#Quick_tutorial)
* [Input files](#Input_files)
* [Contact](#Contact)

## Introduction
Most variants identified from genome-wide association studies (GWAS) are non-coding and function by regulating gene expression. Yet, many GWAS loci fail to colocalize with expression quantitative trait loci, possibly due to limited power for GWAS and eQTL analysis and cellular heterogeneity of effect sizes. Population-scale single-cell (sc) RNASeq datasets with hundreds of individuals are beginning to emerge where we can group cells into pseudo-bulks to study eQTLs for different cell types (sc-eQTLs). Yet, compared to eQTL data from bulk tissues (bk-eQTLs), sc-eQTL datasets are much smaller. Here, we propose a JOint model viewing Bk-eQTLs as a weighted sum of Sc-eQTLs (JOBS) from constituent cell types. JOBS borrows strength from large sample sizes of bk-eQTLs to improve sc-eQTLs analysis, with improvements bigger for more common cell types. With more accurate sc-eQTLs effect estimates, all integrative analyses can be improved.
 
It is developed and maintained by Lida Wang at [Dajiang Liu's Group](https://dajiangliu.blog).

## Installation
The package is hosted on github, which allows installation and update to be very easy. First, make sure you have the dplyr, data.table, and tidyr packages installed.

```
install.packages("devtools")
library(devtools)
```

Then you could install JOBS from the repository here.

```
devtools::install_github("LidaWangPSU/JOBS/JOBS")
library(JOBS)
```
Here we go.

## Quick tutorial
### Prepare bulk and single cell eQTLs summary statistics (effect size and s.e.). We have example input data [here](https://github.com/LidaWangPSU/JOBS/blob/main/example_data/). 
Data were subsetted from brain bulk and single cell eQTLs as an example to run the script.

Input includes
* Effect size: A matrix of eqtls across bulk and single cell, first col: gene-snp pair; second col: bulk effect size; 3+ cols: cell type specific eqtls
* S.E.: A matrix of eqtls standard deviation across bulk and single cell, first col: gene-snp pair; second col: bulk effect size se; 3+ col: cell type specific eqtls se. S.E. file should be 1-1 match with effect size
```
library(data.table)
beta <- as.data.frame(fread("~/example_beta_chr22.txt.gz"))
beta <- as.data.frame(fread("~/example_beta_se_chr22.txt.gz"))
```
  
### Run JOBS
```
jobs_eqtls <- jobs(beta,se)
```

### Output results
We have three outputs

The weight output includes:
* Cell type weights: weights for each cell types
* Refined effect size: in the same format as effect size input
* Refined S.E:in the same format as s.e. input


## Contact
Lida Wang [lida.wang.96@gmail.com](lida.wang.96@gmail.com)
