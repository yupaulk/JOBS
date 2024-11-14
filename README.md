# JOBS
JOint model viewing Bk-eQTLs as a weighted sum of Sc-eQTLs (JOBS)


## Table of contents
* [Introduction](#Introduction)
* [Installation](#Installation)
* [Quick tutorial](#Quick_tutorial)
* [Input files](#Input_files)
* [Contact](#Contact)

## Introduction
Here, we propose a JOint model viewing Bk-eQTLs as a weighted sum of Sc-eQTLs (JOBS) from constituent cell types. JOBS borrows strength from large sample sizes of bk-eQTLs to improve sc-eQTLs analysis, with improvements bigger for more common cell types. With more accurate sc-eQTLs effect estimates, all integrative analyses can be improved.
 
It is developed and maintained by Lida Wang at [Dajiang Liu's Group](https://dajiangliu.blog).

## Installation
The package is hosted on github, which allows installation and update to be very easy. First, make sure you have the dplyr, data.table, tidyr, R.utils, pracma, and quadprog packages installed.

```
install.packages("devtools")
library(devtools)
```

Then you could install JOBS from the repository here.

```
devtools::install_github("LidaWangPSU/JOBS/JOBS")
library(JOBS)
```


## Quick tutorial
### Prepare bulk and single cell eQTLs summary statistics (effect size and s.e.). 

We have example input data [here](https://github.com/LidaWangPSU/JOBS/blob/main/example_data/). Data were subsetted from brain bulk and single cell eQTLs as an example to run the script.

### Input Data Format

#### 1. Effect Size Matrix
This matrix contains eQTL effect sizes for both bulk and single-cell data.

- **Columns**:
  - **Column 1**: Gene-SNP pair identifiers.
  - **Column 2**: Bulk effect size.
  - **Columns 3+**: Cell type-specific eQTL effect sizes.
**Example**:
|    Gene-snp pair    |      Bulk     | Cell type 1  |  Cell type 2 | ...... |  Cell type k |
| ------------------- |      ----     | -----------  |  ----------- | ------ |  ----------- |
| ENSG00000XXXXX-snp1 |      0.65     |      0.55    |      0.30    | ...... |      0.75    |  
| ENSG00000XXXXX-snp2 |     -0.50     |     -0.62    |     -0.35    | ...... |      0.02    |


#### 2. Standard Error Matrix (S.E.)
This matrix has the same dimensions as the effect size matrix and represents the standard errors associated with the eQTL effect sizes.

- **Columns**:
  - **Column 1**: Gene-SNP pair identifiers (should match the effect size matrix).
  - **Column 2**: Bulk effect size standard error.
  - **Columns 3+**: Standard errors for the cell type-specific eQTL effect sizes.
**Example**:
|    Gene-snp pair    |      Bulk     | Cell type 1  |  Cell type 2 | ...... |  Cell type k |
| ------------------- |      ----     | -----------  |  ----------- | ------ |  ----------- |
| ENSG00000XXXXX-snp1 |    0.0071     |    0.0316    |    0.0316    | ...... |    0.0316    |  
| ENSG00000XXXXX-snp2 |    0.0071     |    0.0316    |    0.0316    | ...... |    0.0316    |

(It is acceptable to have missing data in the input. We will analyze all gene-SNP pairs unless sc-eQTL data is entirely missing for all cell types. In other cases, we will simply ignore the missing values and keep "NA" in the output.)

```
library(data.table)
beta <- as.data.frame(fread("~/example_beta_chr22.txt.gz"))
se <- as.data.frame(fread("~/example_beta_se_chr22.txt.gz"))
```
  
### Run JOBS
#### JOBS Step1: cell type weights
Here, we used Non-negative least squares to estimate cell type weights.
You can also specify the weights by yourself, e.g. estimated from scRNAseq data or other methods
```
weight <-jobs.nnls.weights(beta,se)
```


#### JOBS Step2: refine eQTLs edtimation 
```
jobs_eqtls <- jobs.eqtls(beta,se,weight)
```

### Output results
The output includes:
* Refined effect size: in the same format as effect size input
* Refined S.E: in the same format as s.e. input


## Contact
Lida Wang [lida.wang.96@gmail.com](lida.wang.96@gmail.com)
