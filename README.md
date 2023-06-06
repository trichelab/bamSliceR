Rapid expressed variant and allelic bias detection for rare variants and
rare diseases
================
5/10/2022

## Contact Infomation

## Installing

<!-- If you're putting `mutationScreenR` on CRAN, it can be installed with

    install.packages("mutationScreenR") -->

The pre-release version of the package can be pulled from GitLab using
the [devtools](https://github.com/hadley/devtools) package:

    # install.packages("devtools")
    url <- ""
    devtools::install_git(url = url)

The package also could be loaded and further developed using the
build-in git functionality in Rstudio:
<https://pjnewcombe.wordpress.com/2014/05/31/using-rstudio-to-create-and-maintain-an-r-package-on-github/>

## Data processing

## Example on how to create the target genes template

#### Only keep the core histone genes

#### Get all the granges for all the exons region of all genes

#### Keep the exons region only for the core histone genes we interested in

## read test files of 9 patients that have samples from multi-time points: Normal vs. Diagnosis vs. Relapse

#### Check the results. The results could be used for downstream analysis including Fish exact test and VEP

## An exmple of 327 acute-myeloid Leukemia (AML) patients

The example usage of this software is focused on genetic mutations
screening among the 327 acute-myeloid Leukemia (AML) patients on the
histones proteins. Most of these 327 patients has samples from normal
tissue, blood samples(contains cancer cells) from both diagnosis and
relapse. The results could be viewed:

## Non-negative Matrix (NMF) for unsupervised clustering of AML patients based on expression data.

This pacakge integrate the Non-negative matrix factorization algorithm
to unsupervised cluster the patients (with same disease or similar
pathologic condition) based on the expression data from RNA-sequencing
technique.

#### Load Expression data

#### Based on the rank survey, we can select the k for NMF algorithm

#### run NMF with k = 8, time consuming step

#### Plotting the results of NMF
