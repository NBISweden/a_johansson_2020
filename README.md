# a_johansson_2020
AJ partner project 2020-2022

<!-- badges: start -->
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Codecov test coverage](https://codecov.io/gh/NBISweden/a_johansson_2020/branch/master/graph/badge.svg)](https://codecov.io/gh/NBISweden/a_johansson_2020?branch=master)
[![R build status](https://github.com/NBISweden/a_johansson_2020/workflows/R-CMD-check/badge.svg)](https://github.com/NBISweden/a_johansson_2020/actions)
![build-and-push-to-DockerHub](https://github.com/NBISweden/a_johansson_2020/workflows/build-and-push-to-DockerHub/badge.svg)
<!-- badges: end -->

## Objectives
To develop a set of simulation tools and a pipeline to, given whole-genome or genotyping data for a large cohort, simulate phenotypes according to a genetic architecture with specified parameters:
* threshold between rare and common variants,
* the degree of contribution (effect size and direction) contributed by rare and by the common variants,
* degree of population structure,
* varying number of contributing loci with common and rare variants.

## Requirements
* snakemake (>= 6.9.1),
* singularity (>= 3.8.5-2.el7),
* R (>= 3.6.0)

## Input
The following input is expected:
* parameters that describe genetic architecture one wants to simulate (population size, common/rare maf threshold, distribution parameters for effect sizes, number of contributing loci),
* a bed file that specifies genomic regions to be used in simulations, e.g., CDS regions,
* a file containing a list of individuals to be kept in the dataset (in PLINK2 format),
* a trio of files (bed, bim, fam) that contains genotyping/sequencing data in PLINK format.

## Stage 1 - Phemulator: select simulation parameters using interactive tool
At this stage, one can use our interactive tool, Phemulator to select the most appropriate/interesting paramaters for the simulation stage.
Phemulator can be accessed via a web browser and is available at [https://wwlc3x-marcin0kierczak.shinyapps.io/pheno_sim2/](https://wwlc3x-marcin0kierczak.shinyapps.io/pheno_sim2/). 

![](assets/phemulator.png?raw=true)

## Stage 2 - prepare genomic data for the simulation
At this stage, the genomic data are filtered, pre-processed and prepared for being used at the later stages of the simulation.

![](assets/a_johanssonPP_stage2.drawio.png?raw=true)

## Stage 3 - simulate phenotypes
This stage is the actual simulation where the data prepared at stage 2, together with parameters selected at stage 1 are used to simulate phenotypes.
```
Usage: R/optparse_test.R [options]


Options:
        --chr-name=CHR-NAME
                name of the chromosome top process, e.g. chr22

        --data-path=DATA-PATH
                path to data files

        --regions-bed=REGIONS-BED
                bed file defining functional regions of interest, eg. CDS

        --afreq=AFREQ
                an afreq file with allele frequencies computed by Plink2

        --fam=FAM
                fam file

        --plink-prefix=PLINK-PREFIX
                prefix (no extension) of the plink2 data files

        --min-maf=MIN-MAF
                minimal maf to consider (used to remove ultra-rare/fixed markers)

        --rare-maf=RARE-MAF
                markers above this threshold will be considered common

        --num-mrk=NUM-MRK
                number of markers per region to be used for simulating phenotype

        --num-mrk-neg=NUM-MRK-NEG
                per region number of markers with negative effect

        --mean-eff=MEAN-EFF
                mean effect of a rare allele

        --sd-eff=SD-EFF
                std. dev. for the distribution of rare allele effects

        --mean-err=MEAN-ERR
                mean error (residuals)

        --sd-err=SD-ERR
                std. dev. for the distribution errors (residuals)

        --num-sim=NUM-SIM
                number of simulations to run (regions to use)

        -h, --help
                Show this help message and exit
```

![](assets/a_johanssonPP_stage3.drawio.png?raw=true)

## Stage 4 - perform SVA and gene-based test using SAIGE 
This is probably the most exciting stage, where the phenotypes simulated at stage 3 are used to test the behavior of standard single-variant association as well as gene-based tests.

## Stage 5 - visualize the results
```
R -e "rmarkdown::render('script.Rmd',output_file='output.html')"
```


## Possible applications:
1. To validate sensitivity and specificity of various GWA algorithms in the landscape of varying effect sizes, directionalities and mafs. 

