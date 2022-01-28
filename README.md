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

## Stage 1 - select simulation parameters
At this stage, one can use our interactive tool, Phemulator to select the most appropriate/interesting paramaters for the simulation stage.
Phemulator can be accessed via a web browser and is available at [https://wwlc3x-marcin0kierczak.shinyapps.io/pheno_sim2/](https://wwlc3x-marcin0kierczak.shinyapps.io/pheno_sim2/). 

## Stage 2 - prepare genomic data for the simulation
At this stage, the genomic data are filtered, pre-processed and prepared for being used at the later stages of the simulation.

## Stage 3 - simulate phenotypes
This stage is the actual simulation where the data prepared at stage 2, together with parameters selected at stage 1 are used to simulate phenotypes.

## Stage 4 - perform SVA and gene-based test using SAIGE 
This is probably the most exciting stage, where the phenotypes simulated at stage 3 are used to test the behavior of standard single-variant association as well as gene-based tests.

## Possible applications:
1. To validate sensitivity and specificity of various GWA algorithms in the landscape of varying effect sizes, directionalities and mafs. 
  
![](assets/Simulations_diagram.png?raw=true)
