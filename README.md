# a_johansson_2020
Åsa Johansson partner project 2020

<!-- badges: start -->
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Codecov test coverage](https://codecov.io/gh/NBISweden/a_johansson_2020/branch/master/graph/badge.svg)](https://codecov.io/gh/NBISweden/a_johansson_2020?branch=master)
[![R build status](https://github.com/NBISweden/a_johansson_2020/workflows/R-CMD-check/badge.svg)](https://github.com/NBISweden/a_johansson_2020/actions)
<!-- badges: end -->

## Objectives
To develop an R package for testing the scope of applicability of different GWA methodologies, esp. with respect ot varying:
* threshold between rare and common variants,
* the degree of contribution (effect size and direction) contributed by rare and common variants,
* degree of population structure,
* varying amount of contributing loci.

## Using
To use the package on Bianca:
* the `gwasim` package is automatically built into a docker container upon every push,
* container is called gwasim-latest and is stored in `quiestrho` account on DockerHub,
* ssh to Rackham, do `singularity pull --docker-login docker://quiestrho/gwasim-latest`
* transfer the `gwasim_latest.sif` file into Bianca's wharf via sftp,
* move the file from wharf to your project library,
* `docker run gwasim-latest.sif < script.r` to run an R script within the container 

## Input
The following input parameters are expected from the user:
* a VCF file with variants coming from a population under studies,
* a bed file containing a list of functional regions (e.g. gene list) along with their coordinates,
* a minor allele frequency threshold for cut-off between rare and common variants,
* the number of rare and common variants that contribute,
* distribution of effect sizes as a function of allele frequency (separately for the rare and the common alleles),
* percentage of alleles with negative effect (for both the common and the rare variants),
* parameters (mean, standard deviation) of the error term,
* type of the trait (continuous or binary, a cut-off value for the binary trait),

## To add in future:
* kinship matrix for the studied population,
* fixed effects,

## Output
* A text file or a tibble with simulated phenotypes along with the effect of the loci used for simulation.

## Possible applications:
1. To validate sensitivity and specificity of various GWA algorithms in the landscape of varying effect sizes, directionalities and mafs. 

![](assets/Simulations_diagram.png?raw=true)
