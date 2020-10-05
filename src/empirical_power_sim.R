library(renv)
library(vcfR)
library(tidyverse)
library(ggplot2)
data(vcfR_example)
source('R/get_effects.R')
source('R/plot_effects.R')
source('R/impute_G.R')
source('R/get_genotypes.R')
source('R/fix_allele_encoding.R')

perc_negative_common <- 0
perc_negative_rare <- 0.2
n_common <- 1
n_rare <- 10
thr_common_rare <- 0.1
beta_params_common <- c(.5, .5)
beta_params_rare <- c(1, 25)

maf <- vcfR::maf(vcf)[,'Frequency']
rare <- get_effects(maf = maf, thr = thr_common_rare,
                      N = n_rare,
                      shape12 = beta_params_rare,
                      below = T,
                      perc_negative = perc_negative_rare)

G <- get_genotypes(x = vcf, marker_names = names(rare$marker_idx)) %>%
  fix_allele_encoding() %>%
  impute_G(maf = maf)

# Create a test
tmp <- G[,1]
tmp[tmp == 0] <- 3
tmp[tmp == 2] <- 0
tmp[tmp == 3] <- 2
G[,1] <- tmp


y <- G_imp %*% rare$effects + rnorm(n = dim(G_imp)[1], mean = 0, sd = 1)
plot(y)



