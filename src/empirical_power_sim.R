library(renv)
library(vcfR)
library(tidyverse)
library(ggplot2)
data(vcfR_example)
source('R/get_effects.R')
source('R/plot_effects.R')
source('R/impute_G.R')

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

G <- as_tibble(vcfR::vcfR2genind(vcf)@tab[,rare$marker_idx])
G <- G %>% mutate_all(~ stringr::str_count(string = ., "1"))
G <- t(as.matrix(G))
# Impute
G_imp <- impute_G(G = G, maf = maf)
y <- G_imp %*% rare$effects + rnorm(n = dim(G_imp)[1], mean = 0, sd = 1)
plot(y)


