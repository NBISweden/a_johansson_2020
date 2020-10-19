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
source('R/tibble_to_raw_genotypes.R')
source('R/tibble_to_gwaa.R')

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
#tmp <- G[,1]
#tmp[tmp == 0] <- 3
#tmp[tmp == 2] <- 0
#tmp[tmp == 3] <- 2
#G[,1] <- tmp


y <- G %*% rare$effects + rnorm(n = dim(G)[1], mean = 0, sd = 1)
plot(y, 1:dim(G)[1], las = 1, cex.axis = .7, ylab="Individual", xlab="Phenotype")
abline(v=0, col="grey")
grid()

geno <- as_tibble(t(G)) %>% add_column(snp=colnames(G), .before = 1) %>%
    add_column(chr = c(1), .after = 1) %>%
    add_column(pos = c(1:dim(G)[2]), .after = 2) %>%
    add_column(strand = c(0), .after = 3)
data <- tibble_to_gwaa(x = geno, sex = rep(0, times=18), trait = y)
qt <- qtscore(y~1, data = data)
plot(qt)
