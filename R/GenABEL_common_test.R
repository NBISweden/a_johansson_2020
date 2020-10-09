library('GenABEL')
source('R/get_effects.R')
source('R/plot_effects.R')
source('R/impute_G.R')
source('R/get_genotypes.R')
source('R/fix_allele_encoding.R')
source('R/tibble_to_raw_genotypes.R')
source('R/tibble_to_gwaa.R')

data(srdta)
nids(srdta)
nsnps(srdta)

tmp <- summary(srdta@gtdata)
raf <- ((2 * tmp$P.11) + tmp$P.12) / (2 * tmp$NoMeasured)
maf <- pmin(raf, 1 - raf)
common <- get_effects(maf = maf, thr = 0.05,
                              N = 1,
                              shape12 = c(.5, .5),
                              below = F,
                              perc_negative = 0)
genos <- srdta@gtdata[,common$marker_idx]
as.double(genos) %>% fix_allele_encoding() %>% impute_G(maf = maf[common$marker_idx])

