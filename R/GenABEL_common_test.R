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
names(maf) <- colnames(srdta@gtdata)
common <- get_effects(maf = maf, thr = 0.05,
                              N = 1,
                              shape12 = c(.1,.1),
                              below = F,
                              perc_negative = 0)
genos <- srdta@gtdata[,common$marker_idx]
G <- as.double(genos) %>% fix_allele_encoding() %>% impute_G(maf = maf[common$marker_idx])
y <- G %*% common$effects + rnorm(n = dim(G)[1], mean = 0, sd = 1)
plot(y, 1:dim(G)[1], las = 1, cex.axis = .7, ylab="Individual", xlab="Phenotype")
abline(v=0, col="grey")
grid()

colnames(srdta@gtdata)[common$marker_idx]
srdta@phdata[,'qt1'] <- y
(simulated <- names(common$effects))
#snpsubset = !(snpnames(srdta) %in% c('rs2935'))
qt <- qtscore(qt1~sex, data = srdta)
plot(qt)
abline(v=qt[simulated,]$Position, col='red')

