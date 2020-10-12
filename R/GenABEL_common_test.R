library('GenABEL')
source('R/get_effects.R')
source('R/plot_effects.R')
source('R/impute_G.R')
source('R/get_genotypes.R')
source('R/fix_allele_encoding.R')
source('R/tibble_to_raw_genotypes.R')
source('R/tibble_to_gwaa.R')

plot_pheno <- function(y, G) {
  data <- tibble(id = rownames(G), y = y)
  ggplot(data, mapping = aes(x=y)) + geom_histogram() + theme_bw()
}

data(srdta)
nids(srdta)
nsnps(srdta)

# Extract marker summary data, compute reference allele frequency from allele counts and compute minor allele frequency.
tmp <- summary(srdta@gtdata)
raf <- ((2 * tmp$P.22) + tmp$P.12) / (2 * tmp$NoMeasured)
maf <- pmin(raf, 1 - raf)
names(maf) <- colnames(srdta@gtdata)

# Get effects for common alleles
common <- get_effects(maf = maf, thr = 0.05,
                              N = 5,
                              shape12 = c(.1,.1),
                              rare = F,
                              frac_negative = 0)

rare <- get_effects(maf = maf, thr = 0.01,
                      N = 20,
                      shape12 = c(1,25),
                      rare = T,
                      frac_negative = .2)

# Get genotypes, convert them to double, make sure the genotypes matrix contains the minor allele count, impute missing
genos_common <- srdta@gtdata[,common$marker_idx]
G_comm <- as.double(genos_common) %>% fix_allele_encoding() %>% impute_G(maf = maf[common$marker_idx])
genos_rare <- srdta@gtdata[,rare$marker_idx]
G_rare <- as.double(genos_rare) %>% fix_allele_encoding() %>% impute_G(maf = maf[rare$marker_idx])

# Compute phenotype values and plot them
y_comm <- G_comm %*% common$effects + rnorm(n = dim(G_comm)[1], mean = 0, sd = 1)
y_rare <- G_rare %*% rare$effects + rnorm(n = dim(G_rare)[1], mean = 0, sd = 1)

plot_pheno(y_comm, G_comm)
plot_pheno(y_rare, G_rare)

# Inject the simulated phenotype to the gwaa.data object
colnames(srdta@gtdata)[rare$marker_idx]
srdta@phdata[,'qt1'] <- y_rare
(simulated <- names(rare$effects))
#snpsubset = !(snpnames(srdta) %in% c('rs2935'))

# Do GWAS
top <- G[,2]
qt <- qtscore(qt1~sex, data = srdta)
plot(qt, las=1, cex.axis = .7, col = 'darkgreen', cex=.5)
grid()
col <- as.numeric(rare$effects > 0)
col[col == 1] <- 'red'
col[col == 0] <- 'blue'
#abline(v=qt[simulated,]$Position, col=col)
points(qt[simulated,]$Position, -log10(qt[simulated,]$P1df), col=col, pch=19)

# Predicted effect vs. simulated effect
plot(rare$effects, qt[simulated, ]$effB, cex = .5, pch = 19, las=1, cex.axis=.7, xlab='simulated effect', ylab='predicted effect')
grid(col='grey')
abline(a = 0, b=1, col='lightgrey')
