library('gwasim')

plot_pheno <- function(y, G) {
	data <- tibble(id = rownames(G), y = y)
	ggplot(data, mapping = aes(x=y)) + geom_histogram() + theme_bw()
}

data_chr1 <- load.gwaa.data(
	pheno = "/proj/sens2016007/nobackup/julia/GENABEL/generalDATA/pheno-pea3-genabel.dat",
	geno = "/proj/sens2016007/nobackup/marcin/biallelic-chr1.raw"
)

srdta <- data_chr1

# Extract marker summary data, compute reference allele frequency from allele counts and compute minor allele frequency.
tmp <- summary(srdta@gtdata)
raf <- ((2 * tmp$P.22) + tmp$P.12) / (2 * tmp$NoMeasured)
maf <- pmin(raf, 1 - raf)
names(maf) <- colnames(srdta@gtdata)

run_tests <- function(data, maf, times=3, N=1, thr=0.01, shape12=c(.1,.1), rare=F, frac_negative=0, e = c(0,1)) {
  result <- matrix(NA, nrow = times, ncol = N + 1)
  	i <- 1
    for (i in 1:times) {
    	print(paste0('Iteration ', i, ' in progress...'))
		effects <- get_effects(maf=maf, thr=thr, N=N, shape12=shape12, rare=rare, frac_negative=frac_negative)
		genotypes <- data@gtdata[ , effects$marker_idx]
		G <- as.double(genotypes) %>% fix_allele_encoding() %>% impute_G(maf = maf[effects$marker_idx])
		y <- G %*% effects$effects + rnorm(n = dim(G)[1], mean = e[1], sd = e[2])
		data@phdata[,'IL.8'] <- y
		genabel_model <- qtscore(IL.8 ~ sex + age, data = data, clambda = F)
		result[i,] <- c(genabel_model@results$P1df[effects$marker_idx], min(genabel_model@results$P1df))
		i <- i + 1
	}
	return(result)
}

size <- c(nids(srdta), nsnps(srdta))
result <- run_tests(times = 10000)
save(size, result, file='common_run_results.Rda')

#skat_null_model <- SKAT_Null_Model(y~sex, data = data@phdata, y=y)
#skat_model <- SKAT(G, skat_null_model, method = 'SKATO')


