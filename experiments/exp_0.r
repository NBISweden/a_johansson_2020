# This is the initial test. We want to make sure that gwasim works on Rckham.

library('stringr')
library('gwasim')
library('seqminer')
library('SKAT')
library('purrr')
library('tidyverse')
library("doFuture")
library("doRNG")

# Enable parallel processing
registerDoFuture()
plan("multisession")

# Data source http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
# User variables
bed_file <- "/home/marcin/hadrien_proj/nobackup/marcin/data/promoters_chr22.bed"
vcf_file <- "/home/marcin/hadrien_proj/nobackup/marcin/ALL.chr22.freebayes.20200518.snps_indels.NYhc.GRCh38.vcf.gz"
#vcf_file <- "/home/marcin/hadrien_proj/nobackup/marcin/chr22_sample.vcf.gz"
my_chr <- "chr22"
min_n_markers <- 5
n_simulations <- 100

regions <- read_regions_bed(bed_file) %>% filter(chr == my_chr)
regions$region <- paste0("chr", regions$region)
regions_maf_counts <- count_by_maf(regions, vcf_file = vcf_file, maf_threshold = 0.01)

focal_regions <- regions_maf_counts %>%
  filter(n_loci != 0) %>%
  filter(n_leq >= min_n_markers) %>%
  select(region) %>%
  as_vector() %>%
  sample(size = n_simulations, replace = T)
phenos <- sim_y_in_regions(focal_regions, vcf_file, n_markers = min_n_markers, e_sd = 0.001)

formula <- as.formula("y ~ 1")
results <- vector(mode = "list", length = length(phenos))
i <- 1
for (pheno in phenos) {
  focal_region <- pheno$region
  y <- pheno$y
  null_model <- SKAT::SKAT_Null_Model(formula)
  results[[i]] <- scan_CommonRare(focal_region = focal_region, regions = regions, null_model = null_model, vcf_file = vcf_file)
  i <- i + 1
}
save.image(file = 'tidying_testbed_lowest_e.Rdata')


ap_data <- unlist(lapply(results, count_ap))
hist(ap_data/nrow(regions))

midpoints <- regions$start + (regions$end - regions$start)/2
names(midpoints) <- regions$region
plot(x = midpoints, y = rep(0, times=length(midpoints)), pch=19, cex=.5, ylim = c(-.1,1.1), las=1, cex.axis=.7, bty='n',
     xlab='coords', ylab = 'fraction of associated positives')
grid()
points(x = midpoints, y = rep(0, times=length(midpoints)), pch=19, cex=.5)
names(ap_data) <- focal_regions
x <- midpoints[names(ap_data)]
points(x, ap_data/dim(regions)[1], type='h')
points(x, rep(0, length(x)), col='red', pch=19, cex=.5)

