library(tidyverse)
library(seqminer)
library("optparse")
source('R/gp_remap.R')
source('R/fix_encoding.R')
source('R/impute.R')
source('R/read_afreq.R')
source('R/read_regions.R')
source('R/scan_regions.R')

option_list = list(
  make_option(c("--chr-name"), type="character", default="chr22",
              help="name of the chromosome top process, e.g. chr22"),
  make_option(c("--data-path"), type="character", default="/home/marcin/ExomeSeq2ndRel/marcin/",
              help="path to data files"),
  make_option(c("--regions-bed"), type="character", default="GRCh38_chr22_cds.bed",
              help="bed file defining functional regions of interest, eg. CDS"),
  make_option(c("--afreq"), type="character", default="sntst_GRCh38_norel_rnd10000_chr22.afreq",
              help="an afreq file with allele frequencies computed by Plink2"),
  make_option(c("--fam"), type="character", default="sntst_GRCh38_norel_rnd10000_chr22.fam",
              help="fam file"),
#  make_option(c("--bim"), type="character", default="sntst_GRCh38_norel_rnd10000_chr22.bim",
#              help="bim file"),
  make_option(c("--plink-prefix"), type="character", default="sntst_GRCh38_norel_rnd10000_chr22",
              help="prefix (no extension) of the plink2 data files"),
  make_option(c("--min-maf"), type="numeric", default = 0.005,
              help = "minimal maf to consider (used to remove ultra-rare/fixed markers)"),
  make_option(c("--rare-maf"), type="numeric", default = 0.01,
              help = "markers above this threshold will be considered common"),
  make_option(c("--num-mrk"), type="numeric", default = 5,
              help = "number of markers per region to be used for simulating phenotype"),
  make_option(c("--num-mrk-neg"), type="numeric", default = 0,
              help = "per region number of markers with negative effect"),
  make_option(c("--mean-eff"), type="numeric", default = 1,
              help = "mean effect of a rare allele"),
  make_option(c("--sd-eff"), type="numeric", default = 0.01,
              help = "std. dev. for the distribution of rare allele effects"),
  make_option(c("--mean-err"), type="numeric", default = 0,
              help = "mean error (residuals)"),
  make_option(c("--sd-err"), type="numeric", default = 0.05,
              help = "std. dev. for the distribution errors (residuals)"),
  make_option(c("--num-sim"), type="numeric", default = 10,
              help = "number of simulations to run (regions to use)"),
  make_option(c("--output-path"), type="character", default="./",
	      help = "specification of the output path")
);

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

my_chr <- opt$`my-chr` # chromosome to run simulations on
path <- opt$`data-path`  # path to where the raw data files live
regions_file <- paste0(path, opt$`regions-bed`)
afreq_file <- paste0(path, opt$`afreq`)
fam_file <- paste0(path, opt$`fam`)
# bim_file <- paste0(path, opt$`bim`)
plink_file <- paste0(path, opt$`plink-prefix`)

threshold_low <- opt$`min-maf`   # maf has to be above this threshold
threshold_upp <- opt$`rare-maf`  # maf has to be equal to or below this threshold

# Simulation parameters
N <- opt$`num-mrk`              # number of markers to be used per simulation
N_neg <- opt$`num-mrk-neg`      # number of markers with negative effect size
eff_mu <- opt$`mean-eff`        # mean for effects distribution (effect sizes will be sampled from N(eff_mu, eff_sd))
eff_sd <- opt$`sd-eff`          # std. dev. for phenotype distribution
e_mu <- opt$`mean-err`          # mean for residuals distribution
e_sd <- opt$`sd-err`            # std. dev. for residuals distribution
N_sim <- opt$`num-sim`          # number of traits to simulate

# Read file with regions (colllapse to one CDS per region)
cds_uniq <- read_regions(file = regions_file)

# Read allele frequencies and remove all monomorphic markers 
# as well as all the markers that are
# above the maf upper threshold, i.e. the common markers
afreq_all <- read_afreq(file = afreq_file)
afreq <- afreq_all %>% filter(freq_alt > threshold_low & freq_alt <= threshold_upp) %>%
  filter(ref != "I" & ref != "D" & alt != "I" & alt != "D")

# Read the fam file
print(paste0("Reading family information from: ", fam_file))
fam_data <- read.table(fam_file, header=F, col.names=c('FID', 'IID', 'sire', 'dame', 'gender', 'pheno')) %>%
  as_tibble()
no_samples <- dim(fam_data)[1]

# Select the regions containing at least N valid markers for simulation
valid_scan_data <- scan_regions(regions = cds_uniq, freqs = afreq, min_nmrk = N)

############################## Simulate phenotypes
sim_i <- 1
sim_regions <- valid_scan_data[sample(1:nrow(valid_scan_data), size = N_sim, replace = T),]
tbl_colnames <- c('region', 'gene', 'num_valid_mrk', paste0('snp_', 1:N), paste0('ref',1:N), paste0('alt',1:N), paste0('eff_', 1:N), paste0('maf_',1:N))
metadata <- tbl_colnames %>% purrr::map_dfc(setNames, object = list(character()))
phenomatrix <- matrix(NA, nrow = length(fam_data$IID), ncol = N_sim + 1)

for (i in 1:nrow(sim_regions)) {
  sim_region <- sim_regions[i, ]
  cds <- cds_uniq %>% filter(gene == sim_region$gene)
  print(paste0("Processing ", cds$gene, " (round ", sim_i, " of ", N_sim,  ") ..."))
  focal_region <- paste0(cds$start, "-", cds$stop)
  valid_snps <- valid_markers_lst[[as.numeric(sim_region$index)]]
  metadata[sim_i, 'region'] <- focal_region
  metadata[sim_i, 'gene'] <- cds$gene
  metadata[sim_i, 'num_valid_mrk'] <- as.character(sim_region$num_valid_mrk)
  print(paste0("- I found ", sim_region$num_valid_mrk, " valid loci in ", cds$gene, " ", focal_region, "."))

  # Randomly select the loci to use
  selected_loci <- sample(valid_snps$pos, size = N, replace = F)
  selected_snps <- valid_snps %>% filter(pos %in% selected_loci) %>%
    arrange(match(pos, selected_loci))
  metadata[sim_i, 4:(3+N)] <- as.list(selected_loci)
  metadata[sim_i, (4+N):(3+N+N)] <- as.list(selected_snps$ref)
  metadata[sim_i, (4+N+N):(3+N+N+N)] <- as.list(selected_snps$alt)
  to_extract <- afreq[which(afreq$pos %in% as.numeric(selected_loci)), 'index'] %>% unlist()
  G <- seqminer::readPlinkToMatrixByIndex(plink_file, sampleIndex=seq(1:no_samples), markerIndex=to_extract)

  # Impute and fix encoding
  G <- array_branch(G, margin = 2) %>%
    # map(., fix_encoding) %>%
    map(., impute) %>%
    unlist(.) %>%
    matrix(., nrow = dim(G)[1])

  # Re-map genotype-phenotype
  G <- gp_remap(G, map=c(0, 1, 2))

  maf <- colSums(G)/(2*dim(G)[1])
  metadata[sim_i, (4+N+N+N):(3+N+N+N+N)] <- as.list(as.character(maf))

  # Randomly assign effect signs to loci
  eff_signs <- rep(1, times = N)
  if (N_neg > 0) {
    eff_signs[sample(x = 1:N, size = N_neg, replace = F)] <- -1
  }

  # Draw effects from a given distribution
  eff <- abs(rnorm(n = N, mean = eff_mu, sd = eff_sd))
  betas <- eff * eff_signs
  metadata[sim_i, (4+N+N+N+N):(3+N+N+N+N+N)] <- as.list(as.character(betas))

  # Add error term
  e <- rnorm(n = dim(G)[1], mean = e_mu, sd = e_sd)

  # Simulate phenotype
  y <- round(G %*% betas + e, digits = 3)

  ###### Store focal region metadata
  phenomatrix[,sim_i + 1] <- y
  sim_i <- sim_i + 1
}

colnames(phenomatrix) <- c('IID', paste0('y', 1:N_sim))
phenomatrix[,1] <- fam_data$IID


# Write the results
write.table(metadata, file=paste0(opt$`output-path`,'simulated_metadata.txt'), quote=F, row.names=F)
write.table(phenomatrix, file=paste0(opt$`output-path`, 'simulated_phenotypes.txt'), quote=F, row.names=F)
