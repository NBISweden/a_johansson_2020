library(tidyverse)
library(seqminer)

my_chr <- 'chr22' # chromosome to run simulations on
path <- "/home/marcin/ExomeSeq2ndRel/marcin/"  # path to where the raw data files live
regions_file <- paste0(path, "GRCh38_chr22_cds.bed")
afreq_file <- paste0(path, "sntst_GRCh38_norel_rnd10000_chr22.afreq")
fam_file <- paste0(path, "sntst_GRCh38_norel_rnd10000_chr22.fam")
bim_file <- paste0(path, "sntst_GRCh38_norel_rnd10000_chr22.bim")
plink_file <- paste0(path, "sntst_GRCh38_norel_rnd10000_chr22")

threshold_low <- 0.05   # maf has to be above this threshold
threshold_upp <- 0.5    # maf has to be equal to or below this threshold

# Simulation parameters
N <- 1              # number of markers to be used per simulation
N_neg <- 0          # number of markers with negative effect size
eff_mu <- 1         # mean for effects distribution (effect sizes will be sampled from N(eff_mu, eff_sd))
eff_sd <- 0.01      # std. dev. for phenotype distribution
e_mu <- 0           # mean for residuals distribution
e_sd <- 0.05        # std. dev. for residuals distribution
N_sim <- 10         # number of traits to simulate

# Read regions specification, every region (randomly selected) can than be used for a separate simulation
# Make sure regions that are present multiple times are collapsed to one that stretches from min to max coord
# (e.g., from the beginning of exon1 to the end of the last exon)
# Also, add computed region length to the table.
cat("Reading regions from: ", regions_file)
cds <- read_delim(regions_file, delim = '\t', col_names=c('chr', 'start', 'stop', 'gene'))
cds_uniq <- cds %>% group_by(gene) %>%
  summarise(chr=chr, start=min(start), stop=max(stop)) %>%
  distinct() %>%
  mutate(chr = my_chr) %>%
  mutate(length = stop - start)

# Read index file and remove unwanted columns
#cat("Reading index from: ", bim_file)
#index <- read_delim(bim_file, delim = '\t', col_names = c('chr','to_remove','to_remove2','pos','ref','alt')) %>%
#  select(-to_remove, -to_remove2)

# Read allele frequencies from PLINK-generated file, compute maf in a new column,
# remove monomorphic markers
cat("Reading allele frequencies from: ", afreq_file)
afreq_all <- read_delim(afreq_file, delim = '\t') %>%
  separate(ID, into=c('chr','pos','ref','alt'), sep=':') %>%
  select(chr, pos, ref, alt, freq_alt=ALT_FREQS, no_obs=OBS_CT) %>%
  rowid_to_column("index") %>%
  mutate(maf = ifelse(freq_alt > 0.5, 1 - freq_alt, freq_alt)) %>%
  filter(maf > 0)

# Remove all monomorphic markers as well as all the markers that are
# above the maf upper threshold, i.e. the common markers
afreq <- afreq_all %>% filter(freq_alt > threshold_low & freq_alt <= threshold_upp) %>%
          filter(ref != "I" & ref != "D" & alt != "I" & alt != "D")

# Read the fam file
print(paste0("Reading family information from: ", fam_file))
fam_data <- read.table(fam_file, header=F, col.names=c('FID', 'IID', 'sire', 'dame', 'gender', 'pheno')) %>%
  as_tibble()
no_samples <- dim(fam_data)[1]

impute <- function(v, mode='naive') {
  # fix encoding by swaping codes if the alternative allele is not the minor one
  if (sum(v == 2) > sum(v == 0)) {
    w <- v
    w[v == 0] <- 2
    w[v == 2] <- 0
    v <- w
  }

  N <- 2 * sum(v %in% c(0,1,2))
  missing_ids <- which(!(v %in% c(0,1,2)))

  if (mode == 'HW') {
    f_a1 <- (2 * sum(v == 0) + sum(v == 1)) / N
    f_a2 <- (2 * sum(v == 2) + sum(v == 1)) / N
    p2 <- f_a1 * f_a1
    pq <- f_a1 * f_a2
    q2 <- f_a2 * f_a2
    probs <- c(p2,pq,q2)
    probs[is.na(probs)] <- 0
    if (all(probs) == 0) { probs = c(.33, .33, .33) }
    v[missing_idx] <- sample(c(0,1,2),
                             size = length(missing_idx),
                             replace = T,
                             prob = probs)
  } else if (mode == 'naive') {
    v[missing_ids] <- 2
  } else {
    warning('Imputation mode not supported. Assuming naive!')
    v[missing_ids] <- 2
  }
  return(v)
}

# Scan through all the regions and check which ones are valid, i.e. which ones have at least N markers that fulfill
# specified criteria and thus are suitable for simulation
num_regions <- dim(cds_uniq)[1]
scan_data <- data.frame(index=1:num_regions, gene=rep("",times=num_regions), num_valid_mrk=rep(0, times=num_regions))
valid_markers_lst <- vector("list", length = num_regions)
for (region_idx in 1:num_regions) {
  cds  <- cds_uniq[region_idx,]
  valid_markers <- afreq %>%
    filter(pos >= cds$start & pos <= cds$stop)
  scan_data[region_idx,] <- c(region_idx, cds$gene, dim(valid_markers)[1])
  valid_markers_lst[[region_idx]] <- valid_markers
}
scan_data <- scan_data %>% as_tibble()

valid_scan_data <- scan_data %>% filter(num_valid_mrk >= N)

############################## Simulate phenotypes
sim_i <- 1
sim_regions <- valid_scan_data[sample(1:nrow(valid_scan_data), size = N_sim, replace = T),]
tbl_colnames <- c('region', 'gene', 'num_valid_mrk', paste0('snp_', 1:N), paste0('eff_', 1:N), paste0('maf_',1:N))
metadata <- tbl_colnames %>% purrr::map_dfc(setNames, object = list(character()))
phenomatrix <- matrix(NA, nrow = length(fam_data$IID), ncol = N_sim + 1)

for (i in 1:nrow(sim_regions)) {
  sim_region <- sim_regions[i, ]
  cds <- cds_uniq %>% filter(gene == sim_region$gene)
  print(paste0("Pvrocessing ", cds$gene, " (round ", sim_i, " of ", N_sim,  ") ..."))
  focal_region <- paste0(cds$start, "-", cds$stop)
  valid_snps <- valid_markers_lst[[as.numeric(sim_region$index)]]
  metadata[sim_i, 'region'] <- focal_region
  metadata[sim_i, 'gene'] <- cds$gene
  metadata[sim_i, 'num_valid_mrk'] <- as.character(sim_region$num_valid_mrk)
  print(paste0("- I found ", sim_region$num_valid_mrk, " valid loci in ", cds$gene, " ", focal_region, "."))

  # Randomly select the loci to use
  selected_loci <- sample(valid_snps$pos, size = N, replace = F)
  metadata[sim_i, 4:(3+N)] <- as.list(selected_loci)
  to_extract <- afreq[which(afreq$pos %in% as.numeric(selected_loci)), 'index'] %>% unlist()
  G <- seqminer::readPlinkToMatrixByIndex(plink_file, sampleIndex=seq(1:no_samples), markerIndex=to_extract)

  # Impute and fix encoding
  G <- array_branch(G, margin = 2) %>%
    map(., impute) %>%
    unlist(.) %>%
    matrix(., nrow = dim(G)[1])

  # If needed, remap G to get a desired genotype-phenotype map, e.g.:
  # G[G == 1] <- 2
  # to say othat ne allele is enough (dominance model)
  maf <- colSums(G)/(2*dim(G)[1])
  metadata[sim_i, (4+N+N):(3+N+N+N)] <- as.list(as.character(maf))

  # Randomly assign effect signs to loci
  eff_signs <- rep(1, times = N)
  if (N_neg > 0) {
    eff_signs[sample(x = 1:N, size = N_neg, replace = F)] <- -1
  }

  # Draw effects from a given distribution
  eff <- abs(rnorm(n = N, mean = eff_mu, sd = eff_sd))
  betas <- eff * eff_signs
  metadata[sim_i, (4+N):(3+N+N)] <- as.list(as.character(betas))

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
write.table(metadata, file='simulated_metadata.txt', quote=F, row.names=F)
write.table(phenomatrix, file='simulated_phenotypes.txt', quote=F, row.names=F)