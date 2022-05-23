read_afreq <- function(file) {
  # Read allele frequencies from PLINK-generated file, compute maf in a new column,
  # remove monomorphic markers
  cat("Reading allele frequencies from: ", file)
  afreq_all <- read_delim(file, delim = '\t') %>%
    separate(ID, into=c('chr','pos','ref','alt'), sep=':') %>%
    select(chr, pos, ref, alt, freq_alt=ALT_FREQS, no_obs=OBS_CT) %>%
    rowid_to_column("index") %>%
    mutate(maf = ifelse(freq_alt > 0.5, 1 - freq_alt, freq_alt)) %>%
    filter(maf > 0)
  return(afreq_all)
}
