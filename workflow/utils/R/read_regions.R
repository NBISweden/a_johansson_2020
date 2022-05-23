read_regions <- function(file, chr_prefix = "chr") {
  # Read regions specification, every region (randomly selected) can than be used for a separate simulation
  # Make sure regions that are present multiple times are collapsed to one that stretches from min to max coord
  # (e.g., from the beginning of exon1 to the end of the last exon)
  # Also, add computed region length to the table.
  cat("Reading regions from: ", file)
  cds <- read_delim(file, delim = '\t', col_names=c('chr', 'start', 'stop', 'gene'))
  cds_uniq <- cds %>% group_by(gene) %>%
    summarise(chr=chr, start=min(start), stop=max(stop)) %>%
    distinct() %>%
    mutate(chr = paste0(chr_prefix, chr)) %>%
    mutate(length = stop - start)
  return(cds_uniq)
}
