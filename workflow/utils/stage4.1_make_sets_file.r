library(tidyverse)
library(optparse)
source('R/read_regions.R')

#setwd('~/Dropbox/WABI/Projects/Johansson_PP/region_file_dev/')
#afreq_filename <- "reginfo/sntst_GRCh38_norel_rnd1000_chr22.afreq"
#regions_filename <- "reginfo/GRCh38_cds.bed"
#output_filename <- "regions.txt"
#progress_step <- 5 

option_list = list(
  make_option(c("--afreq-filename"), type="character", default=".",
              help="path to file with allele frequencies"),
  make_option(c("--regions-filename"), type="character", default=".",
              help="path to bed file with region definitions"),
  make_option(c("--output-filename"), type="character", default=".",
              help="path to output file"),
  make_option(c("--chr"), type="character", default="chr22",
              help="chromosome to make gene regions file for"),
  make_option(c("--progress-step"), type="numeric", default = 0,
              help = "show progress every N regions, if N=0 no progress bar will be shown"),
);

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
afreq_filename <- opt$`afreq-filename`
regions_filename <- opt$`regions-filename`
output_filename <- opt$`output-filename`
progress_step <- opt$`progress-step`  
my_chr <- opt$`chr`

afreq_data <- read.delim(afreq_filename) %>%
              select(ID) %>%
              mutate(ID, tmp=ID) %>%
              separate(tmp, into=c('chr','pos','ref','alt'), sep=':')

regions <- read_regions(regions_filename, chr_prefix = '')  
regions <- filter(regions, chr == my_chr)
if (progress_step > 0) {
  progress_thr <- floor((progress_step * nrow(regions)) / 100)
}

for (i in 1:nrow(regions)) {
  if (progress_step > 0) {
    if (i %% progress_thr == 0) {
      cat(paste0(i, "/", nrow(regions), " done..."))
    }
  }
  region <- regions[i,]
  set_ids <- afreq_data %>% 
    filter(paste0('chr', chr) == region$chr) %>%
    filter(pos >= region$start & pos <= region$stop) %>%
    select(ID) %>% 
    as_vector()
  if (length(set_ids) > 0) {
    row <- paste0("set_", region$gene, "\t", "var", "\t", paste0(set_ids, collapse="\t"))
    cat(paste0(row, '\n'), file = output_filename, append = T)
    row2 <- paste0("set_", region$gene, "\t", "anno", "\t", paste0(rep('x', times=length(set_ids)), collapse="\t"))
    cat(paste0(row2, '\n'), file = output_filename, append = T)
  }
}
