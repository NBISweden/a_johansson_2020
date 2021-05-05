#' @title Read regions information from bed file.
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description  Read in regions from a bed file (no header allowed).
#' @details Only the first three columns will be read in.
#' chromosome names should begin with 'chr'. The returned coordinates are in the VCF coordinates (1-based).
#' @param bed_file path to the standard headerless bed file
#' @return a tibble with region name,chromosome, start coord, end coord (in VCF 1-based coords) and region size
#'
read_regions_bed <- function(bed_file) {
  bed_data <- read.table(bed_file, header = F) %>%
    as_tibble() %>%
    select(V1, V2, V3) %>%
    rename(chr = V1, start = V2, end = V3) %>%
    mutate(start = start + 1) %>%
    mutate(region = paste0(str_remove(chr, 'chr'), ':', start, '-', end)) %>%
    mutate(size = end - start) %>%
    relocate(region)
  return(bed_data)
}
