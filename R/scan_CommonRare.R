#' @title Run SKAT-O CommonRare model on all regions
#' @author Marcin Kierczak <marcin.kierczak@scilifelab.se>
#' @description Given a phenotype simulated for a focal region, scan for
#' associations in all supplied regions.
#' @param focal_region - name of the focal region
#' @param regions - list of regions to scan
#' @param null_model - SKAT null model
#' @param vcf_file - path to vcf file with genotypes
#' @return a tibble with association p-value for every region and info whether
#' the region is the focal region
#' @export
#'
scan_CommonRare <- function(focal_region, regions, null_model, vcf_file) {
  regions <- regions$region
  tmp <- foreach(x = regions) %dopar% {
    is_focal <- F
    if (x == focal_region) {
      is_focal <- T
    }
    gt <- read_region_vcf(x, vcf_file)
    if (dim(gt$G)[2] > 0) { # some regions are not genotyped and need to be skipped
      skat <- SKAT::SKAT_CommonRare(Z = gt$G, obj = null_model)
      tmp <- list(region = x, p_value = skat$p.value, is_focal = is_focal)
    } else {
      tmp <- list(region = x, p_value = NA, is_focal = is_focal)
    }
    tmp
  }
  tib <- map(tmp, ~data.frame(.)) %>%
    map_dfr(~mutate_all(., as.character)) %>%
    as_tibble() %>%
    mutate_at(2, as.double) %>%
    mutate_at(3, as.logical)
  return(tib)
}
