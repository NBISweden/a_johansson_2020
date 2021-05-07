fp_CommonRare <- function(focal_region, regions, null_model, vcf_file) {
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
