#' @title Process focal region to see whether it is valid
#' @description Given a list of regions, the function randomly selects one, checks whether it is valid and
#' returns the updated list of valid regions and the genotyp
#' @author Marcin Kierczak <marcin.kierczak__INSERT_AT__scilifelab.se>
#' @param  code the code to be silenced
#' @return silenced result of code execution
#'
process_focal_region <- function(vcf_file, regions, N_markers, valid_regions, arch_fun) {
  is_selected <- FALSE
  while (!is_selected) {
    # Randomly draw a region to base simulation upon:
    # select random region from bed by shuffling valid regions and taking the first one
    set.seed(NULL)
    rnd <- sample(which(valid_regions == 1), size = 1)
    locus <- regions[rnd]

    # Extract the region
    focal_region <- hush(seqminer::readVCFToMatrixByRange(vcf_file, range = locus))[[1]]
    G <- focal_region %>%
      t() %>%
      gwasim::impute_G() %>%
      gwasim::fix_allele_encoding()

    if (!('gwasim' %in% class(G)))
      class(G) <- c(class(G), 'gwasim')

    maf <- gwasim::get_maf(G)
    #effects <- gwasim::get_effects(maf=maf, thr=0.01, N=N_markers, get_betas_fun = dbeta,
    #get_betas_args = list(shape1 = .1, shape2 = .1), rare=T, frac_negative=0, seed = 2)
    effects <- arch_fun(maf = maf, N = N_markers)

    # Filter: Check if the selected region contains at least N_markers rare markers
    if (length(effects$marker_idx) == N_markers) {
      is_selected = TRUE
    } else {
      print(paste0('Marking ', locus, ' invalid.'))
      valid_regions[rnd] <- 0
    }
  }
  result <- list(region = locus, maf = maf, effects = effects, G = G, valid_regions = valid_regions)
  return(result)
}
