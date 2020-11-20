#' @title Convert tibble with genotypes along with a vector of phenotypes to GenABEL internal raw gwaa.data format.
#' @author Marcin Kierczak <marcin.kierczak_NO_SPAM_scilifelab.se>
#' @details Convert a tibble with the following shape:
#' test_data <- tibble(
#'   chr = c(1),
#'   pos = c(1001,1002,1005,1017),
#'   snp = c('snp1', 'snp2', 'snp3', 'snp4'),
#'   strand = c('u'),
#'   ind1 = c(0,1,1,2),
#'   id2 = c(2,NA,1,1),
#'   john = c(0,0,1,0),
#'   jill = c(1,1,1,2),
#'   elon = c(2,0,NA,1)
# )
#' where columns from 5 to the end are genotypes and names of the columns are sample/individual names.
#' the `snp` column contains names of the snp markers,
#' `chr` is an integer denoting chromosome,
#' `pos` is the genome-wide position (or chromosome-based position),
#' `strand` is currently ignored, but required as a separate column.
#' Genotypes have to be coded as counts of the MINOR ALLELE:
#' AA - 0, Aa - 1 and aa - 2, NA values are allowed.
#' @return gwaa.data-class object
#' @param x tibble with genotypes as described in the Details section
#' @param sex a vector of sex values: 0 - female, 1 - male
#' @param trait a vector of trait values
#' @param progress report progress every N individuals
#' @import tidyverse magrittr
#' @export
#'
tibble_to_gwaa <- function(x, sex, trait, progress = 1) {
  #   test_data <- tibble(
  #     chr = c(1),
  #     pos = c(1001,1002,1005,1017),
  #     snp = c('snp1', 'snp2', 'snp3', 'snp4'),
  #     strand = c('u'),
  #     ind1 = c(0,1,1,2),
  #     id2 = c(2,NA,1,1),
  #     john = c(0,0,1,0),
  #     jill = c(1,1,1,2),
  #     elon = c(2,0,NA,1)
  #   )
  #
  # x <- test_data
  # progress <- 1
  # sex <- c(0,1,1,0,1)
  # trait <- c(1.3,3.4,4.4,2.2,1.7)

  ids <- colnames(x)[-c(1:4)]
  n_ids <- length(ids)
  marker_names <- x$snp
  # TODO: ADD CHECKS IN FUTURE!!!
  x <- x %>%
    mutate_at(.vars = ids, function(x) { x + 1} ) %>%
    mutate_at(.vars = ids, replace_na, 0)
  n_markers <- length(marker_names)
  chrom <- x$chr
  pos <- x$pos
  nbytes <- ceiling(n_ids/4)
  coding <- as.raw(rep(1, length(pos)))
  strand <- as.raw(rep(0, length(pos)))
  class(coding) <- 'snp.coding'
  class(strand) <- 'snp.strand'

  rdlst <- vector(mode = "list", length = n_markers)
  rdta <- raw(nbytes)
  i <- 1
  for (marker in marker_names) {
    gtin <- x %>% filter(snp == marker) %>% dplyr::select(all_of(ids))
    gchk <- (gtin == 0 | gtin == 1 | gtin == 2 | gtin == 3)
    if (!all(gchk)) {
      cat("Wrong genotype codes:\nCODE\tID\tSNP\n")
      wlst <- which(!gchk)
      for (j in 1:length(wlst)) cat(gtin[wlst[j]], "\t",
                                    ids[wlst[j]], "\t", marker_names[i], "\n")
      stop("execution terminated")
    }
    #start <- (nbytes * (i - 1)) + 1
    #stop <- start + nbytes - 1
    rdta <- GenABEL:::put.snps(gtin)
    rdlst[[i]] <- rdta
    if (progress && round(i/progress) == i/progress)
      cat("Converted", i, "markers...\n")
    i <- i + 1
  }
  pheno_dta = data.frame(id = ids, sex = sex, y = trait)
  g <- snp.data(nids = n_ids, rawdata = unlist(rdlst), idnames = ids,
           snpnames = marker_names,
           chromosome = as.factor(chrom),
           map = pos, coding = coding, strand = strand, male = sex)
  output <- new("gwaa.data", phdata = pheno_dta, gtdata = g)
  rm(g, pheno_dta)
  gc(verbose = FALSE)
  return(output)
  }
