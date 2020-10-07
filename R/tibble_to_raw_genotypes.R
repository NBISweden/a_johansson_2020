#' @title Convert tibble with genotypes to GenABEL internal raw genotype format and writ it into a file.
#' @author Marcin Kierczak <marcin.kierczak_NO_SPAM_scilifelab.se>
#' @details Convert a tibble with the following shape:
#'  test_data <- tibble(
#'     id = c('indiv_1', 'indiv_2', 'indiv_3', 'indiv_4'),
#'     chr = c(1),
#'     pos = c(1001,1002,1005,1017),
#'     strand = c('u'),
#'     snp1 = c(0,1,1,2),
#'     snp2 = c(2,NA,1,1),
#'     snp3 = c(0,0,1,0),
#'     snp4 = c(1,1,1,2)
#'   )
#' where columns from 5 to the end are genotypes and names of the columns are marker names.
#' the `id` column contains names of the individuals,
#' `chr` is an integer denoting chromosome,
#' `pos` is the genome-wide position (or chromosome-based position),
#' `strand` is currently ignored, but required as a separate column.
#' Genotypes have to be coded as counts of the MINOR ALLELE:
#' AA - 0, Aa - 1 and aa - 2, NA values are allowed.
#' @return Nothing is returned
#' @param x tibble with genotypes as described in the Details section
#' @param output name of the file to write output to
#' @param progress report progress every N individuals
#' @export
tibble_to_raw_genotypes <- function(x, output = 'genotypes.raw', progress = 1) {
  #   test_data <- tibble(
  #     id = c('indiv_1', 'indiv_2', 'indiv_3', 'indiv_4'),
  #     chr = c(1),
  #     pos = c(1001,1002,1005,1017),
  #     strand = c('u'),
  #     snp1 = c(0,1,1,2),
  #     snp2 = c(2,NA,1,1),
  #     snp3 = c(0,0,1,0),
  #     snp4 = c(1,1,1,2)
  #   )
  #
  # x <- test_data
  # output <- "genotypes.raw"

  ids <- x$id
  nids <- length(x$id)
  mnams <- colnames(x)[-c(1:4)]
  # TODO: ADD CHECKS IN FUTURE!!!
  x <- x %>%
    mutate_at(.vars = mnams, function(x) { x + 1} ) %>%
    mutate_at(.vars = mnams, replace_na, 0)
  nmrk <- length(mnams)
  chrom <- x$chr
  pos <- x$pos
  nbytes <- ceiling(nids/4)
  coding <- as.raw(rep(1, length(pos)))
  strand <- as.raw(rep(0, length(pos)))

  ofile <- file(output, 'w')
  cat(file = ofile, "#GenABEL raw data version 0.1\n")
  cat(file = ofile, ids, "\n")
  cat(file = ofile, mnams, "\n")
  cat(file = ofile, chrom, "\n")
  cat(file = ofile, pos, "\n")
  cat(file = ofile, coding, "\n")
  cat(file = ofile, strand, "\n")
  rdta <- raw(nbytes)
  for (i in 1:nmrk) {
    #i <- 1
    gtin <- x %>% pull(mnams[i])
    gchk <- (gtin == 0 | gtin == 1 | gtin == 2 | gtin == 3)
    if (!all(gchk)) {
      cat("Wrong genotype codes:\nCODE\tID\tSNP\n")
      wlst <- which(!gchk)
      for (j in 1:length(wlst)) cat(gtin[wlst[j]], "\t",
                                    ids[wlst[j]], "\t", mnams[i], "\n")
      stop("execution terminated")
    }
    rdta <- GenABEL:::put.snps(gtin)
    cat(file = ofile, rdta, "\n")
    if (progress && progress(i/progress) == i/progress)
      cat("Converted", i, "records...\n")
  }
  close(ofile)
}
