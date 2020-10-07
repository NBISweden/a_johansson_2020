#' @title Convert tibble with genotypes to GenABEL internal raw genotype format and writ it into a file.
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
#' @return Nothing is returned
#' @param x tibble with genotypes as described in the Details section
#' @param output name of the file to write output to
#' @param progress report progress every N individuals
#' @export
tibble_to_raw_genotypes <- function(x, output = 'genotypes.raw', progress = 1) {
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
  # output <- "genotypes.raw"
  # progress <- 1

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
  nbytes <- ceiling(num_ids/4)
  coding <- as.raw(rep(1, length(pos)))
  strand <- as.raw(rep(0, length(pos)))

  ofile <- file(output, 'w')
  cat(file = ofile, "#GenABEL raw data version 0.1\n")
  cat(file = ofile, ids, "\n")
  cat(file = ofile, marker_names, "\n")
  cat(file = ofile, chrom, "\n")
  cat(file = ofile, pos, "\n")
  cat(file = ofile, coding, "\n")
  cat(file = ofile, strand, "\n")
  rdta <- raw(nbytes)
  for (i in 1:n_ids) {
    #i <- 1
    gtin <- x %>% pull(ids[i])
    gchk <- (gtin == 0 | gtin == 1 | gtin == 2 | gtin == 3)
    if (!all(gchk)) {
      cat("Wrong genotype codes:\nCODE\tID\tSNP\n")
      wlst <- which(!gchk)
      for (j in 1:length(wlst)) cat(gtin[wlst[j]], "\t",
                                    ids[wlst[j]], "\t", marker_names[i], "\n")
      stop("execution terminated")
    }
    rdta <- GenABEL:::put.snps(gtin)
    cat(file = ofile, rdta, "\n")
    if (progress && round(i/progress) == i/progress)
      cat("Converted", i, "records...\n")
  }
  close(ofile)
}
