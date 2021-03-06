% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tibble_to_gwaa.R
\name{tibble_to_gwaa}
\alias{tibble_to_gwaa}
\title{Convert tibble with genotypes along with a vector of phenotypes to GenABEL internal raw gwaa.data format.}
\usage{
tibble_to_gwaa(x, sex, trait, progress = 1)
}
\arguments{
\item{x}{tibble with genotypes as described in the Details section}

\item{sex}{a vector of sex values: 0 - female, 1 - male}

\item{trait}{a vector of trait values}

\item{progress}{report progress every N individuals}
}
\value{
gwaa.data-class object
}
\description{
Convert tibble with genotypes along with a vector of phenotypes to GenABEL internal raw gwaa.data format.
}
\details{
Convert a tibble with the following shape:
test_data <- tibble(
  chr = c(1),
  pos = c(1001,1002,1005,1017),
  snp = c('snp1', 'snp2', 'snp3', 'snp4'),
  strand = c('u'),
  ind1 = c(0,1,1,2),
  id2 = c(2,NA,1,1),
  john = c(0,0,1,0),
  jill = c(1,1,1,2),
  elon = c(2,0,NA,1)
where columns from 5 to the end are genotypes and names of the columns are sample/individual names.
the `snp` column contains names of the snp markers,
`chr` is an integer denoting chromosome,
`pos` is the genome-wide position (or chromosome-based position),
`strand` is currently ignored, but required as a separate column.
Genotypes have to be coded as counts of the MINOR ALLELE:
AA - 0, Aa - 1 and aa - 2, NA values are allowed.
}
\author{
Marcin Kierczak <marcin.kierczak_NO_SPAM_scilifelab.se>
}
