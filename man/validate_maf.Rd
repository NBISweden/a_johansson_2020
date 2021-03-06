% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_maf.R
\name{validate_maf}
\alias{validate_maf}
\title{Validate minor allele frequency vector}
\usage{
validate_maf(x)
}
\arguments{
\item{x}{maf vector}
}
\value{
a list with problematic elements, the reason why they
are problematic and fixed maf values (if possible). If no problems,
reason and problems will be NULL.
}
\description{
Validate minor allele frequency vector
}
\author{
Marcin Kierczak <marcin.kierczak@scilifelab.se>
}
