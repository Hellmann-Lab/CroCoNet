% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{dist}
\alias{dist}
\title{Distance measures of the pruned modules}
\format{
A named list with 12 elements containing the distance measures per module. Each element is a data frame with 72 rows and 7 columns:
\describe{
\item{regulator}{Character, transcriptional regulator.}
\item{module_size}{Integer, the numer of target genes assigned to a regulator.}
\item{clone1, clone2}{Character the names of the clones compared.}
\item{species1, species2}{Character, the names of the species \code{clone1} and \code{clone2} belong to, respectively.}
\item{dist}{Numeric, distance measure ranging from 0 to 1, calculated based on the correlation of intramodular connectivities.}
}
}
\usage{
dist
}
\description{
Distance measures per clone pair for the original (i.e. not jackknifed) pruned modules.
}
\details{
The distance measures were calculated based on the correlation of intramodular connectivities: \deqn{dist = \frac{1 - cor.kIM}{2}} (a correlation of 1 corresponds to a distance of 0, whereas a correlation of -1 corresponds to a distance of 1). Each element of the list corresponds to a module and contains the distance measures between all possible pairs of clones for this module.
}
\keyword{datasets}
