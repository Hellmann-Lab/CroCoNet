% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalizeEdgeWeights.R
\name{normalizeEdgeWeights}
\alias{normalizeEdgeWeights}
\alias{normaliseEdgeWeights}
\title{Normalize edge weights between 0 and 1}
\usage{
normalizeEdgeWeights(
  network_list,
  signed = FALSE,
  min_weight = NULL,
  max_weight = NULL,
  n_cores = 1L
)
}
\arguments{
\item{network_list}{A named list of \code{\link{igraph}} objects containing the networks of all clones.}

\item{signed}{Logical indicating whether a signed network is desired (default: FALSE, see \code{Details}).}

\item{min_weight, max_weight}{Numeric, the theoretical minimum and maximum values of the edge weights. If set to NULL (default), the normalization is performed using the empirical minimum and maximum. For correlation-based edge weights, please set \code{min_weight} and \code{max_weight} to -1 and 1, respectively.}

\item{n_cores}{Integer, the number of cores (default: 1).}
}
\value{
A named list of \code{\link{igraph}} objects containing the networks of all clones, with the edge weights normalized between 0 and 1. If originally there were both positive and negative edge weights present in the data, a new edge attribute is added to all \code{\link{igraph}} objects:
\describe{
\item{direction}{Character, the direction of the interaction between the 2 genes that form the edge ("+" or "-" depending on the sign of the edge weight before normalization).}
}
}
\description{
Transforms the edge weights of all networks to be in the range [0, 1].
}
\details{
Normalizing the edge weights between 0 and 1 makes them interpretable as adjacencies and ensures that network concepts such as connectivity are applicable.

There are 2 approaches for the normalization:
\itemize{
 \item{If \code{signed} is set to FALSE (unsigned network, the default), gene pairs with high negative edge weights are considered as connected as gene pairs with high positive edge weights. Therefore the negative edge weights are first replaced by their absolute values, then all edge weights are scaled by the maximum weight across all networks: \deqn{w_{new} =  \frac{|w|}{max(|w|)}} After the transformation, the edge weights/adjacencies around 0 correspond to the former low positive and low negative values, while the edge weights/adjacencies around 1 correspond to the former high positive and high negative values.}
 \item{If \code{signed} is set to TRUE (signed network), gene pairs with high negative edge weights are considered unconnected. Therefore all edge weights are transformed between 0 and 1 using a min-max normalization: \deqn{w_{new} = \frac{w - min(w)}{max(w) - min(w)}} After the transformation, the edge weights/adjacencies around 0 correspond to the former high negative values and the edge weights/adjacencies around 1 correspond to the former high positive values.}
 }

If the theoretical minimum and maximum edge weights are known, these can be provided using the parameters \code{min_weight} and \code{max_weight}. For example, if the networks were inferred by calculating correlations between gene expression profiles, \code{min_weight} should be set to -1 and \code{max_weight} should be set to 1. If \code{min_weight} and \code{max_weight} are left at NULL, the minimum and maximum edge weights are calculated empirically using the data.

After normalization by either of the approaches above, it is not possible to tell anymore which edges had a positive and which edges had a negative weight originally. Since this can be a useful piece of information (it can specify the mode of regulation: activation or repression), the information about the sign is stored as a new edge attribute "direction" in the output \code{\link{igraph}} objects ("+" if the original edge weight was positive and "-" is the original edge weight was negative). If all original edge weights were positive, the edge attribute "direction" is not added.
}
\examples{
network_list_norm <- normalizeEdgeWeights(network_list_raw)
}
\seealso{
\href{https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/}{"Signed or unsigned: which network type is preferable?" by Peter Langfelder}
}
