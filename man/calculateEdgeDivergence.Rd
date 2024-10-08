% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotNetworks.R
\name{calculateEdgeDivergence}
\alias{calculateEdgeDivergence}
\title{Calculate the divergence of intramodular edges}
\usage{
calculateEdgeDivergence(
  module_names,
  pruned_modules,
  consensus_network,
  network_list,
  clone2species,
  N = Inf,
  n_cores = 1L
)
}
\arguments{
\item{module_names}{Character or character vector, the name(s) of the module(s) of interest.}

\item{pruned_modules}{Data frame of the pruned modules, required columns:
\describe{
\item{regulator}{Character, transcriptional regulator.}
\item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
}}

\item{consensus_network}{\code{\link{igraph}} object, the consensus network across all species and clones.}

\item{network_list}{A list of \code{\link{igraph}} objects containing the networks per clone.}

\item{clone2species}{A data frame that specifies which species each clone belongs to, required columns:
\describe{
\item{clone}{Character, name of the clone.}
\item{species}{Character, name of the species.}
}}

\item{N}{Integer, the number of strongest edges to subset. If set to Inf (default), all edges in the module are considered.}

\item{n_cores}{Integer, the number of cores (default: 1).}
}
\value{
A data frame of the selected edges with 5 columns:
\describe{
\item{regulator}{Character, transcriptional regulator.}
\item{from, to}{Character, the 2 member genes of the regulator's module that form the edge.}
\item{consensus_weight}{Numeric, consensus edge weight/adjacency (the weighted average of clonewise adjacencies).}
\item{f_statistic}{Numeric, measue of edge divergence. It is calculated as the F-statistic from the ANOVA of edge weights with species as groups.}
\item{p-value}{Numeric, the p-value of the F-statistic.}
}
}
\description{
Selects the top \code{N} strongest edges of the chosen module(s) based on the consensus network and calculates a divergence score for each edge.
}
\details{
In the CroCoNet approach, networks are reconstructed per clone and combined into a single phylogeny-aware consensus network which is the basis of the module assignment.

This function selects the intramodular edges of the input module(s) in the consensus network. If \code{N} is set to \code{Inf}, all intramodular edges are considered for the divergence calculation, if \code{N} is smaller than the module size, the edges are ordered by their consensus edge weight and only the top \code{N} edges are kept per module.

For each edge that was kept, an edge divergence score is calculated based on its edge weights in the networks of individual clones. The edge weights are compared across species using an ANOVA. and the F-statistic (i.e. the variation across the species means / variation within the species) and the p-value of this F-statistic are output as the measures of edge divergence.
}
\examples{
POU5F1_mod_edge_divergence <- calculateEdgeDivergence("POU5F1",
                                                                 pruned_modules,
                                                                 consensus_network,
                                                                 network_list,
                                                                 clone2species)
}
