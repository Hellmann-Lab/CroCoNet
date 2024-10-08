% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pruneModules.R
\name{pruneModules}
\alias{pruneModules}
\title{Prune modules}
\usage{
pruneModules(
  initial_modules,
  method = c("UIK_adj", "UIK_adj_kIM", "topN"),
  consensus_network = NULL,
  min_module_size = 20L,
  max_frac_modules_lost = 0.02,
  exponent = 1L,
  N = 50L
)
}
\arguments{
\item{initial_modules}{Data frame of initial modules, required columns:
\describe{
\item{regulator}{Character, transcriptional regulator.}
\item{target}{Character, member gene of the regulator's initial module.}
\item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
}}

\item{method}{Character, the pruning method, one of "UIK_adj", "UIK_adj_kIM", "topN".}

\item{consensus_network}{\code{\link{igraph}} object, the consensus network across all species and clones.}

\item{min_module_size}{Integer, the lower threshold of module size in case of the methods "UIK_adj" and "UIK_adj_kIM" (default: 20).}

\item{max_frac_modules_lost}{Numeric, the upper threshold of the fraction of removed modules in case of the methods "UIK_adj" and "UIK_adj_kIM" (default: 0.02).}

\item{exponent}{Integer, the exponent the regulator-target adjacency and intramodular connectivity is raised to the power of during the cumulative sum curve calculation in case of the methods "UIK_adj" and "UIK_adj_kIM" (default: 1, i.e. the regulator-target adjacencies and intramodular connectivities stay unchanged).}

\item{N}{An integer or a named integer vector specifying the desired pruned module size(s) in case of the method "topN" (default: 50).}
}
\value{
Data frame of the pruned modules with the following columns:
\describe{
\item{regulator}{Character, transcriptional regulator.}
\item{module_size}{Integer, the numer of genes assigned to a regulator.}
\item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
\item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
}
Additional columns present in \code{initial_modules} will also be preserved in \code{pruned_modules}.
}
\description{
Prunes the initial modules by keeping only the best targets of each transcriptional regulator.
}
\details{
3 methods are implemented to choose the best targets:
\itemize{
 \item{topN: Takes a fixed number of targets per regulator with the highest regulator-target adjacencies (for details see \code{\link{pruneModules_topN}}).}
 \item{UIK_adj: Applies a dynamic stepwise pruning based on the regulator-target adjacencies (for details see \code{\link{pruneModules_UIK_adj}}).}
 \item{UIK_adj_kIM: Applies a dynamic stepwise pruning based on the regulator-target adjacencies and intramodular connectivities (for details see \code{\link{pruneModules_UIK_adj_kIM}}).}
 }
}
\examples{
pruned_modules <- pruneModules(initial_modules, "topN", N = 30)
pruned_modules <- pruneModules(initial_modules, "UIK_adj")
pruned_modules <- pruneModules(initial_modules, "UIK_adj_kIM", consensus_network)
}
