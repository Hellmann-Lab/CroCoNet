#' Add directionality of interactions
#'
#' Determines the directionality of each network edge (positively correlated/coexpressed or negatively correlated/anti-coexpressed) based on a modified Spearman's correlation between the expression profiles of the 2 genes that form the edge.
#'
#' If the networks were inferred using a method that does not distinguish coexpressed and anti-coexpressed gene pairs, it might be useful to add this information for lines of analysis where it makes sense to separate the activated and repressed target genes of a regulator (e.g. for the calculation of eigengenes, see \code{\link{calculateEigengenes}}). If the network inference method output edges with both positive and negative edge weights in the first place (e.g. correlation-based methods), the edge attribute "direction" is already created during the step \code{\link{normalizeEdgeWeights}} and does not have to be calculated again.
#'
#' Here the directionality of a geneA-geneB edge refers to the characteristic whether geneA and geneB are coexpressed or anti-coexpressed and NOT whether geneA regulates geneB or geneB regulates geneA. The network remains undirected in a graph theoretical sense.
#'
#' The calculation of directionality relies on the approximate version of the Spearman's rho, significance testing and blocking implemented by \code{\link{correlatePairs}}. The results are summarized as 3 new edge attributes in the \code{\link{igraph}} object: rho (approximate Spearman's correlation coefficient), p.adj (BH-corrected approximate p-value) and direction ("+" or "-").
#'
#' @param network An \code{\link{igraph}} object containing the consensus network or the network of a clone.
#' @param sce A \code{\link{SingleCellExperiment}} object containing the expression data either for all clones (in case \code{network} is the consensus network) or for the clone of interest (in case \code{network} is a clonewise network). If \code{sce} contains the expression data of all clones, it is also expected to have a metadata column "clone" specifying which clone each cell belongs to; this will be used to define the blocking levels.
#' @param assay Character, the name of the assay in \code{sce} that should be used for the calculation of gene-gene correlations (default: "logcounts").
#' @param n_cores Integer, the number of cores (default: 1).
#'
#' @return An \code{\link{igraph}} object, the input \code{network} extended by the information about the direction of interactions. In addition to the original edge attributes, it contains 3 new attributes:
#' \describe{
#' \item{rho}{Numeric, the approximate Spearman's correlation coefficient between the expression profiles of the 2 genes that form the edge.}
#' \item{p.adj}{Numeric, BH-corrected approximate p-value of rho.}
#' \item{direction}{Character, the direction of the interaction between the 2 genes that form the edge ("+" = positively correlated/coexpressed or "-" = negatively correlated/anti-coexpressed).}
#' }
#'
#' @export
#'
#' @examples
#' consensus_network <- network_list %>%
#'  createConsensus(clone2species, tree) %>%
#'  addDirectionality(sce)
addDirectionality <- function(network, sce, assay = "logcounts", n_cores = 1L) {

  # check input data
  if (!inherits(network, "igraph"))
    stop("The argument \"network\" should be of class \"igraph\".")

  if (!inherits(sce, "SingleCellExperiment"))
    stop("The argument \"sce\" should be of class \"SingleCellExperiment\".")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") && !(inherits(n_cores, "numeric") && n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  if ((!inherits(assay, "character") || length(assay) != 1))
    stop("The argument \"assay\" should be a string.")

  if (!assay %in% names(SummarizedExperiment::assays(sce)))
    stop("The argument \"sce\" dowes not contain the specified assay.")

  # set to NULL due to NSE notes in R CMD check
  . = from = to = rho = p.adj = FDR = direction = p.value = NULL

  # convert igraph to data table
  dt <- igraph::as_data_frame(network, "edges") %>%
    data.table::as.data.table()

  # expression matrix
  expr <- SummarizedExperiment::assay(sce, assay)

  # connections in the network
  pairings <- as.matrix(dt[, .(from, to)])

  # block info (if present)
  block <- sce$clone
  if (!is.null(block)) block <- as.factor(block)

  # set seed and schedule the restoration of the old seed
  old_seed <- globalenv()$.Random.seed
  on.exit(suspendInterrupts({
    if (is.null(old_seed)) {
      rm(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      assign(".Random.seed", value = old_seed, envir = globalenv(), inherits = FALSE)
    }
  }), add = TRUE)
  set.seed(0)

  # calculate correlations
  corr = suppressWarnings(data.table::as.data.table(scran::correlatePairs(expr, pairings = pairings, block = block, BPPARAM = BiocParallel::MulticoreParam(workers = n_cores))))

  # add correlations to the network data table
  dt_withDir <- dt[corr, on = c(from = "gene1", to = "gene2")][
    , `:=` (p.adj = data.table::fifelse(is.na(FDR), 0, FDR),
            direction = data.table::fcase(rho > 0, "+",
                                          rho < 0, "-",
                                          default = NA_character_),
            p.value = NULL, FDR = NULL)]

  # convert data table back to igraph
  igraph::graph_from_data_frame(dt_withDir, directed = FALSE, vertices = data.frame(vertex = V(network)$name))

}
