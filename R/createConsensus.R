#' Create consensus network
#'
#' Integrates networks across different clones and different species into a single consensus network in a phylogeny-aware manner.
#'
#' The input networks should all contain the same nodes (genes). The output consensus network contains the same nodes and all edges that were detected in at least 1 of the clones.
#'
#' For each edge, the consensus edge weight (adjacency) is calculated as the weighted mean of clonewise edge weights. The weighted mean corrects for 1) the phylogenetic distances between species (if the phylogenetic tree is provided) and 2) the different numbers of clones per species. As a result, the approach downweighs the edge weights of the clones that 1) belong to closely related species or 2) belong to species with many clones, so that an imbalanced sampling across the phylogenetic tree does not bias the consensus network.
#'
#' If an edge is not present in one of the clones, the edge weight in that clone is regarded as 0 for the calculation of the weighted mean. The number of clones and the names of clone where an edge was detected are saved as 2 new edge attributes ("n_supporting_clones" and "supporting_clones", respectively) in the output consensus network.
#'
#' In the next steps of the pipeline, this consensus network can be used to assign modules jointly for all species while avoiding species bias (see \code{\link{assignInitialModules}} and \code{\link{pruneModules}}).
#'
#' @param network_list A named list of \code{\link{igraph}} objects containing the networks of all clones.
#' @param clone2species A data frame specifying which species each clone belongs to, required columns:
#' \describe{
#' \item{clone}{Character, name of the clone.}
#' \item{species}{Character, name of the species.}
#' }
#' @param tree Object of class \code{\link{phylo}}, the phylogenetic tree of the species.
#'
#' @return Consensus network in an \code{\link{igraph}} format with the following edge attributes:
#' \describe{
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted average of clonewise edge weights.}
#' \item{n_supporting_clones}{Integer, the number of clones where the edge was detected.}
#' \item{supporting_clones}{Character, the list of clones where the edge was detected.}
#' }
#' @export
#'
#' @examples
#' consensus_network <- createConsensus(network_list, clone2species, tree)
createConsensus <- function(network_list, clone2species, tree = NULL) {

  # check input data
  if (!inherits(network_list, "list"))
    stop("The argument \"network_list\" should be a named list.")

  if (is.null(names(network_list)))
    stop("The argument \"network_list\" should be a named list.")

  if (any(!sapply(network_list, function(net) {inherits(net, "igraph")})))
    stop("All elements of \"network_list\" should be of class \"igraph\".")

  if (!is.data.frame(clone2species))
    stop("The argument \"clone2species\" should be a data frame.")

  if (any(!(c("clone", "species") %in% colnames(clone2species))))
    stop("The argument \"clone2species\" should contain the columns \"clone\" and \"species\".")

  if (any(sort(unique(clone2species$clone)) != sort(names(network_list))))
    stop("The names of \"network_list\" should match the clone names in the column \"clone\" of \"clone2species\".")

  if (!is.null(tree)) {

    if (any(!clone2species$species %in% tree$tip.label))
      stop("One or more species in \"clone2species\" are not present in the phylogenetic tree. Please make sure the the species names in the column \"species\" of \"clone2species\" match the tip labels of \"tree\".")

  }

  # avoid NSE notes in R CMD check
  . = weight2 = factor = weight = clone = n_supporting_clones = supporting_clones = from = to = NULL

  # convert igraphs to data tables
  dt_list <- convertToDT(network_list)

  # if a phylogenetic tree is not provided, correct only for the number of clones per species
  if (is.null(tree)) {

    # calculate the weight of each clone for the weighted mean (inversely proportional to the number of clones that belong to the same species)
    factors <- clone2species %>%
      dplyr::group_by(.data[["species"]]) %>%
      dplyr::mutate(n_clones = length(.data[["clone"]])) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(n_species = length(unique(.data[["species"]])),
                    factor = 1/(.data[["n_clones"]]*.data[["n_species"]])) %>%
      dplyr::select(.data[["clone"]], .data[["factor"]]) %>%
      data.table::as.data.table()

  # if a phylogenetic tree is provided, correct for both the number of clones per species and the phylogenetic distance between the species
  } else {

    # get phylogenetic similarity matrix based on the tree
    sim.matrix = ape::cophenetic.phylo(tree) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("species1") %>%
      tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
      dplyr::right_join(clone2species %>% dplyr::rename(species1 = .data[["species"]], clone1 = .data[["clone"]]), by = "species1", relationship =
                   "many-to-many") %>%
      dplyr::right_join(clone2species %>% dplyr::rename(species2 = .data[["species"]], clone2 = .data[["clone"]]), by = "species2", relationship =
                   "many-to-many") %>%
      dplyr::mutate(similarity = 1 - .data[["distance"]]/max(.data[["distance"]]))

    # calculate the weight of each clone for the weighted mean (inversely proportional to the sum of phylogenetic similarities between the given clone and all others)
    factors <- sim.matrix %>%
      dplyr::group_by(.data[["clone1"]]) %>%
      dplyr::summarise(factor = 1 / sum(.data[["similarity"]])) %>%
      dplyr::mutate(factor = .data[["factor"]] / sum(.data[["factor"]])) %>%
      dplyr::rename(clone = .data[["clone1"]]) %>%
      data.table::as.data.table()

    factors <- factors[match(clone2species$clone, factors$clone), ]

  }

  # calculate weighted mean per edge across all clones
  consensus <- data.table::rbindlist(dt_list, idcol = "clone")[
    factors, on = "clone"][
      , weight2 := factor * weight][
        , .(weight = sum(weight2), n_supporting_clones = .N, supporting_clones = paste(clone, collapse = ",")), by = .(from, to)][
          order(from, to)]

  # convert data table to igraph
  igraph::graph_from_data_frame(consensus, directed = FALSE, vertices = data.frame(vertex = V(network_list[[1]])$name))

}
