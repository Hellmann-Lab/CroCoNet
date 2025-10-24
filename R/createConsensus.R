#' Create consensus network
#'
#' Integrates networks across different replicates and different species into a single consensus network in a phylogeny-aware manner.
#'
#' The input networks should all contain the same nodes (genes). The output consensus network contains the same nodes and all edges that were detected in at least 1 of the replicates.
#'
#' For each edge, the consensus edge weight (adjacency) is calculated as the weighted mean of replicate-wise edge weights. The weighted mean corrects for 1) the phylogenetic distances between species (if the phylogenetic tree is provided) and 2) the different numbers of replicates per species. As a result, the approach downweighs the edge weights of the replicates that 1) belong to closely related species or 2) belong to species with many replicates, so that an imbalanced sampling across the phylogenetic tree does not bias the consensus network.
#'
#' If an edge is not present in one of the replicates, the edge weight in that replicate is regarded as 0 for the calculation of the weighted mean. If there are more than 10\% such missing/zero-weight edges, the number of replicates and the names of replicate where an edge was detected are saved as 2 new edge attributes ("n_supporting_replicates" and "supporting_replicates", respectively) in the output consensus network.
#'
#' In the next steps of the pipeline, this consensus network can be used to assign modules jointly for all species while avoiding species bias (see \code{\link{assignInitialModules}} and \code{\link{pruneModules}}).
#'
#' @param network_list A named list of \code{\link{igraph}} objects containing the networks of all replicates.
#' @param replicate2species A data frame specifying which species each replicate belongs to, required columns:
#' \describe{
#' \item{replicate}{Character, name of the replicate.}
#' \item{species}{Character, name of the species.}
#' }
#' @param tree Object of class \code{\link{phylo}}, the phylogenetic tree of the species.
#'
#' @return Consensus network in an \code{\link{igraph}} format with the following edge attributes:
#' \describe{
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted average of replicate-wise edge weights.}
#' \item{n_supporting_replicates}{Integer, the number of replicates where the edge was detected (only added if more than 10\% of the edges are not present in all replicates).}
#' \item{supporting_replicates}{Character, the list of replicates where the edge was detected (only added if more than 10\% of the edges are not present in all replicates).}
#' }
#' @export
#'
#' @examples
#' consensus_network <- createConsensus(network_list, replicate2species, tree)
createConsensus <- function(network_list, replicate2species, tree = NULL) {

  # check input data
  if (!inherits(network_list, "list"))
    stop("The argument \"network_list\" should be a named list.")

  if (is.null(names(network_list)))
    stop("The argument \"network_list\" should be a named list.")

  if (any(!sapply(network_list, function(net) {inherits(net, "igraph")})))
    stop("All elements of \"network_list\" should be of class \"igraph\".")

  if (!is.data.frame(replicate2species))
    stop("The argument \"replicate2species\" should be a data frame.")

  if (any(!(c("replicate", "species") %in% colnames(replicate2species))))
    stop("The argument \"replicate2species\" should contain the columns \"replicate\" and \"species\".")

  if (any(sort(unique(replicate2species$replicate)) != sort(names(network_list))))
    stop("The names of \"network_list\" should match the replicate names in the column \"replicate\" of \"replicate2species\".")

  if (!is.null(tree)) {

    if (!inherits(tree, "phylo"))
      stop("The argument \"tree\" should be a phylo object.")

    if (any(!replicate2species$species %in% tree$tip.label))
      stop("One or more species in \"replicate2species\" are not present in the phylogenetic tree. Please make sure the the species names in the column \"species\" of \"replicate2species\" match the tip labels of \"tree\".")

  }

  # avoid NSE notes in R CMD check
  . = phylo_corr_weight = phylo_weight_corr_dir = factor = weight = direction = replicate = n_supporting_replicates = supporting_replicates = from = to = NULL

  # convert igraphs to data tables
  dt_list <- convertToDT(network_list)

  # if a phylogenetic tree is not provided, correct only for the number of replicates per species
  if (is.null(tree)) {

    # calculate the weight of each replicate for the weighted mean (inversely proportional to the number of replicates that belong to the same species)
    factors <- replicate2species %>%
      dplyr::group_by(.data[["species"]]) %>%
      dplyr::mutate(n_replicates = length(.data[["replicate"]])) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(n_species = length(unique(.data[["species"]])),
                    factor = 1/(.data[["n_replicates"]]*.data[["n_species"]])) %>%
      dplyr::select(.data[["replicate"]], .data[["factor"]]) %>%
      data.table::as.data.table()

  # if a phylogenetic tree is provided, correct for both the number of replicates per species and the phylogenetic distance between the species
  } else {

    # get phylogenetic similarity matrix based on the tree
    sim.matrix = ape::cophenetic.phylo(tree) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("species1") %>%
      tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
      dplyr::right_join(replicate2species %>% dplyr::rename(species1 = .data[["species"]], replicate1 = .data[["replicate"]]), by = "species1", relationship =
                   "many-to-many") %>%
      dplyr::right_join(replicate2species %>% dplyr::rename(species2 = .data[["species"]], replicate2 = .data[["replicate"]]), by = "species2", relationship =
                   "many-to-many") %>%
      dplyr::mutate(similarity = 1 - .data[["distance"]]/max(.data[["distance"]]))

    # calculate the weight of each replicate for the weighted mean (inversely proportional to the sum of phylogenetic similarities between the given replicate and all others)
    factors <- sim.matrix %>%
      dplyr::group_by(.data[["replicate1"]]) %>%
      dplyr::summarise(factor = 1 / sum(.data[["similarity"]])) %>%
      dplyr::mutate(factor = .data[["factor"]] / sum(.data[["factor"]])) %>%
      dplyr::rename(replicate = .data[["replicate1"]]) %>%
      data.table::as.data.table()

    factors <- factors[match(replicate2species$replicate, factors$replicate), ]

  }

  # combine all replicates and add the weights for the weighted mean
  dt_combined <- data.table::rbindlist(dt_list, idcol = "replicate")[factors, on = "replicate"]

  # are there more than 10% edges not present in all replicates?
  n_genes <- igraph::vcount(network_list[[1]])
  n_edges <- nrow(dt_combined)
  n_replicates <- length(dt_list)
  n_edges_expected <- n_genes * (n_genes - 1) / 2 * n_replicates
  full <- n_edges > 0.9*n_edges_expected

  # calculate consensus edge weight and consensus direction (if the column is present) per edge across all replicates
  if ("direction" %in% colnames(dt_combined)) {

    if (full) {

      consensus <- dt_combined[
        , `:=`(phylo_corr_weight = factor * weight,
               phylo_weight_corr_dir = factor * weight * ifelse(direction == "+", 1, -1))][
                 , .(weight = sum(phylo_corr_weight),
                     direction = ifelse(sum(phylo_weight_corr_dir) > 0, "+", "-")),
                 by = .(from, to)][
                   order(from, to)]

    } else {

      consensus <- dt_combined[
        , `:=`(phylo_corr_weight = factor * weight,
               phylo_weight_corr_dir = factor * weight * ifelse(direction == "+", 1, -1))][
                 , .(weight = sum(phylo_corr_weight),
                     direction = ifelse(sum(phylo_weight_corr_dir) > 0, "+", "-"),
                     n_supporting_replicates = .N,
                     supporting_replicates = paste(replicate, collapse = ",")),
                 by = .(from, to)][
                   order(from, to)]

    }

  } else {

    if (full) {

      consensus <- dt_combined[
        , phylo_corr_weight := factor * weight][
          , .(weight = sum(phylo_corr_weight)),
          by = .(from, to)][
            order(from, to)]

    } else {

      consensus <- dt_combined[
        , phylo_corr_weight := factor * weight][
          , .(weight = sum(phylo_corr_weight),
              n_supporting_replicates = .N,
              supporting_replicates = paste(replicate, collapse = ",")),
          by = .(from, to)][
            order(from, to)]

    }

  }

  # convert data table to igraph
  igraph::graph_from_data_frame(consensus, directed = FALSE, vertices = data.frame(vertex = V(network_list[[1]])$name))

}
