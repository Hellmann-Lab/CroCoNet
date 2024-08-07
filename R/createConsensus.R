#' Create consensus network
#'
#' Integrates networks across different clones and different species into a single consensus network. The input networks should all contain the same nodes (genes). For each edge, the consensus adjacencies are calculated as the weighted average of clonewise adjacencies using weights that correct for 1) the phylogenetic distances between species (if the phylogenetic tree is provided) and 2) the different numbers of clones per species.
#'
#' @param network_list A named list that contains the networks of each clone in an igraph format.
#' @param clone2species A data frame with columns 'clone' and 'species' that specifies which species each clone belongs to. The names of clones and species should match the names of 'network_list' and the tip labels of 'tree', respectively.
#' @param tree Object of class 'phylo' that gives the phylogenetic tree of the species.
#'
#' @return Consensus network in an igraph format.
#' @export
#'
#' @examples
#' consensus_network <- createConsensus(network_list_scaled_filt, clone2species, tree)
createConsensus <- function(network_list, clone2species, tree = NULL) {

  . = weight2 = factor = weight = clone = n_supporting_clones = supporting_clones = from = to = NULL # due to NSE notes in R CMD check

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

  dt_list <- convertToDT(network_list)

  if (is.null(tree)) {

    factors <- clone2species %>%
      dplyr::group_by(.data[["species"]]) %>%
      dplyr::mutate(n_clones = length(.data[["clone"]])) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(n_species = length(unique(.data[["species"]])),
                    factor = 1/(.data[["n_clones"]]*.data[["n_species"]])) %>%
      dplyr::select(.data[["clone"]], .data[["factor"]]) %>%
      data.table::as.data.table()

  } else {

    # phylogenetic similarity matrix
    sim.matrix = ape::cophenetic.phylo(tree) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("species1") %>%
      tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
      dplyr::right_join(clone2species %>% dplyr::rename(species1 = .data[["species"]], clone1 = .data[["clone"]]), by = "species1", relationship =
                   "many-to-many") %>%
      dplyr::right_join(clone2species %>% dplyr::rename(species2 = .data[["species"]], clone2 = .data[["clone"]]), by = "species2", relationship =
                   "many-to-many") %>%
      dplyr::mutate(similarity = 1 - .data[["distance"]]/max(.data[["distance"]]))

    # calculate factors (inversely proportional to the distance from the centroid)
    factors <- sim.matrix %>%
      dplyr::group_by(.data[["clone1"]]) %>%
      dplyr::summarise(factor = 1 / sum(.data[["similarity"]])) %>%
      dplyr::mutate(factor = .data[["factor"]] / sum(.data[["factor"]])) %>%
      dplyr::rename(clone = .data[["clone1"]]) %>%
      data.table::as.data.table()

  }

  consensus <- data.table::rbindlist(dt_list, idcol = "clone")[
    factors, on = "clone"][
      , weight2 := factor * weight][
        , .(weight = sum(weight2), n_supporting_clones = .N, supporting_clones = paste(clone, collapse = ",")), by = .(from, to)]

  igraph::graph_from_data_frame(consensus, directed = FALSE, vertices = data.frame(vertex = V(network_list[[1]])$name))

}
