#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


#' Convert a list of data frames to a list of igraphs
#'
#' @description Takes care of the conversion between data table and igraph formats of the networks.
#' @param dt_list A list of data tables that contain the columns 'from', 'to' and 'weight'. Each row should correspond to an edge in the network with 'from' and 'to' as the end nodes and 'weight' as the edge weight. Additional columns will be converted to edge attributes in the igraphs.
#' @param network_genes Character vector of all genes in the network.
#' @param n_cores Number of cores.
#'
#' @return A list of igraphs.
#' @export
convertToGraph <- function(dt_list, network_genes, n_cores = 1L) {

  dt = NULL

  doParallel::registerDoParallel(n_cores)

  network_list <- foreach::foreach(dt = dt_list) %dopar% {

                                   igraph::graph_from_data_frame(dt, directed = FALSE, vertices = data.frame(vertex = network_genes))

                                 }

  doParallel::stopImplicitCluster()

  names(network_list) <- names(dt_list)

  network_list

}


#' Convert a list of igraphs to a list of data tables
#'
#' @description Takes care of the conversion between the igraph and data table formats of the networks.
#' @param network_list A list of igraph objects
#'
#' @return A list of data tables with the columns 'from', 'to', 'weight' and any additional edge attributes the input igraphs contain.
#' @export
convertToDT <- function(network_list) {

  dt_list <- lapply(network_list, function(network) {

    igraph::as_data_frame(network, "edges") %>%
      data.table::as.data.table()

  })

}


#' Convert a tree to a data frame
#'
#' @description Converts a tree to a data frame where each row corresponds to a branch in the tree.
#' @param tree Object of class [phylo] with clones as tips. It is expected to have a component 'species' that specifies which species each tip belongs to.
#'
#' @return A data frame of tree branches.
#' @export
getTreeDf <- function(tree) {

  # make sure that the edges are in the same order in the data frame as in the tree (will make life easier later on)
  data.frame(branch.length = tree$edge.length) %>%
    dplyr::left_join(tidytree::as_tibble(tree) %>%
                       tidyr::drop_na(.data[["branch.length"]]),
                     by = "branch.length",
                     relationship = "many-to-many") %>%
    dplyr::distinct() %>%
    dplyr::left_join(data.frame(label = tree$tip.label,
                                species = tree$species),
                     by = "label")

}


#' Check validity of colors
#'
#' Checks whether the elements of a character vector are valid colors.
#'
#' @param x A character vector.
#'
#' @return A logical value indicating whether the input \code{x} contains valid colors.
#' @noRd
areColors <- function(x) {

  vapply(x, function(X) {

    tryCatch(is.matrix(grDevices::col2rgb(X)),
             error = function(e) FALSE)

  }, TRUE)

}


#' Wraps a string into several lines at "_" characters
#'
#' @param names Character vector containing the strings to wrap.
#'
#' @return Character vector after wrapping.
#' @noRd
wrapLongNames <- function(names) {

  sapply(names, function(name) {

    gsub(" ", "", paste0(strwrap(gsub("_", "_ ", name), width = 20), collapse = "\n"))

  })

}
