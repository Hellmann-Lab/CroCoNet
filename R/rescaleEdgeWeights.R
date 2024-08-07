#' Rescale edge weights between 0 and 1
#'
#' @param network_list A named list that contains the networks of each clone in an igraph format.
#' @param n_cores Number of cores.
#'
#' @return A named list that contains the networks of each clone in an igraph format, with the edge weights rescaled between 0 and 1.
#' @export
#'
#' @examples
#' network_list_scaled <- rescaleEdgeWeights(network_list_raw)
rescaleEdgeWeights <- function(network_list, n_cores = 1L) {

  if (!inherits(network_list, "list"))
    stop("The argument \"network_list\" should be a named list.")

  if (is.null(names(network_list)))
    stop("The argument \"network_list\" should be a named list.")

  if (any(!sapply(network_list, function(net) {inherits(net, "igraph")})))
    stop("All elements of \"network_list\" should be of class \"igraph\".")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") & !(inherits(n_cores, "numeric") & n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")


  weight = NULL

  # convert to data table
  dt_list <- convertToDT(network_list)

  # get all interaction scores from all networks
  all_weights <- data.table::rbindlist(dt_list)$weight

  # scale all values (if all positive: simply divide by the max, if it contains negative values: squish between 0 and 1)
  if (all(all_weights >= 0)) {

    dt_list_rescaled <- lapply(dt_list, function(dt) {

      dt[, weight := weight / max(all_weights)]

    })

  } else {

    warning("Negative edge weights found. Rescaling between 0 and 1 will remove the edge with the highest negative edge weight.")

    dt_list_rescaled <- lapply(dt_list, function(dt) {

      min_all_weights <- min(all_weights)
      max_all_weights <- max(all_weights)
      dt[, weight := (weight - min_all_weights) / (max_all_weights - min_all_weights)]

    })

  }

  # convert back to igraph
  convertToGraph(dt_list_rescaled, V(network_list[[1]])$name, n_cores)

}
