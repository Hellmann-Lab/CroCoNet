#' Plot the networks of modules
#'
#' Plots the top \code{N} strongest connections of the chosen module(s) based on the consensus network. The edge width indicates the consensus edge weights, while the edge color indicates the degree of edge divergence across species.
#'
#' @description In the CroCoNet approach, networks are reconstructed per clone and combined into a single phylogeny-aware consensus network which is the basis of the module assignment.
#'
#' This function selects the intramodular edges of the input module(s) in the consensus network, orders them by the consensus edge weight and keeps the top \code{N} edges per module for plotting. This means that only those module member genes will appear on the plot that are involved in the top \code{N} connections, the regulator itself might also be omitted if it is not particularly well-connected.
#'
#' The edges that were kept are plotted using \code{\link{ggraph}} per module. The default layout is "gem", other reasonable options include "stress", "circle", "nicely", "dh", "graphopt", "mds", "fr", "kk" and "drl" (see also \code{\link{ggraph}}, \code{\link{layout_tbl_graph_stress}} and \code{\link{layout_tbl_graph_igraph}}). The width of the edges represents the consensus edge weights, the range of widths can be set using \code{edge_width}.
#'
#' For each edge, an edge divergence score is calculated based on its edge weights in the networks of individual clones. The edge weights are compared across species using an ANOVA, and the F-statistic (i.e. between-species variability / within-species variability) is regarded as the measure of edge divergence and used to color the edges on the plot.
#'
#' If the aim is to plot the networks of several modules together, the input \code{module_names} should be a vector of module names. The \code{\link{ggraph}} plots are in this case combined together into a single \code{\link{patchwork}} object with the titles of the subplots matching the elements of \code{module_names}. All subplots have the same scale so that the edge weights and divergence scores are comparable across modules.
#'
#' @param module_names Character or character vector, the name(s) of the module(s) of interest.
#' @param pruned_modules Data frame of the pruned modules, required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' }
#' @param consensus_network \code{\link{igraph}} object, the consensus network across all species and clones.
#' @param network_list A named list of \code{\link{igraph}} objects containing the networks of all clones.
#' @param clone2species A data frame that specifies which species each clone belongs to, required columns:
#' \describe{
#' \item{clone}{Character, name of the clone.}
#' \item{species}{Character, name of the species.}
#' }
#' @param N Integer, the number of edges to plot (default: 200).
#' @param n_cores Integer, the number of cores (default: 1).
#' @param layout Character, graph layout (default: "kk", other reasonable options are "auto", "stress", "star", "circle", "nicely", "dh", "gem", "graphopt", "grid", "mds", "fr", "drl", "eigen", "fabric", "linear" and "unrooted"). See also \code{\link{ggraph}}, \code{\link{layout_tbl_graph_stress}} and \code{\link{layout_tbl_graph_igraph}}.
#' @param seed Integer, the seed for setting up the graph layout (default: 0, only relevant for certain layouts such as "gem", "nicely", "dh", "graphopt", "fr" and "drl").
#' @param colors Character vector, the colors to visualize the edge divergences. The vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}.
#' @param font_size Numeric, font size (default: 14).
#' @param edge_width Numeric vector of length 2, the range of edge widths for plotting the graph edges (default: c(0.2, 1.2)).
#' @param ncol Integer, the number of columns the subplots should be organized into if several modules are input. If NULL (default), the dimensions of the grid will follow the default of \code{\link{wrap_plots}}.
#'
#' @return A \code{\link{ggraph}} object in case \code{module_names} is a single module name and a \code{\link{patchwork}} object in case \code{module_names} is a vector of module names.
#' @export
#'
#' @examples plotNetworks("POU5F1", pruned_modules, consensus_network, network_list, clone2species)
plotNetworks <- function(module_names, pruned_modules, consensus_network, network_list, clone2species, N = 200L, n_cores = 1L, layout = "kk", seed = 0, colors = NULL, font_size = 14, edge_width = c(0.2, 1.2), ncol = NULL) {

  if (!layout %in% graph_layouts)
    stop("The chosen value of \"layout\" is not among the implemented options. Please choose one out of \"kk\", \"auto\", \"stress\", \"star\", \"circle\", \"nicely\", \"dh\", \"gem\", \"graphopt\", \"grid\", \"mds\", \"fr\", \"drl\", \"eigen\", \"fabric\", \"linear\" and \"unrooted\".")

  if (!inherits(seed, "numeric"))
    stop("The argument \"seed\" should be a numeric.")

  if (!is.null(colors) && (!inherits(colors, "character") || any(!areColors(colors))))
    stop("The argument \"colors\" should be a character vector of valid color representations.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  if (!inherits(edge_width, "numeric") || length(edge_width) != 2 || any(edge_width <= 0))
    stop("The argument \"edge_width\" should be a numeric vector of length 2 with positive values.")

  if (!is.null(ncol) && (length(ncol) != 1 || (!inherits(ncol, "integer") && !(inherits(ncol, "numeric") && ncol == round(ncol))) || ncol < 1))
    stop("The argument \"ncol\" should be a positive integer.")

  # get edges
  edges <- calculateEdgeDivergence(module_names, pruned_modules, consensus_network, network_list, clone2species, N)

  if (is.null(colors))
    colors <- edge_div_colors

  if (length(module_names) == 1) {

    return(
      plotNetwork(edges, layout, seed, colors, font_size, edge_width)
    )

  } else {

    lim_edge_weight = range(edges$consensus_weight)

    max_abs_div <- max(abs(log10(edges$f_statistic)))
    lim_edge_divergence <- c(-max_abs_div, max_abs_div)

    network_plots <- lapply(module_names, function(module) {

      plotNetwork(edges %>% dplyr::filter(.data[["regulator"]] == module), layout, seed, colors, font_size, edge_width, lim_edge_weight, lim_edge_divergence) +
        ggplot2::ggtitle(module)

    })

    network_plots[2:length(network_plots)] <- lapply(network_plots[2:length(network_plots)], function(plot) {plot + ggplot2::theme(legend.position = "none")})

    patchwork::wrap_plots(network_plots, guides = "collect", ncol = ncol)

  }


}


#' Calculate the divergence of intramodular edges
#'
#' Selects the top \code{N} strongest edges of the chosen module(s) based on the consensus network and calculates a divergence score for each edge.
#'
#' In the CroCoNet approach, networks are reconstructed per clone and combined into a single phylogeny-aware consensus network which is the basis of the module assignment.
#'
#' This function selects the intramodular edges of the input module(s) in the consensus network. If \code{N} is set to \code{Inf}, all intramodular edges are considered for the divergence calculation, if \code{N} is smaller than the module size, the edges are ordered by their consensus edge weight and only the top \code{N} edges are kept per module.
#'
#' For each edge that was kept, an edge divergence score is calculated based on its edge weights in the networks of individual clones. The edge weights are compared across species using an ANOVA. and the F-statistic (i.e. the variation across the species means / variation within the species) and the p-value of this F-statistic are output as the measures of edge divergence.
#'
#' @param module_names Character or character vector, the name(s) of the module(s) of interest.
#' @param pruned_modules Data frame of the pruned modules, required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' }
#' @param consensus_network \code{\link{igraph}} object, the consensus network across all species and clones.
#' @param network_list A list of \code{\link{igraph}} objects containing the networks per clone.
#' @param clone2species A data frame that specifies which species each clone belongs to, required columns:
#' \describe{
#' \item{clone}{Character, name of the clone.}
#' \item{species}{Character, name of the species.}
#' }
#' @param N Integer, the number of strongest edges to subset. If set to Inf (default), all edges in the module are considered.
#' @param n_cores Integer, the number of cores (default: 1).
#'
#' @return A data frame of the selected edges with 5 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{from, to}{Character, the 2 member genes of the regulator's module that form the edge.}
#' \item{consensus_weight}{Numeric, consensus edge weight/adjacency (the weighted average of clonewise adjacencies).}
#' \item{f_statistic}{Numeric, measue of edge divergence. It is calculated as the F-statistic from the ANOVA of edge weights with species as groups.}
#' \item{p-value}{Numeric, the p-value of the F-statistic.}
#' }
#' @export
#'
#' @examples POU5F1_mod_edge_divergence <- calculateEdgeDivergence("POU5F1",
#'                                                                  pruned_modules,
#'                                                                  consensus_network,
#'                                                                  network_list,
#'                                                                  clone2species)
calculateEdgeDivergence <- function(module_names, pruned_modules, consensus_network, network_list, clone2species, N = Inf, n_cores = 1L) {

  if (!inherits(module_names, "character"))
    stop("The argument \"module_names\" should be a character vector.")

  if (!is.data.frame(pruned_modules))
    stop("The argument \"pruned_modules\" should be a data frame.")

  if (any(!c("regulator", "target") %in% colnames(pruned_modules)))
    stop("The argument \"pruned_modules\" should contain the columns \"regulator\" and \"target\".")

  if (!inherits(consensus_network, "igraph"))
    stop("The argument \"consensus_network\" should be of class \"igraph\".")

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

  if (length(N) != 1 || (!inherits(N, "integer") && !(inherits(N, "numeric") && N == round(N))) || N < 1)
    stop("The argument \"N\" should be a positive integer.")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") && !(inherits(n_cores, "numeric") && n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  module = . = clone = consensus_weight = from = module = p.value = species = statistic = to = weight = NULL

  consensus_edges <- data.table::as.data.table(igraph::as_data_frame(consensus_network))
  clonewise_edges <- lapply(network_list, function(network) data.table::as.data.table(igraph::as_data_frame(network))) %>%
    data.table::rbindlist( idcol = "clone")

  doParallel::registerDoParallel(n_cores)

  edge_divergence <- foreach::foreach(module = module_names,
                                      .combine = dplyr::bind_rows,
                                      .multicombine = TRUE) %dopar% {

                                        # get module genes
                                        module_genes <- c(pruned_modules %>%
                                                            dplyr::filter(.data[["regulator"]] == module) %>%
                                                            dplyr::pull(.data[["target"]]),
                                                          module)

                                        module_genes_not_found <- module_genes[!module_genes %in% V(consensus_network)$name]

                                        if (length(module_genes_not_found) > 0)
                                          warning(paste0("The following genes in module ", module, " were not found in the provided consensus network: ", paste(module_genes_not_found, collapse = ", "), "."))

                                        # get consensus edges
                                        consensus_module_edges <- consensus_edges[
                                          from %in% module_genes & to %in% module_genes
                                        ][order(-weight)][1:min(N, .N), .(from, to, consensus_weight = weight)]

                                         consensus_module_edges <- data.table::as.data.table(clone2species)[
                                           , (consensus_module_edges), by = .(clone, species)]

                                        # get the adjacencies for all interactions between these genes in each clone
                                        clonewise_edges[
                                          consensus_module_edges, on = .(clone, from, to)][
                                            , .(species, clone, from, to, weight = data.table::fifelse(is.na(weight), 0, weight), consensus_weight)][
                                              , stats::oneway.test(weight ~ species, var.equal = TRUE)[c("statistic", "p.value")] %>% as.list(),
                                              by = .(from, to, consensus_weight)][
                                                order(from, to)][
                                                  , .(from, to, consensus_weight, f_statistic = statistic, p_value = p.value, regulator = module)] %>%
                                          as.data.frame()

                                      }

  doParallel::stopImplicitCluster()

  edge_divergence

}


#' Plot the network of a single module
#'
#' @param edges with 5 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{to}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' \item{consensus_weight}{Numeric, consensus edge weight/adjacency (the weighted average of clonewise adjacencies).}
#' \item{f_statistic}{Numeric, measue of edge divergence. It is calculated as the F-statistic from the ANOVA of edge weights with species as groups.}
#' \item{p-value}{Numeric, the p-value of the F-statistic.}
#' }
#' @param layout Character, graph layout (default: "gem", other reasonable options include "stress", "circle", "nicely", "dh", "graphopt", "mds", "fr", "kk" and "drl"). See also \code{\link{ggraph}}, \code{\link{layout_tbl_graph_stress}} and \code{\link{layout_tbl_graph_igraph}}.
#' @param seed Integer, the seed for setting up the graph layout (default: 0, only relevant for certain layouts such as "gem", "nicely", "dh", "graphopt", "fr" and "drl").
#' @param colors Character vector, the colors to visualize the distances. The vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}.
#' @param font_size Numeric, font size (default: 14).
#' @param edge_width Numeric vector of length 2, the range of edge widths for plotting the graph edges (default: c(0.2, 1.2)).
#' @param lim_edge_weight Numeric vector of length 2, the range of the edge weights to plot.
#' @param lim_edge_divergence Numeric vector of length 2, the range of the edge divergences to plot.
#'
#' @return A \code{\link{ggraph}} object.
#' @noRd
plotNetwork <- function(edges, layout = "kk", seed = 0, colors, font_size = 14, edge_width = c(0.2, 1.2), lim_edge_weight = NULL, lim_edge_divergence = NULL) {

  . = NULL

  # reconstruct graph
  graph <- igraph::graph_from_data_frame(edges)

  # set seed and schedule the restoration of the old seed
  old_seed <- globalenv()$.Random.seed
  on.exit(suspendInterrupts({
    if (is.null(old_seed)) {
      rm(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      assign(".Random.seed", value = old_seed, envir = globalenv(), inherits = FALSE)
    }
  }), add = TRUE)
  set.seed(seed)

  # layout
  gr <- ggraph::ggraph(graph, layout = layout)

  # mark regulator
  gr$data$is_regulator <- gr$data$name == unique(edges$regulator)

  # calculate min and max coordinates for plotting
  coord <- gr$data
  xmin <- min(coord$x)
  xmax <- max(coord$x)
  ymin <- min(coord$y)
  ymax <- max(coord$y)
  xmin_new <- xmin - (xmax - xmin)*0.1
  xmax_new <- xmax + (xmax - xmin)*0.1
  ymin_new <- ymin - (ymax - ymin)*0.05
  ymax_new <- ymax + (ymax - ymin)*0.05

  if (is.null(lim_edge_weight)) {

    lim_edge_weight = range(edges$consensus_weight)

  }

  if (is.null(lim_edge_divergence)) {

    max_abs_div <- max(abs(log10(edges$f_statistic)))
    lim_edge_divergence <- c(-max_abs_div, max_abs_div)

  }

  # plot network and color edges based on f statistic
  set.seed(seed)
  gr +
    ggraph::geom_edge_link(ggplot2::aes(edge_width = .data[["consensus_weight"]], color = log10(.data[["f_statistic"]]))) +
    ggraph::scale_edge_width(range = edge_width, limits = lim_edge_weight, guide = "none") +
    ggraph::scale_edge_color_gradientn(colors = colors,
                                       limits = lim_edge_divergence,
                                       breaks = lim_edge_divergence*0.9,
                                       labels = c("low", "high"),
                                       name = "edge\ndivergence") +
    ggraph::geom_node_point(size = 6, shape = 21, fill = "transparent", color = "transparent") +
    ggraph::geom_node_label(data = . %>% dplyr::filter(!.data[["is_regulator"]]), ggplot2::aes(label = .data[["name"]]), label.padding = ggplot2::unit(0.08, "lines"), label.size = 0.001, size = font_size / 5, fill = "grey90") +
    ggraph::geom_node_label(data = . %>% dplyr::filter(.data[["is_regulator"]]), ggplot2::aes(label = .data[["name"]]), label.padding = ggplot2::unit(0.11, "lines"), label.size = 0.005, size = font_size / 4, fill = "grey40", color = "white") +
    ggplot2::theme_void(base_size = font_size) +
    ggplot2::guides(edge_color = ggraph::guide_edge_colorbar()) +
    ggplot2::theme(legend.ticks = ggplot2::element_blank()) +
    ggplot2::xlim(xmin_new, xmax_new) +
    ggplot2::ylim(ymin_new, ymax_new)

}
