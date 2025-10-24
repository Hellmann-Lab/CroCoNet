#' Plot the tree representations of modules
#'
#' Plots the tree representations of one or more modules.
#'
#' The module trees are reconstructed based on pairwise preservation scores between replicates, both within and across species (see \code{\link{reconstructTrees}}). The tips of the resulting tree represent the replicates and the branch lengths represent the dissimilarity of module connectivity patterns between the networks of 2 replicates. Therefore, if 2 tips fall close to each other within a module tree, it means that the module topology is similar between the corresponding replicates.
#'
#' These trees carry information about the cross-species conservation of a module. If a module is diverged between 2 species (i.e. it is poorly preserved between the species but well-preserved within each species), we expect the tips of the tree to separate according to species and the 2 species to be connected by a long branch. In contrast, if a module is conserved between 2 species (i.e. it is similarly preserved both within and across species), we expect the tips of the different species to be well-mixed within the tree with no systematic separation.
#'
#' The module trees are plotted using \code{\link{ggtree}} as unrooted trees with 'daylight' layout and the tips are colored based on their species identity. If \code{show_labels} is set to TRUE, the tips will be in addition labelled as the corresponding replicate.
#'
#' If the aim is to plot trees of several modules together, the input should be a named list of module trees. The \code{\link{ggtree}} plots are in this case combined together into a single \code{\link{patchwork}} object with the titles of the subplots matching the names of the input list. All subplots have the same scale so that the branch lengths are comparable across modules.
#'
#' @param trees \code{\link{phylo}} object or a named list of \code{\link{phylo}} objects, the tree representations(s) of one or more modules.
#' @param show_labels Logical, if TRUE (default), the each tip will be labelled as the corresponding replicate, if FALSE, tips are not labelled.
#' @param species_colors Character vector, tip colors per species.
#' @param font_size Numeric, font size (default: 14).
#' @param tip_size Numeric, the size of the labels/circles at the tips (default: 1).
#' @param branch_width Numeric, the widths of the tree branches (default: 0.4).
#' @param ncol Integer, the number of columns the subplots should be organized into if several modules are input. If NULL (default), the dimensions of the grid will follow the default of \code{\link{wrap_plots}}.
#'
#' @return A \code{\link{ggtree}} object in case \code{trees} is a single \code{\link{phylo}} object and a \code{\link{patchwork}} object in case \code{trees} is a list of \code{\link{phylo}} objects.
#' @export
#'
#' @examples plotTrees(trees)
plotTrees <- function(trees, show_labels = TRUE, species_colors = NULL, font_size = 14, tip_size = 1, branch_width = 0.4, ncol = NULL) {

  if (inherits(trees, "list")) {

    if (any(!sapply(trees, function(tree) {inherits(tree, "phylo")})))
      stop("One or more list elements of \"trees\" are not data frames. \"trees\" should be either a data frame or a named list of data frames.")

    if (is.null(names(trees)))
      stop("The input object \"trees\" is a list without names. \"trees\" should be either a data frame or a named list of data frames.")

    for (id in names(trees)) {

      tree <- trees[[id]]

      if (is.null(tree$species))
        stop(paste0("No species information found for tree \"", id, "\". Please add a component \"species\" to all trees specifying which species each tip belongs to."))

      if (any(tree$edge_lengths < 0))
        warning(paste0("Negative branch lengths found in tree \"", id, "\"."))

    }

    species <- unique(trees[[1]]$species)

  } else if (inherits(trees, "phylo")) {

    if (is.null(trees$species))
      stop("No species information found for \"trees\". Please add a component \"species\" to the tree specifying which species each tip belongs to.")

    if (any(trees$edge_lengths < 0))
      warning("Negative branch lengths found in \"trees\".")

    species <- unique(trees$species)

  } else {

    stop("The argument \"trees\" should be an object of class \"phylo\" or a named list of objects of class \"phylo\".")

  }

  if (!inherits(show_labels, "logical") || length(show_labels) != 1)
    stop("The argument \"show_labels\" should be a logical value.")

  if (!is.null(species_colors) && (!inherits(species_colors, "character") || any(!areColors(species_colors))))
    stop("The argument \"species_colors\" should be a character vector of valid color representations.")

  if (!is.null(species_colors) && length(species_colors) != length(species))
    stop("The argument \"species_colors\" should contain as many colors as there are unique species in the trees.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  if (!inherits(tip_size, "numeric") || length(tip_size) != 1 || tip_size <= 0)
    stop("The argument \"tip_size\" should be a positive numeric value.")

  if (!inherits(branch_width, "numeric") || length(branch_width) != 1 || branch_width <= 0)
    stop("The argument \"branch_width\" should be a positive numeric value.")

  if (!is.null(ncol) && (length(ncol) != 1 || (!inherits(ncol, "integer") && !(inherits(ncol, "numeric") && ncol == round(ncol))) || ncol < 1))
    stop("The argument \"ncol\" should be a positive integer.")

  if (inherits(trees, "phylo")) {

    if (is.null(species_colors))
      species_colors <- species_color_ramp(seq(0, 1, len = length(unique(trees$species))))

    tree_df <- getTreeDf(trees)

    return(
      plotTree(trees, tree_df, species_colors, show_labels, tip_size, branch_width, font_size)
    )

  } else {

    if (is.null(species_colors))
      species_colors <- species_color_ramp(seq(0, 1, len = length(unique(trees[[1]]$species))))

    tree_plots <- lapply(names(trees), function(name) {

      tree <- trees[[name]]
      tree_df <- getTreeDf(tree)
      plotTree(tree, tree_df, species_colors, show_labels, tip_size, branch_width, font_size) +
        ggplot2::ggtitle(name)

    }) %>% matchScales()

    if (length(tree_plots) > 1)
      tree_plots[2:length(tree_plots)] <- lapply(tree_plots[2:length(tree_plots)], function(plot) {plot + ggplot2::theme(legend.position = "none")})

    patchwork::wrap_plots(tree_plots, guides = "collect", ncol = ncol)

  }

}


#' Plot the tree representation of a single module
#'
#' Creates a \code{\link{ggtree}} plot based on a \code{\link{phylo}} object and the corresponding data frame.
#' @param tree \code{\link{phylo}} object or a named list of \code{\link{phylo}} objects, the tree representations(s) of one or more modules.
#' @param tree_df The data frame format of 'tree' where each row corresponds to a branch.
#' @param species_colors Character vector, tip colors per species.
#' @param show_labels Logical, if TRUE (default), the each tip will be labelled as the corresponding replicate, if FALSE, tips are not labelled.
#' @param tip_size Numeric, the size of the labels/circles at the tips (default: 1).
#' @param branch_width Numeric, the widths of the tree branches (default: 0.4).
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A \code{\link{ggtree}} object.
#' @noRd
plotTree <- function(tree, tree_df, species_colors = NULL, show_labels = TRUE, tip_size = 1, branch_width = 0.4, font_size = 14) {

  if (is.null(species_colors))
    species_colors <- species_color_ramp(seq(0, 1, len = length(unique(tree$species))))

  if (show_labels) {

    light_species_colors <- colorspace::lighten(species_colors, 0.4)
    names(light_species_colors) <- names(species_colors)

    p <- suppressMessages(

      ggtree::ggtree(tree, layout = "daylight", size = branch_width) %<+%
        tree_df +
        ggtree::geom_tiplab(ggplot2::aes(label = .data[["label"]], fill = .data[["species"]]), angle = 0, hjust = 0.5,  geom = "label", color = "black", size = font_size/4, box.padding = ggplot2::unit(tip_size/5, "lines"), label.padding = ggplot2::unit(tip_size/5, "lines")) +
        ggplot2::scale_fill_manual(values = light_species_colors, breaks = levels(tree_df$species)) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(color = "transparent"))) +
        ggplot2::theme(text = ggplot2::element_text(size = font_size),
                       legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
                       plot.title = ggplot2::element_text(size = ggplot2::rel(1.2)))

    )

    return(p)

  } else {

    p <- suppressMessages(

      ggtree::ggtree(tree, layout = "daylight", size = branch_width) %<+%
        tree_df +
        ggtree::geom_tippoint(ggplot2::aes(fill = .data[["species"]]), shape = 21, size = tip_size*7, color = "transparent") +
        ggplot2::scale_fill_manual(values = species_colors, breaks = levels(tree_df$species)) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = tip_size*5))) +
        ggplot2::theme(text = ggplot2::element_text(size = font_size),
                       legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
                       plot.title = ggplot2::element_text(size = ggplot2::rel(1.2)))

    )

    p

  }

}


#' Match the scaling of plots
#'
#' Adjusts the scaling of each plot in \code{plot_list} so that the range of x and y are the same and the plots are centered.
#' @param plot_list A list of \code{\link{ggplot}} objects.
#'
#' @return A list of \code{\link{ggplot}} objects with matched scales.
#' @noRd
matchScales <- function(plot_list) {

  x_ranges <- lapply(plot_list, function(plot) {

    c(ggplot2::layer_scales(plot)$x$range$range[1], ggplot2::layer_scales(plot)$x$range$range[2])

  })

  max_x_range <- max(sapply(x_ranges, function(vec) {

    vec[2] - vec[1]

  }))*1.1

  x_ranges_new <- lapply(x_ranges, function(vec) {

    x_range <- vec[2] - vec[1]
    x_expand <- (max_x_range - x_range) / 2
    c(vec[1] - x_expand, vec[2] + x_expand)

  })

  y_ranges <- lapply(plot_list, function(plot) {

    c(ggplot2::layer_scales(plot)$y$range$range[1], ggplot2::layer_scales(plot)$y$range$range[2])

  })

  max_y_range <- max(sapply(y_ranges, function(vec) {

    vec[2] - vec[1]

  }))*1.1

  y_ranges_new <- lapply(y_ranges, function(vec) {

    y_range <- vec[2] - vec[1]
    y_expand <- (max_y_range - y_range) / 2
    c(vec[1] - y_expand, vec[2] + y_expand)

  })

  plot_list <- lapply(1:length(plot_list), function(i) {

    plot_list[[i]] +
      ggplot2::scale_x_continuous(limits = x_ranges_new[[i]]) +
      ggplot2::scale_y_continuous(limits = y_ranges_new[[i]])

  })

  plot_list

}
