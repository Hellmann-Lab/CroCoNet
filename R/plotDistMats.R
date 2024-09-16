#' Plot the distance matrices of clones based on module connectivity patterns
#'
#' Plots the pairwise distances between clones based on module connectivity patterns for one or more modules.
#'
#' As part of the CroCoNet approach, pairwise module preservation scores are calculated between clones, both within and across species (see \code{\link{calculatePresStats}}) to gain information about the cross-species differences but also about the within-species diversity of the modules. These correlation-based preservation statistics quantify how well the module connectivity patterns are preserved between the networks of two clones. They can be converted into distance measures using the formula \code{dist = (1 - pres)/2} (see \code{\link{convertPresToDist}}).
#'
#' This function plots the distance measures as a tile plot where each tile corresponds to a pair of clones and the color of the tile corresponds to the distance based on module connectivity patterns between this pair of clones. The distance of a clone with itself is always 0 and thus not meaningful, therefore these tiles are colored grey.
#'
#' If the aim is to plot distance matrices for several modules together, the input should be a named list of data frames, each containing the distances for one of the modules. The tile plots are in this case combined together into a single \code{\link{patchwork}} object with the titles of the subplots matching the names of the input list. All subplots have the same scale so that the distances are comparable across modules.
#'
#' @param dist_dfs Data frame or a named list of data frames containing the distance measures per clone pair for one or more modules. Required columns for each data frame:
#' \describe{
#' \item{clone1, clone2}{Character the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species \code{clone1} and \code{clone2} belong to, respectively.}
#' \item{dist}{Numeric, distance measure ranging from 0 to 1, calculated based on the preservation score of the given module between \code{clone1} and \code{clone2}.}
#' }
#' @param colors Character vector, the colors to visualize the distances. The vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}.
#' @param font_size Numeric, font size (default: 14).
#' @param ncol Integer, the number of columns the subplots should be organized into if several modules are input. If NULL (default), the dimensions of the grid will follow the default of \code{\link{wrap_plots}}.
#'
#' @return A \code{\link{ggplot}} object in case \code{dist_dfs} is a single data frame and a \code{\link{patchwork}} object in case \code{dist_dfs} is a list of data frames.
#' @export
#'
#' @examples plotDistMats(dist)
plotDistMats <- function(dist_dfs, colors = NULL, font_size = 14, ncol = NULL) {

  if (inherits(dist_dfs, "list")) {

    if (any(!sapply(dist_dfs, is.data.frame)))
      stop("One or more list elements of \"dist_dfs\" are not data frames. \"dist_dfs\" should be either a data frame or a named list of data frames.")

    if (is.null(names(dist_dfs)))
      stop("The input object \"dist_dfs\" is a list without names. \"dist_dfs\" should be either a data frame or a named list of data frames.")

    if (!all(sapply(dist_dfs, function(df) {

      all(c("clone1", "clone2", "species1", "species2", "dist") %in% colnames(df))

    })))
      stop("All data frames in \"dist_dfs\" are expected to contain the columns \"clone1\", \"clone2\", \"species1\", \"species2\" and \"dist\".")

  } else if (is.data.frame(dist_dfs)) {

    if (any(!c("clone1", "clone2", "species1", "species2", "dist") %in% colnames(dist_dfs)))
      stop("The argument \"dist_dfs\" should contain the columns \"clone1\", \"clone2\", \"species1\", \"species2\" and \"dist\".")

  } else {

    stop("The argument \"dist_dfs\" should be a data frame or a named list of data frames.")

  }

  if (!is.null(colors) && (!inherits(colors, "character") || any(!areColors(colors))))
    stop("The argument \"colors\" should be a character vector of valid color representations.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  if (!is.null(ncol) && (length(ncol) != 1 || (!inherits(ncol, "integer") && !(inherits(ncol, "numeric") && ncol == round(ncol))) || ncol < 1))
    stop("The argument \"ncol\" should be a positive integer.")

  if (is.null(colors))
    colors <- dist_colors

  if (!inherits(dist_dfs, "list")) {

    return(
      plotDistMat(dist_dfs, colors, font_size)
    )

  } else {

    dist_range <- dplyr::bind_rows(dist_dfs) %>%
      dplyr::pull(.data[["dist"]]) %>%
      range()

    dist_plots <- lapply(names(dist_dfs), function(name) {

      suppressMessages(
        plotDistMat(dist_dfs[[name]], colors, font_size) +
        ggplot2::scale_fill_gradientn(colors = colors, na.value = "grey95", limits = dist_range, name = "distance") +
        ggplot2::ggtitle(name)
      )

    })

    if (length(dist_plots) > 1)
      dist_plots[2:length(dist_plots)] <- lapply(dist_plots[2:length(dist_plots)], function(plot) {plot + ggplot2::theme(legend.position = "none")})

    patchwork::wrap_plots(dist_plots, guides = "collect", ncol = ncol)

  }

}

#' Plot the distance matrix of clones based on the connectivity patterns of a single module
#'
#' Plots the pairwise distances between clones based on module connectivity patterns as a tile plot.
#'
#' @param dist_df Data frame containing the distance measures per clone pair, required columns:
#' \describe{
#' \item{clone1, clone2}{Character the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species \code{clone1} and \code{clone2} belong to, respectively.}
#' \item{dist}{Numeric, distance measure ranging from 0 to 1, calculated based on the preservation score of the given module between \code{clone1} and \code{clone2}.}
#' }.
#' @param colors Character vector, the colors to visualize the distances. The vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A \code{\link{ggplot}} object
#' @noRd
plotDistMat <- function(dist_df, colors, font_size = 14) {

  dist_df <- dist_df %>%
    dplyr::bind_rows(cbind(dist_df %>%
                             dplyr::distinct(.data[["clone1"]], .data[["species1"]]),
                           dist_df %>%
                             dplyr::distinct(.data[["clone1"]], .data[["species1"]]) %>%
                             dplyr::rename(clone2 = .data[["clone1"]], species2 = .data[["species1"]])))

  ggplot2::ggplot(dist_df, ggplot2::aes(x = .data[["clone1"]], y = .data[["clone2"]], fill = .data[["dist"]])) +
    ggplot2::geom_tile(colour = "white", size = 0.5) +
    ggplot2::scale_y_discrete(expand = c(0,0), limits = rev, position = "right") +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_gradientn(colors = colors, na.value = "grey95", name = "distance") +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::facet_grid(.data[["species2"]] ~ .data[["species1"]], scales = "free", space = "free", switch = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

}
