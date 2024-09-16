#' Plot the distributions of tree-based statistics
#'
#' Plots the distribution(s) of the chosen tree-based statistic(s) for the actual and random modules. If the module trees represent well the differences of module topology within and across species, most tree-based statistics are expected to be lower for the actual modules than for the random modules.
#'
#' As part of the CroCoNet approach, pairwise module preservation scores are calculated between clones, both within and across species (see \code{\link{calculatePresStats}}) and neighbor-joining trees are reconstructed based on these preservation scores per module (see \code{\link{convertPresToDist}} and \code{\link{reconstructTrees}}). The tips of the resulting tree represent the clones and the branch lengths represent the dissimilarity of module topology between the networks of 2 clones. Various statistics can be defined based on these trees such as total tree length, total within-species diversity, diversity of a species and species-to-other branch length (see \code{\link{calculateTreeStats}}. These tree-based statistics can then be used to identify conserved/diverged modules (see \code{\link{findConservedDivergedModules}}) and pinpoint individual target genes within these modules that contribute the most to conservation/divergence (see \code{\link{findConservedDivergedTargets}}).
#'
#' Compared to the random modules, the actual modules are expected to be better preserved between clones and thus the branch lengths of their module trees are also expected to be shorter (especially between clones of the same species). As a result, most tree-based statistics (especially the ones measuring diversity) are expected to be lower for the actual modules than for the random modules. The function plots the distributions of the chosen tree-based statistics both for the actual and for the random modules and thus allows these expectations to be visually checked. If several different statistics are input, these will be shown as the rows of the faceted plot.
#'
#' @param tree_stats Data frame of the tree statistics for the actual (pruned) modules. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the tree statistic of interest per module and clone pair.}
#' }
#' @param random_tree_stats Data frame of the tree statistics for the random modules. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the tree statistics of interest per module and clone pair.}
#' }
#' @param stats Character or character vector, the name(s) of the column(s) containing the statistics of interest.
#' @param colors Character vector of length 2, the colors for the actual and random modules.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return  A violin plot as a \code{\link{ggplot}} object showing the distributions of the chosen tree-based statistics for both the actual and the random modules.
#' @export
#'
#' @examples plotTreeStatDistributions(tree_stats,
#'                                     random_tree_stats,
#'                                     c("within_species_diversity", "total_tree_length"))
#' @family functions to plot tree-based statistics
plotTreeStatDistributions <- function(tree_stats, random_tree_stats, stats, colors = NULL, font_size = 14) {

  # check input data
  if (!is.data.frame(tree_stats))
    stop("The argument \"tree_stats\" should be a data frame.")

  if (!"regulator" %in% colnames(tree_stats))
    stop("The argument \"tree_stats\" should contain the column \"regulator\".")

  if (!is.data.frame(random_tree_stats))
    stop("The argument \"random_tree_stats\" should be a data frame.")

  if (!"regulator" %in% colnames(random_tree_stats))
    stop("The argument \"random_tree_stats\" should contain the column \"regulator\".")

  if (!inherits(stats, "character"))
    stop("The argument \"stats\" should be a character or character vector.")

  missing_in_tree_stats <- stats[!stats %in% colnames(tree_stats)]

  if (length(missing_in_tree_stats) > 0)
    stop(paste0("The following statistic(s) in \"stats\" were not found in \"tree_stats\": ", paste(missing_in_tree_stats, collapse = ", "), "."))

  missing_in_random_tree_stats <- stats[!stats %in% colnames(random_tree_stats)]

  if (length(missing_in_random_tree_stats) > 0)
    stop(paste0("The following statistic(s) in \"stats\" were not found in \"random_tree_stats\": ", paste(missing_in_random_tree_stats, collapse = ", "), "."))

  if (!is.null(colors) && (!inherits(colors, "character") || any(!areColors(colors)) || length(colors) != 2))
    stop("The argument \"colors\" should be a character vector of valid color representations, with length of 2.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  # if no colors are provided, take the default
  if (is.null(colors))
    colors <- c("#5A81A8", "#FFA500")

  # combine the actual and random modules
  combined_tree_stats <- dplyr::bind_rows("actual\nmodules" = tree_stats,
                                          "random\nmodules" = random_tree_stats,
                                          .id = "module_set")

  # if only 1 statistic has to be plotted, do not facet
  if (length(stats) == 1) {

    # remove modules where any of the scores is NA
    combined_tree_stats <- combined_tree_stats %>%
      tidyr::drop_na(.data[[stats]])

    # initialize plot without facets
    p <- combined_tree_stats %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[["module_set"]], y = .data[[stats]], fill = .data[["module_set"]])) +
      ggplot2::ylab(stats)

  # if several statistics have to be plotted, facet bystatistic (y-axis)
  } else {

    # remove modules from the distribution of a statistic if any of their scores is NA (keep for the other statistics though)
    combined_tree_stats <- combined_tree_stats %>%
      tidyr::pivot_longer(cols = stats, names_to = "statistic", values_to = "value") %>%
      dplyr::mutate(statistic = factor(wrapLongNames(.data[["statistic"]]), wrapLongNames(stats))) %>%
      tidyr::drop_na(.data[["value"]])

    # initialize plot with y facet
    p <- combined_tree_stats %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[["module_set"]], y = .data[["value"]], fill = .data[["module_set"]])) +
      ggplot2::facet_grid(.data[["statistic"]] ~ .) +
      ggplot2::ylab("tree statistic")

  }

  # violin plot colored by actual and random modules
  p + ggplot2::geom_violin(draw_quantiles = 0.5) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = font_size, color = "black"),
                   axis.title.x = ggplot2::element_blank())

}


#' Plot tree-based statistics against each other
#'
#' Plots a chosen tree-based statistic (e.g. total tree length) against another (e.g. within-species diversity) for all pruned modules and if provided, all random modules.
#'
#' As part of the CroCoNet approach, pairwise module preservation scores are calculated between clones, both within and across species (see \code{\link{calculatePresStats}}) and neighbor-joining trees are reconstructed based on these preservation scores per module (see \code{\link{convertPresToDist}} and \code{\link{reconstructTrees}}). The tips of the resulting tree represent the clones and the branch lengths represent the dissimilarity of module topology between the networks of 2 clones.  Various statistics can be defined based on these trees such as total tree length, total within-species diversity, diversity of a species and species-to-other branch length (see \code{\link{calculateTreeStats}}).
#'
#' These tree-based statistics can be used to pinpoint conserved and diverged modules. First, a linear regression model is fit between the total tree length and within-species diversity (in case the focus of interest is conservation and overall divergence) or between the species-to-other branch length and the diversity of a species (in case the focus of interest is the species-specific divergence), and then outlier data points are identified that do not follow the general linear trend (see \code{\link{fitTreeStatsLm}} and \code{\link{findConservedDivergedModules}}).
#'
#' This function visualizes the relationship of 2 tree-based statistics. This can be an informative check for the pairs of statistics that are used to measure cross-species conservation - total-tree length VS within-species diversity (default) or species-to-other branch length VS the diversity of a species - before fitting a linear regression model between them.
#'
#' If \code{random_tree_stats} is provided, the actual and the random modules are plotted together, shown in 2 different colors. Most tree-based statistics (especially the ones measuring diversity) are expected to be lower for the actual modules than for the random modules, therefore the 2 sets of modules are expected to cluster separately on the plot with the actual modules located towards the lower values.
#'
#' If \code{random_tree_stats} is not provided, only the actual modules in \code{tree_stats} are shown on the plot. In this case, if the column "module_size" is present in \code{tree_stats}, the data points are colored by module size.
#'
#' @param tree_stats Data frame of the tree statistics for the actual (pruned) modules. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator (only needed if the data points are desired to be colored by module size).}
#' \item{\{\{nameOfStat\}\}}{Numeric, 2 columns containing the tree statistics to be plotted against each other.}
#' }
#' @param random_tree_stats Data frame of the tree statistics for the random modules (optional). If provided, the following columns are required:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{\{\{nameOfStat\}\}}{Numeric, 2 columns containing the tree statistics to be plotted against each other.}
#' }
#' @param stats Character vector of length 2, the names of the columns containing the 2 statistics that should be plotted against each other (default: c("within_species_diversity", "total_tree_length"), another useful option is to plot the species-to-other branch length against the diversity for a chosen species). The 1st element of the vector will be plotted on the x-axis and the 2nd element of the vector will be plotted on the y-axis.
#' @param colors Character vector, either the colors to visualize the module sizes (if only \code{tree_stats} is provided), or the colors to visualize the actual and random modules (if both \code{tree_stats} and \code{random_tree_stats} are provided). In the first case, the vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}. In the second case, the vector should contain 2 colors for the actual and random modules.
#' @param font_size Numeric, font size (default: 14).
#' @param point_size Numeric, the size of the points (default: 0.3).
#' @param point_alpha Numeric, the opacity of the points (default: 1).
#' @return A scatterplot as a \code{\link{ggplot}} object showing a chosen tree-based statistic (e.g. total tree length) against another (e.g. within-species diversity) for all pruned modules and if provided, all random modules.
#' @export
#'
#' @examples plotTreeStats(tree_stats,
#'                         random_tree_stats,
#'                         c("within_species_diversity", "total_tree_length"))
#' @family functions to plot tree-based statistics
plotTreeStats <- function(tree_stats, random_tree_stats = NULL, stats = c("within_species_diversity", "total_tree_length"), colors = NULL, font_size = 14, point_size = 0.05, point_alpha = 1) {

  # check input data
  if (!is.data.frame(tree_stats))
    stop("The argument \"tree_stats\" should be a data frame.")

  if (any(!c("regulator") %in% colnames(tree_stats)))
    stop("The argument \"tree_stats\" should contain the column \"regulator\".")

  if (!inherits(stats, "character") || length(stats) != 2)
    stop("The argument \"stats\" should be a character vector of length 2.")

  missing_in_tree_stats <- stats[!stats %in% colnames(tree_stats)]

  if (length(missing_in_tree_stats) > 0)
    stop(paste0("The following statistic(s) in \"stats\" were not found in \"tree_stats\": ", paste(missing_in_tree_stats, collapse = ", "), "."))

  if (!is.null(random_tree_stats)) {

    if (!is.data.frame(random_tree_stats))
      stop("The argument \"random_tree_stats\" should be a data frame.")

    if (any(!c("regulator") %in% colnames(random_tree_stats)))
      stop("The argument \"random_tree_stats\" should contain the column \"regulator\".")

    missing_in_random_tree_stats <- stats[!stats %in% colnames(random_tree_stats)]

    if (length(missing_in_random_tree_stats) > 0)
      stop(paste0("The following statistic(s) in \"stats\" were not found in \"random_tree_stats\": ", paste(missing_in_random_tree_stats, collapse = ", "), "."))

  }

  if (!is.null(colors) && (!inherits(colors, "character") || any(!areColors(colors))))
    stop("The argument \"colors\" should be a character vector of valid color representations.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  if (!inherits(point_size, "numeric") || length(point_size) != 1 || point_size <= 0)
    stop("The argument \"point_size\" should be a positive numeric value.")

  if (!inherits(point_alpha, "numeric") || length(point_alpha) != 1 || point_alpha < 0 || point_alpha > 1)
    stop("The argument \"point_alpha\" should be a numeric value between 0 and 1.")

  # if "random_tree_stats" is not provided, plot the actual modules only
  if (is.null(random_tree_stats)) {

    # remove hubs where any of the chosen statistics are NA
    combined_tree_stats <-  tree_stats %>%
      tidyr::drop_na(dplyr::any_of(stats))

    # if the info about module sizes is present, color the data points by module size
    if ("module_size" %in% colnames(combined_tree_stats)) {

      # if no custom color scale was provided, use the default module size colors
      if (is.null(colors))
        colors <- module_size_colors

      # setup plot
      p <- ggplot2::ggplot(combined_tree_stats, ggplot2::aes(x = .data[[stats[1]]], y = .data[[stats[2]]], colour = .data[["module_size"]])) +
        ggplot2::scale_color_gradientn(colours = colors, name = "module size")

    # if the info about module sizes is not present, keep the data points black
    } else {

      # setup plot
      p <- ggplot2::ggplot(combined_tree_stats, ggplot2::aes(x = .data[[stats[1]]], y = .data[[stats[2]]]))

    }

  # if "random_tree_stats" is provided, plot the actual and the random modules together, distinguished by color
  } else {

    # combine the actual and random modules & remove hubs where any of the chosen statistics are NA
    combined_tree_stats <-  dplyr::bind_rows("actual modules" = tree_stats,
                                             "random modules" = random_tree_stats,
                                             .id = "module_set") %>%
      tidyr::drop_na(dplyr::any_of(stats))

    # if no custom colors were provided, use the default colors
    if (is.null(colors))
      colors <- c("#5A81A8", "#FFA500")
    else if (length(colors) != 2)
      stop("If both \"tree_stats\" and \"random_tree_stats\" are provided, the argument \"colors\" should be a character vector of length of 2 specifying the colors for the actual and random modules.")

    # setup plot
    p <- ggplot2::ggplot(combined_tree_stats, ggplot2::aes(x = .data[[stats[1]]], y = .data[[stats[2]]], color = .data[["module_set"]])) +
      ggplot2::scale_color_manual(values = colors, guide = ggplot2::guide_legend(title = NULL, override.aes = list(size = 2)))

  }

  # if the CIs of the statistics are present, add error bars to the data points
  if (all(c(paste0("lwr_", stats[1]), paste0("upr_", stats[1]), paste0("lwr_", stats[2]), paste0("upr_", stats[2])) %in% colnames(combined_tree_stats)) &&
      all(!is.na(c(combined_tree_stats[[paste0("lwr_", stats[1])]], combined_tree_stats[[paste0("upr_", stats[1])]], combined_tree_stats[[paste0("lwr_", stats[2])]], combined_tree_stats[[paste0("upr_", stats[2])]])))) {

    p <- p +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data[[paste0("lwr_", stats[1])]], xmax = .data[[paste0("upr_", stats[1])]]), linewidth = point_size*4, alpha = point_alpha) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[[paste0("lwr_", stats[2])]], ymax = .data[[paste0("upr_", stats[2])]]), linewidth = point_size*4, alpha = point_alpha)

  }

  # scatterplot of the 1st statistic against the 2nd statistic
  p +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::xlab(stats[1]) +
    ggplot2::ylab(stats[2]) +
    ggplot2::theme_bw(base_size = font_size)

}
