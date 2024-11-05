#' Plot the distributions of preservation statistics
#'
#' Plots the distribution(s) of the chosen preservation statistic(s) for the actual and random modules and for within-species and cross-species clone pairs. If a statistic is a good measure of module preservation, the preservation scores are expected to be 1) higher for the actual modules than for the random modules and 2) higher within species than across species.
#'
#' As part of the CroCoNet approach, pairwise module preservation scores are calculated between clones, both within and across species (see \code{\link{calculatePresStats}}) to gain information about the cross-species differences but also about the within-species diversity of the modules. These correlation-based preservation statistics quantify how well the module topology is preserved between the networks of two clones. The statistics are calculated not just for the actual, biologically meaningful modules, but also for random modules with matching sizes.
#'
#' The actual modules are expected to be more preserved than the random modules, and all modules, but especially the actual ones, are expected to be more preserved within species than across species. The function plots the distributions of the within-species and cross-species scores both for the actual and for the random modules and thus allows these expectations to be visually checked. If several different statistics are input, these will be shown as the rows of the faceted plot. When it comes to choosing the best statistic, it is recommended to take the one that follows these expected trends the most.
#'
#' @param pres_stats Data frame of the preservation statistics for the actual (pruned) modules. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{clone1, clone2}{Character, the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species \code{clone1} and \code{clone2} belongs to, respectively.}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the preservation statistic(s) of interest per module and clone pair.}
#' }
#' @param random_pres_stats Data frame of the preservation statistics for the random modules. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{clone1, clone2}{The names of the clones compared.}
#' \item{species1, species2}{The names of the species \code{clone1} and \code{clone2} belongs to, respectively.}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the preservation statistic(s) of interest per module and clone pair.}
#' }
#' @param stats Character or character vector, the name(s) of the column(s) containing the statistic(s) of interest.
#' @param colors Character vector of length 2, the colors for the actual and random modules.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A violin plot as a \code{\link{ggplot}} object showing the distributions of chosen preservation statistics both within species and across species and both for the actual and for the random modules.
#' @export
#'
#' @examples plotPresStatDistributions(pres_stats, random_pres_stats, "cor_kIM")
#' @family functions to plot preservation statistics
plotPresStatDistributions <- function(pres_stats, random_pres_stats, stats, colors = NULL, font_size = 14) {

  # check input data
  if (!is.data.frame(pres_stats))
    stop("The argument \"pres_stats\" should be a data frame.")

  if (any(!c("regulator", "clone1", "clone2", "species1", "species2") %in% colnames(pres_stats)))
    stop("The argument \"pres_stats\" should contain the columns \"regulator\", \"clone1\", \"clone2\", \"species1\" and \"species2\".")

  if (!is.data.frame(random_pres_stats))
    stop("The argument \"random_pres_stats\" should be a data frame.")

  if (any(!c("regulator", "clone1", "clone2", "species1", "species2") %in% colnames(random_pres_stats)))
    stop("The argument \"random_pres_stats\" should contain the columns \"regulator\", \"clone1\", \"clone2\", \"species1\" and \"species2\".")

  if (!inherits(stats, "character"))
    stop("The argument \"stats\" should be a character or character vector.")

  missing_in_pres_stats <- stats[!stats %in% colnames(pres_stats)]

  if (length(missing_in_pres_stats) > 0)
    stop(paste0("The following statistic(s) in \"stats\" were not found in \"pres_stats\": ", paste(missing_in_pres_stats, collapse = ", "), "."))

  missing_in_random_pres_stats <- stats[!stats %in% colnames(random_pres_stats)]

  if (length(missing_in_random_pres_stats) > 0)
    stop(paste0("The following statistic(s) in \"stats\" were not found in \"random_pres_stats\": ", paste(missing_in_random_pres_stats, collapse = ", "), "."))

  if (!is.null(colors) && (!inherits(colors, "character") || any(!areColors(colors)) || length(colors) != 2))
    stop("The argument \"colors\" should be a character vector of valid color representations, with length of 2.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  # if no colors are provided, take the default
  if (is.null(colors))
    colors <- c("#5A81A8", "#FFA500")

  # combine the actual and random modules
  combined_pres_stats <- dplyr::bind_rows("actual\nmodules" = pres_stats,
                                          "random\nmodules" = random_pres_stats,
                                          .id = "module_set") %>%
    # label each score as within-species or cross-species
    dplyr::mutate(category = factor(ifelse(.data[["species1"]] == .data[["species2"]], "within-species", "cross-species"), c("within-species", "cross-species")))

  # if only 1 statistic has to be plotted, facet by within-species/cross-species scores only (x-axis)
  if (length(stats) == 1) {

    # remove modules where any of the scores is NA
    combined_pres_stats <- combined_pres_stats %>%
      dplyr::group_by(dplyr::across(c("module_set", "regulator"))) %>%
      dplyr::filter(sum(is.na(.data[[stats]])) == 0) %>%
      dplyr::ungroup()

    # initialize plot with x facets only
    p <- combined_pres_stats %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[["module_set"]], y = .data[[stats]], fill = .data[["module_set"]])) +
      ggplot2::facet_wrap( ~ .data[["category"]]) +
      ggplot2::ylab(stats)

  # if several statistics have to be plotted, facet by within-species/cross-species scores (x-axis) and statistic (y-axis)
  } else {

    # remove modules from the distribution of a statistic if any of their scores is NA (keep for the other statistics though)
    combined_pres_stats <- combined_pres_stats %>%
      tidyr::pivot_longer(cols = stats, names_to = "statistic", values_to = "value") %>%
      dplyr::mutate(statistic = factor(.data[["statistic"]], stats)) %>%
      dplyr::group_by(dplyr::across(c("module_set", "regulator", "statistic"))) %>%
      dplyr::filter(sum(is.na(.data[["value"]])) == 0) %>%
      dplyr::ungroup()

    # initialize plot with x and y facets
    p <- combined_pres_stats %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[["module_set"]], y = .data[["value"]], fill = .data[["module_set"]])) +
      ggplot2::facet_grid(.data[["statistic"]] ~ .data[["category"]]) +
      ggplot2::ylab("preservation score")

  }

  # violin plot colored by actual and random modules
  p +
    ggplot2::geom_violin(draw_quantiles = 0.5) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = font_size - 2, color = "black"))

}


#' Plot cross-species VS within-species preservation statistics per species pair
#'
#' Plots the cross-species VS within-species preservation scores for all possible species pairs. If a statistic is a good measure of module preservation, the within-species scores are expected to be higher than the cross-species scores for all species pairs, and the difference is expected to increase with increasing phylogenetic distance between the two species compared.
#'
#' As part of the CroCoNet approach, pairwise module preservation scores are calculated between clones, both within and across species (see \code{\link{calculatePresStats}}) to gain information about the cross-species differences but also about the within-species diversity of the modules. These correlation-based preservation statistics quantify how well the module topology is preserved between the networks of two clones.
#'
#' For each possible speciesA-speciesB pair in the data, the function first subsets the within-species scores for speciesA, the within-species scores for speciesB and the cross-species scores between speciesA and speciesB, then calculates the mean of the within-species scores as well as the mean of the cross-species scores per module.
#'
#' For example, if there are 2 species, human and gorilla, with 3 human clones (\emph{h1}, \emph{h2} and \emph{h3}), and 2 gorilla clones (\emph{g1} and \emph{g2}), for a single module there will be 3 preservation scores within the human clones (\eqn{p_{h1-h2}}, \eqn{p_{h1-h3}} and \eqn{p_{h2-h3}}), 1 preservation score within the gorilla clones (\eqn{p_{g1-g2}}) and 6 preservation scores across the 2 species (\eqn{p_{h1-g1}}, \eqn{p_{h1-g2}}, \eqn{p_{h2-g1}}, \eqn{p_{h2-g2}}, \eqn{p_{h3-g1}} and \eqn{p_{h3-g2}}). For this module and the human-gorilla species pair, the summarized within-species score will then be the mean of the 4 within-species scores (\eqn{mean(p_{h1-h2}, p_{h1-h3}, p_{h2-h3}, p_{g1-g2})}) and the cross-species score will be the mean of the 6 cross-species scores (\eqn{mean(p_{h1-g1}, p_{h1-g2}, p_{h2-g1}, p_{h2-g2}, p_{h3-g1}, p_{h3-g2})}).
#'
#' After calculating the summarized scores for all modules and species pairs, the cross-species scores are plotted against the within-species scores faceted by species pair. The within-species scores are expected to be higher than the cross-species scores, thus the data points are expected to fall under the diagonal in all facets. However, the more phylogenetically distant two species are, the bigger deviation we expect from the diagonal.
#'
#' If \code{random_pres_stats} is provided, the actual and the random modules are plotted together, shown in 2 different colors. The actual modules are expected to better preserved than the random modules (especially within species), therefore the 2 sets of modules are expected to cluster separately on the plot, with the actual modules located towards the higher (within-species) scores.
#'
#' If \code{random_pres_stats} is not provided, only the actual modules in \code{pres_stats} are shown on the plot. In this case, if the column "module_size" is present in \code{pres_stats}, the data points are colored by module size.
#'
#' If several different statistics are input, these will be shown as the rows of the faceted plot. When it comes to choosing the best statistic, it is recommended to take the one that follows the most the expected trends in terms of within-species VS cross-species scores, phylogeny, and actual VS random modules.
#'
#' @param pres_stats Data frame of the preservation statistics for the actual (pruned) modules. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator (only needed if the data points are desired to be colored by module size).}
#' \item{clone1, clone2}{Character, the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species \code{clone1} and \code{clone2} belongs to, respectively.}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the preservation statistic of interest per module and clone pair.}
#' }
#' @param random_pres_stats Data frame of the preservation statistics for the random modules (optional). If provided, the following columns are required:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{clone1, clone2}{The names of the clones compared.}
#' \item{species1, species2}{The names of the species \code{clone1} and \code{clone2} belongs to, respectively.}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the preservation statistic(s) of interest per module and clone pair.}
#' }
#' @param stats Character or character vector, the name(s) of the column(s) containing the statistics of interest.
#' @param colors Character vector, either the colors to visualize the module sizes (if only \code{pres_stats} is provided), or the colors to visualize the actual and random modules (if both \code{pres_stats} and \code{random_pres_stats} are provided). In the first case, the vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}. In the second case, the vector should contain 2 colors for the actual and random modules.
#' @param font_size Numeric, font size (default: 14).
#' @param point_size Numeric, the size of the points (default: 0.3).
#' @param point_alpha Numeric, the opacity of the points (default: 0.7).
#' @return A scatterplot as a \code{\link{ggplot}} object showing the cross-species VS within-species preservation scores per species pair for all pruned modules and if provided, all random modules.
#' @export
#'
#' @examples plotPresStats(pres_stats, stats = "cor_kIM")
#' @family functions to plot preservation statistics
plotPresStats <- function(pres_stats, random_pres_stats = NULL, stats, colors = NULL, font_size = 14, point_size = 0.1, point_alpha = 0.7) {

  # check input data
  if (!is.data.frame(pres_stats))
    stop("The argument \"pres_stats\" should be a data frame.")

  if (any(!c("regulator", "clone1", "clone2", "species1", "species2") %in% colnames(pres_stats)))
    stop("The argument \"pres_stats\" should contain the columns \"regulator\", \"clone1\", \"clone2\", \"species1\" and \"species2\".")

  if (!inherits(stats, "character"))
    stop("The argument \"stats\" should be a character or character vector.")

  missing_in_pres_stats <- stats[!stats %in% colnames(pres_stats)]

  if (length(missing_in_pres_stats) > 0)
    stop(paste0("The following statistic(s) in \"stats\" were not found in \"pres_stats\": ", paste(missing_in_pres_stats, collapse = ", "), "."))

  if (!is.null(colors) && (!inherits(colors, "character") || any(!areColors(colors))))
    stop("The argument \"colors\" should be a character vector of valid color representations.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  if (!inherits(point_size, "numeric") || length(point_size) != 1 || point_size <= 0)
    stop("The argument \"point_size\" should be a positive numeric value.")

  if (!inherits(point_alpha, "numeric") || length(point_alpha) != 1 || point_alpha < 0 || point_alpha > 1)
    stop("The argument \"point_alpha\" should be a numeric value between 0 and 1.")

  ## Data preparation & filtering

  # if "random_pres_stats" is not provided, take the actual modules only
  if (is.null(random_pres_stats)) {

    combined_pres_stats <-  pres_stats

  # if "random_pres_stats" is provided, combine the actual and random modules
  } else {

    combined_pres_stats <-  dplyr::bind_rows("actual modules" = pres_stats,
                                             "random modules" = random_pres_stats,
                                             .id = "module_set")

  }

  # if only 1 statistic has to be plotted, remove modules where any of the scores for this 1 statistic is NA
  if (length(stats) == 1) {

    combined_pres_stats <- combined_pres_stats %>%
      dplyr::rename(value = .data[[stats]]) %>%
      dplyr::group_by(dplyr::across(dplyr::any_of(c("module_set", "regulator")))) %>%
      dplyr::filter(sum(is.na(.data[["value"]])) == 0) %>%
      dplyr::ungroup()

    # if several statistics have to be plotted, remove modules from the distribution of a statistic if any of their scores is NA (keep for the other statistics though)
  } else {

    combined_pres_stats <- combined_pres_stats %>%
      tidyr::pivot_longer(cols = stats, names_to = "statistic", values_to = "value") %>%
      dplyr::mutate(statistic = factor(.data[["statistic"]], stats)) %>%
      dplyr::group_by(dplyr::across(dplyr::any_of(c("module_set", "regulator", "statistic")))) %>%
      dplyr::filter(sum(is.na(.data[["value"]])) == 0) %>%
      dplyr::ungroup()

  }

  ## Summarization of scores

  # divide into categories (within/cross-species) and species pairs, then summarize
  combined_pres_stats <- combined_pres_stats %>%
    dplyr::mutate(category = factor(ifelse(.data[["species1"]] == .data[["species2"]], "within_species", "cross_species"), c("within_species", "cross_species"))) %>%
    groupIntoSpeciesPairs() %>%
    dplyr::group_by(dplyr::across(dplyr::any_of(c("module_set", "regulator", "module_size", "statistic", "species_compared", "category")))) %>%
    dplyr::summarise(value = mean(.data[["value"]])) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = "category", values_from = "value")

  # minimum and maximum score
  min_score <- min(c(combined_pres_stats$within_species, combined_pres_stats$cross_species))
  max_score <- max(c(combined_pres_stats$within_species, combined_pres_stats$cross_species))

  ## Basic plot

  # scatterplot of cross-species VS within-species scores
  p <- ggplot2::ggplot(combined_pres_stats, ggplot2::aes(x = .data[["within_species"]], y = .data[["cross_species"]])) +
    ggplot2::xlim(min_score, max_score) +
    ggplot2::ylim(min_score, max_score) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30") +
    ggplot2::theme_bw(base_size = font_size)

  ## Color

  # if "random_pres_stats" is not provided, color the data points either by module size or nothing
  if (is.null(random_pres_stats)) {

    # if the info about module sizes is present, color by module size
    if ("module_size" %in% colnames(combined_pres_stats)) {

      # if no custom color scale was provided, use the default module size colors
      if (is.null(colors))
        colors <- module_size_colors

      p <- p +
        ggplot2::geom_point(ggplot2::aes(colour = .data[["module_size"]]), size = point_size, alpha = point_alpha) +
        ggplot2::scale_color_gradientn(colours = colors, name = "module size")

    # if the info about module sizes is not present, keep the data points black
    } else {

      p <- p +
        ggplot2::geom_point(size = point_size, alpha = point_alpha)

    }

   # if "random_pres_stats" is provided, color the data points by the module set (actual/random)
  } else {

    # if no custom color scale was provided, use the default colors
    if (is.null(colors))
      colors <- c("#5A81A8", "#FFA500")

    p <- p +
      ggplot2::geom_point(ggplot2::aes(colour = .data[["module_set"]]), size = point_size, alpha = point_alpha) +
      ggplot2::scale_color_manual(values = colors, guide = ggplot2::guide_legend(title = NULL, override.aes = list(size = 2)))

  }

  ## Faceting & axis labels

  # if only 1 statistic has to be plotted, facet by species pair only (x-axis) and use the name of the statistic as axis labels
  if (length(stats) == 1) {

   p <- p +
     ggplot2::facet_wrap(~ .data[["species_compared"]]) +
     ggplot2::xlab(parse(text = paste0(stats, '["within-species"]'))) +
     ggplot2::ylab(parse(text = paste0(stats, '["cross-species"]')))

  # if several statistics have to be plotted, facet by species pair (x-axis) and statistic (y-axis) and use a generic axis label
  } else {

    p <-  p +
      ggplot2::facet_grid(.data[["statistic"]] ~ .data[["species_compared"]]) +
      ggplot2::xlab("within-species preservation score") +
      ggplot2::ylab("cross-species preservation score")

  }

  p

}

#' Group preservation statistics into species pairs
#'
#' For each possible speciesA-speciesB pair in the data, the function subsets the within-species scores for speciesA, the within-species scores for speciesB and the cross-species scores between speciesA and speciesB, then combines these subsets across all pairs.
#'
#' The output allows the differences between within-species and cross-species scores to be compared between species pairs of varying evolutionay distances.
#'
#' In the final data frame, each within-species score appears several times, once for each pair the given species is part of, while each cross-species score appears only once at the corresponding species pair.
#'
#' @param pres_stats Data frame of the preservation statistics. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{clone1, clone2}{Character, the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species \code{clone1} and \code{clone2} belongs to, respectively.}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the preservation statistic of interest per module and clone pair.}
#' }
#'
#' @return Data frame of the preservation statistics reorganized per species pair. In addition to the columns of \code{pres_stats}, it contains 1 new column:
#' \describe{
#' \item{species_compared}{Character, the species comparison for which the preservation score can be used (format: "speciesA VS speciesB"). For each value of this column, the data frame contains the within-species scores for speciesA, the within-species scores for speciesB and the cross-species scores between speciesA and speciesB. }
#' }
#' @noRd
groupIntoSpeciesPairs <- function(pres_stats) {

  species <- unique(pres_stats$species1)
  n_species <- length(species)

  pres_stats_grouped <- NULL

  for (i in 1:(n_species-1)) {

    for (j in (i+1):n_species) {

      speciesA <- species[i]
      speciesB <- species[j]

      pres_stats_grouped <- dplyr::bind_rows(pres_stats_grouped,
                                             pres_stats %>%
                                               dplyr::filter((.data[["species1"]] == speciesA & .data[["species2"]] == speciesA) |
                                                               (.data[["species1"]] == speciesB & .data[["species2"]] == speciesB) |
                                                               (.data[["species1"]] == speciesA & .data[["species2"]] == speciesB) |
                                                               (.data[["species1"]] == speciesB & .data[["species2"]] == speciesA)) %>%
                                               dplyr::mutate(species_compared = paste0(speciesA, " VS ", speciesB)))

    }

  }

  pres_stats_grouped %>%
    dplyr::mutate(species_compared = factor(.data[["species_compared"]], unique(.data[["species_compared"]])))

}
