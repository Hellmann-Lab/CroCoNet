#' Plot the distribution of a preservation statistic
#'
#' Plots the distribution of the chosen preservation statistic(s) for the actual and random modules, for within-species and cross-species clone pairs. If a statistic is a good measure of module preservation, the preservation scores are expected to be 1) higher for the actual modules than for the random modules and 2) higher within species than across species.
#'
#' As part of the CroCoNet approach, pairwise module preservation scores are calculated between clones, both within and across species (see \code{\link{calculatePresStats}}) to gain information about the cross-species differences but also about the within-species diversity of the modules. These correlation-based preservation statistics quantify how well the module connectivity patterns are preserved between the networks of two clones. The statistics are also calculated not just for the actual, biologically meaningful modules, but also for random modules with matching sizes.
#'
#' The actual modules are expected to be more preserved than the random modules, and all modules, but especially the actual ones, are expected tp be more preserved within species than across species. The function plots the distributions of the within-species and cross-species scores both for the actual and for the random modules and thus allows these expectations to be visually checked. If several different statistics are input, these will be shown as the rows of the faceted plot. When it comes to choosing the best statistic, it is recommended to take the one that follow the most these expected trends.
#'
#' @param pres_stats Data frame of the preservation statistics for the actual (pruned) modules. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{clone1, clone2}{Character, the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species 'clone1' and 'clone2' belongs to, respectively.}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the preservation statistic of interest per module and clone pair.}
#' }
#' @param random_pres_stats Data frame of the preservation statistics for the random modules. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{clone1, clone2}{The names of the clones compared.}
#' \item{species1, species2}{The names of the species 'clone1' and 'clone2' belongs to, respectively.}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the preservation statistics of interest per module and clone pair.}
#' }
#' @param stats Character or character vector, the name(s) of the column(s) containing the statistics of interest.
#' @param colors Character vector of length 2, the colors for the within-species and cross-species preservation score distributions.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A \code{\link{ggplot}} object.
#' @export
#'
#' @examples plotPresStatDistributions(pres_stats, random_pres_stats, "cor_kIM")
plotPresStatDistributions <- function(pres_stats, random_pres_stats, stats, colors = NULL, font_size = 14) {

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

  if (!is.null(colors) & (!inherits(colors, "character") || any(!areColors(colors)) || length(colors) != 2))
    stop("The argument \"colors\" should be a character vector of valid color representations, with length of 2.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  if (is.null(colors))
    colors <- c("#FFA500", "#5A81A8")

  combined_pres_stats <- dplyr::bind_rows("actual modules" = pres_stats,
                                          "random modules" = random_pres_stats,
                                          .id = "module_set") %>%
    dplyr::mutate(category = factor(ifelse(.data[["species1"]] == .data[["species2"]], "within-species", "cross-species"), c("within-species", "cross-species"))) %>%
    tidyr::pivot_longer(cols = stats, names_to = "statistic", values_to = "value") %>%
    dplyr::mutate(statistic = factor(.data[["statistic"]], stats)) %>%
    dplyr::group_by(.data[["regulator"]], .data[["statistic"]]) %>%
    dplyr::filter(sum(is.na(.data[["value"]])) == 0) %>%
    dplyr::ungroup()

  combined_pres_stats %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["category"]], y = .data[["value"]], fill = .data[["category"]])) +
    ggplot2::geom_violin(draw_quantiles = 0.5) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::facet_grid(.data[["statistic"]] ~ .data[["module_set"]]) +
    ggplot2::xlab("type of comparison") +
    ggplot2::ylab("preservation score")

}


#' Plot cross-species VS within-species preservation statistics per species pair
#'
#' Plots the cross-species VS within-species preservation scores for all possible species pairs. If a statistic is a good measure of module preservation, the within-species scores are expected to be higher than the cross-species scores for all pairs, and the difference is expected to increase with increasing phylogenetic distance between the two species compared.
#'
#' As part of the CroCoNet approach, pairwise module preservation scores are calculated between clones, both within and across species (see \code{\link{calculatePresStats}}) to gain information about the cross-species differences but also about the within-species diversity of the modules. These correlation-based preservation statistics quantify how well the module connectivity patterns are preserved between the networks of two clones.
#'
#' For each possible speciesA-speciesB pair in the data, the function first subsets the within-species scores for speciesA, the within-species scores for speciesB and the cross-species scores between speciesA and speciesB, then calculates the mean of the within-species scores as well as the mean of the cross-species scores per module. For example, if there are 2 species, human and gorilla, with 3 human clones (\emph{h1}, \emph{h2} and \emph{h3}), and 2 gorilla clones (\emph{g1} and \emph{g2}), for a single module there will be 3 preservation scores within the human clones (\eqn{p_{h1-h2}}, \eqn{p_{h1-h3}} and \eqn{p_{h2-h3}}), 1 preservation score within the gorilla clones (\eqn{p_{g1-g2}}) and 6 preservation scores across the 2 species (\eqn{p_{h1-g1}}, \eqn{p_{h1-g2}}, \eqn{p_{h2-g1}}, \eqn{p_{h2-g2}}, \eqn{p_{h3-g1}} and \eqn{p_{h3-g2}}). For this module and the human-gorilla species pair, the summarized within-species score will then be the mean of the 4 within-species scores (\eqn{mean(p_{h1-h2}, p_{h1-h3}, p_{h2-h3}, p_{g1-g2})}) and the cross-species score will be the mean of the 6 cross-species scores (\eqn{mean(p_{h1-g1}, p_{h1-g2}, p_{h2-g1}, p_{h2-g2}, p_{h3-g1}, p_{h3-g2})}).
#'
#' After calculating the summarized within-species and cross-species scores for all modules and species pairs, the data points are shown faceted by species pair and colored by module size. The within-species scores are expected to be higher than the cross-species scores, thus the data points are expected to fall under the diagonal in all facets. However, the more phylogenetically distant two species are, the bigger deviation we expect from the diagonal. If several different statistics are input, these will be shown as the rows of the faceted plot. When it comes to choosing the best statistic, it is recommended to take the one that follow the most these expected trends.
#'
#' @param pres_stats Data frame of the preservation statistics for the actual (pruned) modules. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{clone1, clone2}{Character, the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species 'clone1' and 'clone2' belongs to, respectively.}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the preservation statistic of interest per module and clone pair.}
#' }
#' @param stats Character or character vector, the name(s) of the column(s) containing the statistics of interest.
#' @param colors Character vector, the colors to visualize the module sizes. The vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}.
#' @param font_size Numeric, font size (default: 14).
#' @param point_size Numeric, the size of the points (default: 0.3).
#' @param point_alpha Numeric, the opacity of the points (default: 0.7).
#' @return A \code{\link{ggplot}} object.
#' @export
#'
#' @examples plotPresStats(pres_stats, "cor_kIM")
plotPresStats <- function(pres_stats, stats, colors = NULL, font_size = 14, point_size = 0.3, point_alpha = 0.7) {

  if (!is.data.frame(pres_stats))
    stop("The argument \"pres_stats\" should be a data frame.")

  if (any(!c("regulator", "module_size", "clone1", "clone2", "species1", "species2") %in% colnames(pres_stats)))
    stop("The argument \"pres_stats\" should contain the columns \"regulator\", \"module_size\", \"clone1\", \"clone2\", \"species1\" and \"species2\".")

  if (!inherits(stats, "character"))
    stop("The argument \"stats\" should be a character or character vector.")

  missing_in_pres_stats <- stats[!stats %in% colnames(pres_stats)]

  if (length(missing_in_pres_stats) > 0)
    stop(paste0("The following statistic(s) in \"stats\" were not found in \"pres_stats\": ", paste(missing_in_pres_stats, collapse = ", "), "."))

  if (!is.null(colors) & (!inherits(colors, "character") || any(!areColors(colors))))
    stop("The argument \"colors\" should be a character vector of valid color representations.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  if (!inherits(point_size, "numeric") || length(point_size) != 1 || point_size <= 0)
    stop("The argument \"point_size\" should be a positive numeric value.")

  if (!inherits(point_alpha, "numeric") || length(point_alpha) != 1 || point_alpha < 0 || point_alpha > 1)
    stop("The argument \"point_alpha\" should be a numeric value between 0 and 1.")

  if (is.null(colors))
    colors <- module_size_colors

  combined_pres_stats <- pres_stats %>%
    dplyr::mutate(category = factor(ifelse(.data[["species1"]] == .data[["species2"]], "within_species", "cross_species"), c("within_species", "cross_species"))) %>%
    tidyr::pivot_longer(cols = stats, names_to = "statistic", values_to = "value") %>%
    dplyr::mutate(statistic = factor(.data[["statistic"]], stats)) %>%
    dplyr::group_by(.data[["regulator"]], .data[["statistic"]]) %>%
    dplyr::filter(sum(is.na(.data[["value"]])) == 0) %>%
    dplyr::ungroup() %>%
    groupIntoSpeciesPairs() %>%
    dplyr::group_by(.data[["regulator"]], .data[["module_size"]], .data[["statistic"]], .data[["species_compared"]], .data[["category"]]) %>%
    dplyr::summarise(value = mean(.data[["value"]])) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = "category", values_from = "value")

  min_score <- min(c(combined_pres_stats$within_species, combined_pres_stats$cross_species))
  max_score <- max(c(combined_pres_stats$within_species, combined_pres_stats$cross_species))

  ggplot2::ggplot(combined_pres_stats, ggplot2::aes(x = .data[["within_species"]], y = .data[["cross_species"]], colour = .data[["module_size"]])) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_color_gradientn(colours = colors, name = "module size") +
    ggplot2::xlim(min_score, max_score) +
    ggplot2::ylim(min_score, max_score) +
    ggplot2::xlab("within-species preservation score") +
    ggplot2::ylab("cross-species preservation score") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30") +
    ggplot2::facet_grid(.data[["statistic"]] ~ .data[["species_compared"]]) +
    ggplot2::theme_bw(base_size = font_size)

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
#' \item{species1, species2}{Character, the names of the species 'clone1' and 'clone2' belongs to, respectively.}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the preservation statistic of interest per module and clone pair.}
#' }
#'
#' @return Data frame of the preservation statistics reorganized per species pair.
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
