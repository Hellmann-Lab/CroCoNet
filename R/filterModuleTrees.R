#' Filter modules based on their tree representations
#'
#' Filters modules using two tree-based statistics (the total tree length and within-species diversity) by comparing the values of the actual modules to the values of the random modules.
#'
#' The within-species diversity measures the variability of module connectivity patterns within species and the total tree length measures the variability of module connectivity patterns both within and across species. The within-species diversity is a subset of the total tree, hence we expect a linear relationship between the two statistics. This line explains the detection robustness: modules that have both low within-species diversity and low total tree length are well-preserved both within and across species, meaning that these modules could be robustly detected in all replicates, whereas modules, for which both metrics are high, are poorly preserved not just across but also within species, indicating a high detection uncertainty ("wobbliness"). Random modules are expected to fall at the wobbly end of the spectrum, while actual modules are expected to fall at the robust end of the spectrum.
#'
#' The function removes the modules that are uncertain to detect, i.e. modules the are too similar to the random modules in terms of their total tree length and within-species diversity. It calculates how probable it is that the statistics of an actual module come from the distribution of all actual modules (\eqn{p_{actual}}) and how probable it is that they come from the distribution of all random modules (\eqn{p_{random}}) using the probability density functions of the two bivariate normal distributions. If \deqn{\frac{p_{actual}}{p_{actual} + p_{random}} > p_{cutoff}} is fulfilled, the module is kept, otherwise it is removed.
#' @param tree_stats Data frame of the tree-based statistics for the actual (pruned) modules. Required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{total_tree_length}{Numeric, total tree length per module (typically the median across all jackknife versions of a module.}
#' \item{within_species_diversity}{Numeric, within-species diveristy per module (typically the median across all jackknife versions of a module.}
#' }
#' @param random_tree_stats Data frame of the tree-based statistics for the random modules. Required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{total_tree_length}{Numeric, total tree length per module (typically the median across all jackknife versions of a module).}
#' \item{within_species_diversity}{Numeric, within-species diveristy per module (typically the median across all jackknife versions of a module).}
#' }
#' @param p_cutoff Numeric, the cutoff for the probability of a module's statistics being drawn from the distribution of the actual modules and not the distribution of random modules. Modules are retained only if this probability exceeds the cutoff.
#'
#' @return Data frame of the tree-based statistics for the actual (pruned) modules after filtering.
#' @export
#' @export
filterModuleTrees <- function(tree_stats, random_tree_stats, p_cutoff = 0.95) {

  if (!is.data.frame(tree_stats))
    stop("The argument \"tree_stats\" should be a data frame.")

  if (!is.data.frame(random_tree_stats))
    stop("The argument \"random_tree_stats\" should be a data frame.")

  if (any(!(c("regulator", "total_tree_length", "within_species_diversity") %in% colnames(tree_stats))))
    stop("The argument \"tree_stats\" should contain the columns \"regulator\", \"total_tree_length\" and \"within_species_diversity\".")

  if (any(!(c("regulator", "total_tree_length", "within_species_diversity") %in% colnames(random_tree_stats))))
    stop("The argument \"random_tree_stats\" should contain the columns \"regulator\", \"total_tree_length\" and \"within_species_diversity\".")

  actual <- tree_stats %>%
    dplyr::select(.data[["total_tree_length"]], .data[["within_species_diversity"]]) %>%
    as.matrix()

  mu_actual <- colMeans(actual)

  sigma_actual <- stats::cov(actual)

  random <- random_tree_stats %>%
    dplyr::select(.data[["total_tree_length"]], .data[["within_species_diversity"]]) %>%
    as.matrix()

  mu_random <- colMeans(random)

  sigma_random <- stats::cov(random)

  regulators_to_keep <- tree_stats %>%
    dplyr::select(.data[["regulator"]],
                  .data[["total_tree_length"]],
                  .data[["within_species_diversity"]]) %>%
    tidyr::pivot_longer(cols = c("total_tree_length", "within_species_diversity"), names_to = "statistic") %>%
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::summarize(d_actual = mvtnorm::dmvnorm(.data[["value"]], mean = mu_actual, sigma = sigma_actual),
                     d_random = mvtnorm::dmvnorm(.data[["value"]], mean = mu_random, sigma = sigma_random),
                     p_actual = .data[["d_actual"]] / (.data[["d_actual"]] + .data[["d_random"]])) %>%
    dplyr::filter(.data[["p_actual"]] > p_cutoff) %>%
    dplyr::pull(.data[["regulator"]])

  tree_stats %>%
    dplyr::filter(.data[["regulator"]] %in% regulators_to_keep)

}
