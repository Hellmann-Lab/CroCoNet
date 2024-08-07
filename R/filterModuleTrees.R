#' Filter modules based on their tree representations
#'
#' Filters modules using two tree-based statistics (the total tree length and within-species diversity) by comparing the values of the actual modules to the values of the random modules.
#'
#' The within-species diversity measures the variability of module connectivity patterns within species and the total tree length measures the variability of module connectivity patterns both within and across species. The within-species diversity is a subset of the total tree, hence we expect a linear relationship between the two statistics. This line explains the detection robustness: modules that have both low within-species diversity and low total tree length are well-preserved both within and across species, meaning that these modules could be robustly detected in all clones, whereas modules, for which both metrics are high, are poorly preserved not just across but also within species, indicating a high detection uncertainty ("wobbliness"). Random modules are expected to fall at the wobbly end of the spectrum, while actual modules are expected to fall at the robust end of the spectrum.
#'
#' The function removes the modules that are uncertain to detect, i.e. modules the are too similar to the random modules in terms of their total tree length and within-species diversity. This similarity is measured by the probability that the statistic of an actual modules comes from the distribution of the statistics of all random modules (one-sample one-sided z-test with the mean and standard deviation of the random module statistics as the expected mean and standard deviation). If the adjusted p-value of a module is below the specified \code{p_adj_cutoff} for either of the two statistics, i.e. the module is significantly different from the random modules, it is kept.
#' @param tree_stats_actual Data frame of the tree-based statistics for the actual (pruned) modules containing (but not restricted to) the following columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{total_tree_length}{Numeric, total tree length per module (typically the median across all jackknife versions of a module.}
#' \item{within_species_diversity}{Numeric, within-species diveristy per module (typically the median across all jackknife versions of a module.}
#' }
#' @param tree_stats_random Data frame of the tree-based statistics for the random modules containing (but not restricted to) the following columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{total_tree_length}{Numeric, total tree length per module (typically the median across all jackknife versions of a module).}
#' \item{within_species_diversity}{Numeric, within-species diveristy per module (typically the median across all jackknife versions of a module).}
#' }
#' @param p_adj_cutoff Numeric, adjusted p-value cutoff (default: 0.1).
#'
#' @return Data frame of the tree-based statistics for the actual (pruned) modules after filtering.
#' @export
filterModuleTrees <- function(tree_stats_actual, tree_stats_random, p_adj_cutoff = 0.1) {

  if (!is.data.frame(tree_stats_actual))
    stop("The argument \"tree_stats_actual\" should be a data frame.")

  if (!is.data.frame(tree_stats_random))
    stop("The argument \"tree_stats_random\" should be a data frame.")

  if (any(!(c("regulator", "total_tree_length", "within_species_diversity") %in% colnames(tree_stats_actual))))
    stop("The argument \"tree_stats_actual\" should contain the columns \"regulator\", \"total_tree_length\" and \"within_species_diversity\".")

  if (any(!(c("regulator", "total_tree_length", "within_species_diversity") %in% colnames(tree_stats_random))))
    stop("The argument \"tree_stats_random\" should contain the columns \"regulator\", \"total_tree_length\" and \"within_species_diversity\".")

  if (!inherits(p_adj_cutoff, "numeric") || length(p_adj_cutoff) != 1 || p_adj_cutoff < 0 || p_adj_cutoff >= 1)
    stop("The argument \"p_adj_cutoff\" should be a numeric value between 0 and 1.")

  regulators_to_keep <- dplyr::bind_rows(actual = tree_stats_actual,
                                         random = tree_stats_random,
                                         .id = "module_set") %>%
    dplyr::select(.data[["module_set"]], .data[["regulator"]],
                  .data[["total_tree_length"]],
                  .data[["within_species_diversity"]]) %>%
    tidyr::pivot_longer(cols = c("total_tree_length", "within_species_diversity"), names_to = "statistic") %>%
    tidyr::pivot_wider(names_from = "module_set", values_from = "value") %>%
    dplyr::group_by(.data[["statistic"]]) %>%
    dplyr::mutate(mean_random = mean(.data[["random"]]),
                  sd_random = stats::sd(.data[["random"]])) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p_random = stats::pnorm(.data[["actual"]], .data[["mean_random"]], .data[["sd_random"]])) %>%
    dplyr::mutate(p_adj_random = stats::p.adjust(.data[["p_random"]], method = "bonferroni")) %>%
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::filter(sum(.data[["p_adj_random"]] < p_adj_cutoff) > 0) %>%
    dplyr::pull(.data[["regulator"]]) %>%
    unique()

  tree_stats_actual %>%
    dplyr::filter(.data[["regulator"]] %in% regulators_to_keep)

}
