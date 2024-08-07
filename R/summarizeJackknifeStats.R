#' Summarize jackknife statistics
#'
#' Calculates the estimate and confidence interval of a module-level statistic based on all values obtained by jackknifing the module.
#'
#' To gain information about the confidence of various statistics, jackknifing can be used: each member gene of a module is removed and the statistics are re-calculated. This way, the median or mean and its confidence interval across the jackknifed versions can be used to estimate a statistic of interest for the module as a whole. For continuous statistics the median, for boolean statistics the mean is recommended as the summary method.
#' @param stats_df Data frame of the statistics per jackknife module version:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulato.r}
#' \item{type}{Character, module type (orig = original or jk = jackknifed).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig').}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type 'orig').}
#' \item{clone}{Character, the name of the clone (optional, only needed if the statistics are calculated per clone).}
#' \item{species}{Character, the name of the species (optional, only needed if the statistics are calculated per clone or per species).}
#' \item{clone1, clone2}{Character, the names of the clones compared (optional, only needed if the statistics are calculated per clone pair).}
#' \item{species1, species2}{Character, the names of the species compared (optional, only needed if the statistics are calculated per clone pair or species pair).}
#' \item{\{\{nameOfStat\}\}}{Numeric, integer or logical, one or more columns containing the values of the statistic(s) specified in 'stats'.}
#' }
#' @param stats Character, the name of the statistic(s) that need to be summarized.
#' @param summary_method Character or character vector, the summary method ("mean" or "median") to be used for the statistics in 'stats'. If only one value is provided, the same summary method will be used for all statistics, if a vector is provided, each element of the vector will be used for the corresponding element in 'stats'. By default, "mean" for statistics measuring monophyleticity and "median" for all other statistics.
#' @param conf_level Numeric, confidence level of the interval (default: 0.95).
#'
#' @return A data frame of the statistics per module:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Module size, the numer of target genes assigned to a regulator.}
#' \item{clone}{The name of the clone (only if the column is present in the input 'stats_df').}
#' \item{species}{The name of the species (only if the column is present in the input 'stats_df').}
#' \item{clone1, clone2}{The names of the clones compared (only if the column is present in the input 'stats_df').}
#' \item{species1, species2}{The names of the species compared (only if the column is present in the input 'stats_df').}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the estimates (mean or median) of the statistic(s) specified in 'stats'.}
#' \item{var_\{\{nameOfStat\}\}}{Numeric, one or more columns containing the variances of the statistic(s) specified in 'stats'.}
#' \item{lwr_\{\{nameOfStat\}\}}{Numeric, one or more columns containing the lower bounds of the confidence interval for the statistic(s) specified in 'stats'.}
#' \item{upr_\{\{nameOfStat\}\}}{Numeric, one or more columns containing the upper bounds of the confidence interval for the statistic(s) specified in 'stats'.}
#' }
#' @export
#'
#' @examples
#' pres_stats <- summarizeJackknifeStats(pres_stats_jk)
#' tree_stats <- summarizeJackknifeStats(tree_stats_jk,
#'                                       c("total_tree_length", "within_species_diversity"))
summarizeJackknifeStats <- function(stats_df, stats = setdiff(colnames(stats_df), c("regulator", "type", "id", "gene_removed", "module_size", "clone", "species", "clone1", "clone2", "species1", "species2")), summary_method = ifelse(endsWith(stats, "monophyl"), "mean", "median"), conf_level = 0.95) {

  if (!is.data.frame(stats_df))
    stop("The argument \"stats_df\" should be a data frame.")

  if (any(!c("regulator", "module_size", "type", "id", "gene_removed") %in% colnames(stats_df)))
    stop("The argument \"stats_df\" should contain the columns \"regulator\", \"module_size\", \"type\", \"id\", and \"gene_removed\".")

  if (!inherits(stats, "character"))
    stop("The argument \"stats\" should be a character or character vector.")

  missing_in_stats_df <- stats[!stats %in% colnames(stats_df)]

  if (length(missing_in_stats_df) > 0)
    stop(paste0("The following statistic(s) in \"stats\" were not found in \"stats_df\": ", paste(missing_in_stats_df, collapse = ", "), "."))

  if (is.null(summary_method) || !all(summary_method %in% c("mean", "median")))
    stop("The argument \"summary_method\" should be a character vector where each element is either \"mean\" or \"median\".")

  if (length(summary_method) != length(stats))
    stop("The arguments \"stats\" and \"summary_method\" should have the same length.")

  if (!inherits(conf_level, "numeric") || length(conf_level) != 1 || conf_level < 0 || conf_level >= 1)
    stop("The argument \"conf_level\" should be a numeric value between 0 and 1.")

  names(summary_method) <- stats

  stats_df[colnames(stats_df) %in% c(stats, c("regulator", "type", "id", "gene_removed", "module_size", "clone", "species", "clone1", "clone2", "species1", "species2"))] %>%
    dplyr::filter(.data[["type"]] != "orig") %>%
    dplyr::select(-.data[["type"]]) %>%
    tidyr::pivot_longer(cols = stats, names_to = "statistic", values_to = "jk") %>%
    dplyr::group_by(dplyr::across(-c("id", "gene_removed", "jk"))) %>%
    dplyr::summarise(summarizeStat(.data[["jk"]], summary_method[unique(.data[["statistic"]])], conf_level)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider_spec(data.frame(statistic = stats) %>%
                              tidyr::expand_grid(.value = c("estimate", "var", "lwr", "upr")) %>%
                              dplyr::mutate(.name = gsub("estimate_", "", paste0(.data[[".value"]], "_", .data[["statistic"]]))))

}



#' Summarize values of a statistic
#'
#' @description Calculates the estimate (mean or median), the confidence interval of the estimate with the specified confidence level and the variance of the data provided.
#' @param values Numeric, integer or logical vector.
#' @param summary_method Character, the measure of central tendency ("mean" or "median") to be used.
#' @param conf_level Numeric, confidence level of the interval (default: 0.95).
#'
#' @return A data frame with 1 row and 4 columns:
#' \describe{
#' \item{estimate}{Numeric, the central tendency (mean or median) of 'values'.}
#' \item{var}{Numeric, the variance of 'values'.}
#' \item{lwr}{Numeric, the lower bound of the confidence interval.}
#' \item{upr}{Numeric, the upper bound of the confidence interval.}
#' }
#' @export
summarizeStat <- function(values, summary_method, conf_level = 0.95) {

  if (summary_method == "mean") {

    n <- length(values)
    moe_values <- -stats::qt(p = (1 - conf_level)/2, df = n - 1) * stats::sd(values) / sqrt(n)
    mean_values <- mean(values)

    data.frame(estimate = mean_values,
               var = stats::var(values),
               lwr = mean_values - moe_values,
               upr = mean_values + moe_values)

  } else {

    median_values <- stats::median(values)

    if (is.na(median_values)) {

      ci <- c(NA, NA)

    } else {

      n <- length(values)
      k <- stats::qbinom(p = (1 - conf_level) / 2, size = n, prob = 0.5, lower.tail = TRUE)
      ci <- sort(values)[c(k, n - k + 1)]
      if (identical(ci, NA_real_)) ci <- c(-Inf, Inf)

    }

    data.frame(estimate = median_values,
               var = stats::var(values),
               lwr = ci[1],
               upr = ci[2])

  }

}

