#' Find conserved and diverged modules
#'
#' Finds conserved and diverged network modules based on a (weighted) linear model between the total tree length and within-species diversity (in case the focus of interest is conservation and overall divergence) or between the species-to-other branch length and diversity of a species (in case the focus of interest is species-specific divergence) of the module trees. Modules are considered conserved or diverged if they fall outside the 95\% prediction interval of the regression line.
#'
#'  As part of the CroCoNet approach, pairwise module preservation scores are calculated between clones, both within and across species (see \code{\link{calculatePresStats}}) and neighbor-joining trees are reconstructed based on these preservation scores per module (see \code{\link{reconstructTrees}}). The tips of the resulting tree represent the clones and the branch lengths represent the dissimilarity of module connectivity patterns between the networks of 2 clones.
#'
#' Various useful statistics can be defined based on these trees (see also \code{\link{calculateTreeStats}}). The total tree length is the sum of all branch lengths in the tree and it measures module variability both within and across species. The diversity of a species is the sum of the branches connecting the clones of this species and it measures module variability within this particular species. The within-species diversity is the sum of the diversity values across all species and it measures module variability within species in general. The species-to-other branch length is defined as the length of the internal branch that connects the subtree of the clones from the species of interest to the subtree of all other clones and it measures module variability between this species and all others.
#'
#' Overall divergence is best represented by the total tree length, while species-specific divergence is best represented by the species-to-other branch length, however the absolute values of these statistics in themselves are difficult to interpret. The diversity (either overall or for a particular species) can be used as an internal reference to determine whether the total tree length or the species-to-other branch length is larger/smaller than expected. The general relationship between total tree length and overall diversity or between the species-to-other branch length and the diversity of the species is captured via linear regression (see also \code{\link{fitTreeStatsLm}}). In term of cross-species conservation, the most interesting modules are then the outliers that do not follow the linear trend.
#'
#' By analyzing the linear regression between total tree length and within-species diversity, we can single out modules that are conserved or diverged overall. The modules that fall far below the trend line are the most conserved ones: they have a lower total tree length than expected based on their within-species diversity, meaning that the species are not well-separated within the tree but rather mixed among each other. In contrast, the modules that are located far above the trend line are the most diverged ones: they have a higher total tree length than expected based on their within-species diversity, meaning that they thave long cross-species branches between one or more pairs of species.
#'
#' By analyzing the linear regression between the species-to-other branch length and the diversity of a species, we can single out modules that are diverged between this particular species and all others: these are the modules that fall far above the trend line, i.e. that have a longer species-to-other branch than expected based on the diversity of the species. The modules that fall below the trend line are in this case not meaningful, since we can only call a module conserved, if it is conserved across all species.
#'
#' For both types of analysis, we quantitatively define which modules are outliers, and thus conserved or diverged, by calculating the prediction interval of the linear fit. A module is considered diverged overall/in a species-specific manner if it has a higher total tree length/species-to-other branch length than the upper boundary of the prediction interval, while a module is considered conserved if it has a lower total tree length than the lower boundary of the prediction interval.
#'
#' The degree of conservation/divergence can be further compared between the modules categorized as conserved/diverged. We recommend 2 measures: 1) the residual, which is the absolute difference between the observed and expected total tree lengths/species-to-other branch lengths at a given diversity value, and 2) the t-score, which is the residual normalized by the standard error of the total tree length/species-to-other brnach length prediction at a given diversity value.
#'
#' The linear regression can be weighted by the error of the data points. To gain this information, we can use jackknifing: each member gene of a module is removed and all statistics are recalculated (see the argument "jackknife" in the functions \code{\link{calculatePresStats}}, \code{\link{reconstructTrees}} and \code{\link{calculateTreeStats}}. The weight of a module in the regression is then defined to be inversely proportional to the variance of the dependent variable (total tree length/species-to-other branch length) across all jackknife versions.
#' @param tree_stats Data frame of the tree-based statistics for the pruned modules.
#'
#' Columns required in case the focus of interest is conservation and overall divergence:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{total_tree_length}{Numeric, total tree length per module (typically the median across all jackknife versions of the module).}
#' \item{var_total_tree_length}{Numeric, the variance of total tree lengths across all jackknifed versions of the module (optional, only needed for weighted regression).}
#' \item{lwr_total_tree_length}{Numeric, the lower bound of the confidence interval of the total tree length calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{upr_total_tree_length}{Numeric, the upper bound of the confidence interval of the total tree length calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{within_species_diversity}{Numeric, within-species diveristy per module (typically the median across all jackknife versions of the module).}
#' \item{lwr_within_species_diversity}{Numeric, the lower bound of the confidence interval of the within-species diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{upr_within_species_diversity}{Numeric, the upper bound of the confidence interval of the within-species diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' }
#' Columns required in case the focus of interest is species-specific divergence:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{\{\{species\}\}_to_other_branch_length}{Numeric, the species-to-other branch length per module (typically the median across all jackknife versions of the module).}
#' \item{var_\{\{species\}\}_to_other_branch_length}{Numeric, the variance of the species-to-other branch lengths across all jackknifed versions of the module (optional, only needed for weighted regression).}
#' \item{lwr_\{\{species\}\}_to_other_branch_length}{Numeric, the lower bound of the confidence interval of the species-to-other branch length calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{upr_\{\{species\}\}_to_other_branch_length}{Numeric, the upper bound of the confidence interval of the species-to-other branch length calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{\{\{species\}\}_diversity}{Numeric, the diveristy of the species of interest per module (typically the median across all jackknife versions of the module)}
#' \item{lwr_\{\{species\}\}_diversity}{Numeric, the lower bound of the confidence interval of the species-specific diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{upr_\{\{species\}\}_diversity}{Numeric, the upper bound of the confidence interval of the species-specific diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' }
#' @param lm_tree_stats An object of class \code{\link{lm}}, the (weighted) linear model between the total tree length and within-species diversity (in case the focus of interest is conservation and overall divergence) or between the species-to-other branch length and diversity of a species (in case the focus of interest is divergence between this particular species and all others).
#' @param conf_level Confidence level of the prediction interval (default: 0.95).
#'
#' @return Data frame of cross-species conservation measures per module. In addition to the relevant columns of the input "tree_stats", it contains 6 or 7 new columns:
#' \describe{
#' \item{focus}{Character, the focus of interest in terms of cross-species conservation, either "overall" if the focus of interest is conservation and overall divergence, or the name of a species if the focus of interest is the divergence between this particular species and all others. Derived from the input linear model (\code{lm_tree_stats}).}
#' \item{fit}{Numeric, the fitted total tree length/species-to_other branch length at the diversity value of the module.}
#' \item{lwr_fit}{Numeric, the lower bound of the prediction interval of the fit.}
#' \item{upr_fit}{Numeric, the upper bound of the prediction interval of the fit.}
#' \item{residual}{Numeric, the residual of the module in the linear model. It is calculated as the difference between the observed and expected (fitted) total tree lengths/species-to-other branch lengths.}
#' \item{weight}{Numeric, the weight of the module in the linear regression, inversely proportional to the variance of total tree lengths/species-to-other branch lengths (only if \code{weighted_lm} is set to TRUE).}
#' \item{t_score}{Numeric, the t-score of the module. It is calculated as the residual normalized by the standard error of the total tree length/species-to-other branch length prediction at the given diversity value.}
#' \item{conservation}{Character, "not_significant" if the module falls inside the prediction interval of the fit, "diverged" if a module has a higher total tree length/species-to-other branch length than the upper boundary of the prediction interval, and "conserved" if a module has a lower total tree length/species-to-other branch length than the lower boundary of the prediction interval.}
#' }
#' @export
#'
#' @examples module_conservation_overall <- findConservedDivergedModules(tree_stats, lm_overall)
findConservedDivergedModules <- function(tree_stats, lm_tree_stats, conf_level = 0.95) {

  if (!is.data.frame(tree_stats))
    stop("The argument \"tree_stats\" should be a data frame.")

  if (!inherits(lm_tree_stats, "lm"))
    stop("The argument \"lm_tree_stats\" should be an object of class \"lm\".")

  # focus of interest
  focus <- lm_tree_stats$focus

  if (is.null(focus))
    stop("Please specify the focus of analysis by adding a component \"focus\" to \"lm_tree_stats\". This should be \"overall\" if the focus of interest is conservation and overall divergence and the name of a species if the focus of interest is species-specific divergence.")

  if (!inherits(focus, "character") || length(focus) != 1)
    stop("The component \"focus\" of \"lm_tree_stats\" should be \"overall\" if the focus of interest is conservation and overall divergence and the name of a species if the focus of interest is species-specific divergence.")

  # variable for which the lm was fit
  y_var <- lm_tree_stats$terms[[2]]
  x_var <- lm_tree_stats$terms[[3]]

  if (focus == "overall") {

    if (x_var != "within_species_diversity" || y_var != "total_tree_length")
      stop("The component \"focus\" of \'lm_tree_stats\" should match the independent and dependent variables in the linear model. If \"focus\" is set to \"overall\", the independent variable should be \"within_species_diversity\" and the dependent variable should be \"total_tree_length\". If \"focus\" is the name of a species, the independent variable should be \"{{species}}_diversity\" and the dependent variable should be \"{{species}}_to_other_branch_length\".")

  } else {

    if (x_var != paste0(focus, "_diversity") || y_var != paste0(focus, "_to_other_branch_length"))
      stop("The component \"focus\" of \'lm_tree_stats\" should match the independent and dependent variables in the linear model. If \"focus\" is set to \"overall\", the independent variable should be \"within_species_diversity\" and the dependent variable should be \"total_tree_length\". If \"focus\" is the name of a species, the independent variable should be \"{{species}}_diversity\" and the dependent variable should be \"{{species}}_to_other_branch_length\".")

  }

  if (!"regulator" %in% colnames(tree_stats))
   stop(paste0("The argument \"tree_stats\" should contain the column \"regulator\"."))

  if (any(!(c(x_var, y_var) %in% colnames(tree_stats))))
    stop(paste0("The argument \"tree_stats\" should contain the columns corresponding to the dependent and independent variables in \"lm_tree_stats\" (\"", x_var, "\" and \"", y_var, "\"). Please include these in the data frame or rerun the function \"fitTreeStatsLm()\" with the variables of interest."))

  if (!inherits(conf_level, "numeric") || length(conf_level) != 1 || conf_level < 0 || conf_level >= 1)
    stop("The argument \"conf_level\" should be a numeric value between 0 and 1.")

  # 95%  prediction interval of the fit
  pred_interval <- suppressWarnings(
    stats::predict(lm_tree_stats,
                   interval = "prediction",
                   level = conf_level,
                   se.fit = TRUE)
  )

  # subset data
  tree_stats <- tree_stats[, colnames(tree_stats) %in% c("regulator", "module_size", x_var, paste0("lwr_", x_var), paste0("upr_", x_var), y_var, paste0("lwr_", y_var), paste0("upr_", y_var))] %>%
    dplyr::mutate(focus = focus, .before = 1)

  # get residuals and t-scores for weighted fit
  if (!is.null(lm_tree_stats$weight)) {

    module_conservation <- cbind(tree_stats,
                                 pred_interval %>%
                                   as.data.frame() %>%
                                   dplyr::select(fit = .data[["fit.fit"]],
                                                 lwr_fit = .data[["fit.lwr"]],
                                                 upr_fit = .data[["fit.upr"]],
                                                 .data[["residual.scale"]],
                                                 .data[["se.fit"]])) %>%
      dplyr::mutate(residual = lm_tree_stats$residuals,
                    weight = lm_tree_stats$weight,
                    se_tot = (.data[["residual.scale"]]^2/.data[["weight"]] + .data[["se.fit"]]^2)^0.5,
                    t_score = .data[["residual"]] / .data[["se_tot"]]) %>%
      dplyr::select(-.data[["residual.scale"]], -.data[["se.fit"]], -.data[["se_tot"]])

  # get residuals and t-scores for unweighted fit
  } else {

    # get conserved and diverged modules
    module_conservation <- cbind(tree_stats,
                                 pred_interval %>%
                                   as.data.frame() %>%
                                   dplyr::select(fit = .data[["fit.fit"]],
                                                 lwr_fit = .data[["fit.lwr"]],
                                                 upr_fit = .data[["fit.upr"]],
                                                 .data[["residual.scale"]],
                                                 .data[["se.fit"]])) %>%
      dplyr::mutate(residual = lm_tree_stats$residuals,
                    se_tot = (.data[["residual.scale"]]^2 + .data[["se.fit"]]^2)^0.5,
                    t_score = .data[["residual"]] / .data[["se_tot"]]) %>%
      dplyr::select(-.data[["residual.scale"]], -.data[["se.fit"]], -.data[["se_tot"]])

  }

  # label conserved and diverge modules in case the focus of interest is conservation and overall divergence
  if (focus == "overall") {

    module_conservation <- module_conservation %>%
      dplyr::mutate(conservation = dplyr::case_when(.data[[y_var]] > .data[["upr_fit"]] ~ "diverged",
                                                    .data[[y_var]] < .data[["lwr_fit"]] ~ "conserved",
                                                    TRUE ~ "not_significant"))

  # label diverged modules in case the focus of interest is species-specific divergence
  } else {

    module_conservation <- module_conservation %>%
      dplyr::mutate(conservation = dplyr::case_when(.data[[y_var]] > .data[["upr_fit"]] ~ "diverged",
                                                    TRUE ~ "not_significant"))

  }

  module_conservation

}
