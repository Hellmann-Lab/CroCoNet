#' Find conserved and diverged target genes within a module
#'
#' Finds target genes within a conserved/diverged module that are particularly responsible for the conservation/divergence.
#'
#' To determine whether a module as whole is conserved, diverged overall or diverged on a specific lineage, the CroConet approach relies on module trees reconstrcuted from pairwise preservation scores between clones and statistics calculated based on these trees (total tree length, diversity, species-to-other branch length). To identify individual genes that contribute the most to conservation/divergence, the same statistics can be used in combination with jackknifing.
#'
#' During jackknifing, each member gene of a module is removed and all statistics are recalculated, including the tree-based statistics that inform us about cross-species conservation (see the parameter \code{jackknife} in the function \code{\link{calculatePresStats}}). Our working hypothesis is that if removing a target from a diverged module makes that module more conserved, then that target was responsible for divergence in the original module, and vice versa, if removing a target from a conserved module makes that module more diverged, then that target was responsible for conservation in the original module.
#'
#' To quantify these effects, the function takes use of the linear models that were fitted across all modules between the total tree length and within-species diversity (in case the focus of interest is conservation and overall divergence) or between the species-to-other branch length and the diversity within that species (in case the focus of interest is the species-specific divergence), similarly to how the conserved and diverged modules were identified in the first place. However, in this case it is not the aggregate statistic for a module as a whole that is compared to the regression line, but the statistic of each jackknife module version separately. The jackknife version that has the highest residual is the most diverged, therefore the corresponding target gene is the most conserved, while the jackknife version that has the lowest residual is the most conserved, therefore the corresponding target gene is the most diverged.
#'
#' It is important to note, that in most cases the divergence/conservation of a module cannot be contributed to a single target gene, but it is rather the combined signal of all targets together. For example, it is expected that most jackknife versions of a diverged module will still fall into the diverged region above the prediction interval of the regression line, even when the most diverged target genes are removed. This also means that the "most conserved" target gene of diverged module is still most likely diverged, just to a lesser extent than the others. Therefore, it mainly makes sense to investigate the most diverged targets of the diverged modules and the most conserved targets of the conserved modules, not the other way around.
#'
#' Based on the linear model provided, the function calculates the predicted value and its prediction interval at the diversity value of each jackknife module version as well as the residual of each jackknife module version. The width of the prediction interval is assumed to be constant for the entire range of diversity values across all jackknife versions of the same module and it is calculated based on the aggregate diversity value and (in case of a weighted linear model) the aggregate weight of the module as a whole (i.e. the same data that was used for fitting \code{lm_tree_stats}).
#' @param module_name Character, the name of the module of interest.
#' @param tree_stats_jk Data frame of tree statistics per jackknife module version across all modules.
#'
#' Columns required in case the focus of interest is conservation and overall divergence:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type (orig = original or jk = jackknifed).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig').}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type 'orig').}
#' \item{total_tree_length}{Numeric, total tree length per jackknife module version.}
#' \item{within_species_diversity}{Numeric, within-species diversity per jackknife module version.}
#' }
#' Columns required in case the focus of interest is species-specific divergence:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type (orig = original or jk = jackknifed).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig').}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type 'orig').}
#' \item{\{\{species\}\}_to_other_branch_length}{Numeric, the branch length between the species of interest and all others per jackknife module version}
#' \item{\{\{species\}\}_diversity}{Numeric, diversity within the species of interest per jackknife module version.}
#' }
#' @param lm_tree_stats Object of class \code{\link{lm}}, the (weighted) linear model fit between the total tree length and within-species diversity (in case the focus of interest is conservation and overall divergence) or between the species-to-other branch length and diversity of a species (in case the focus of interest is the divergence between that particular species and all others).
#' @param conf_level Confidence level of the prediction interval (default: 0.95).
#'
#' @return Data frame of cross-species conservation measures per target gene for the module of interest. In addition to the relevant columns of \code{tree_stats_jk}, it contains 4 new columns:
#' \describe{
#' \item{focus}{Character, the focus of interest in terms of cross-species conservation, either "overall" if the focus of interest is conservation and overall divergence, or the name of a species if the focus of interest is the divergence between that particular species and all others.}
#' \item{fit}{Numeric, the fitted total tree length/species-to-other branch length at the diversity value of a jackknife module version.}
#' \item{lwr_fit}{Numeric, the lower bound of the prediction interval of the fit.}
#' \item{upr_fit}{Numeric, the upper bound of the prediction interval of the fit.}
#' \item{residual}{Numeric, the residual of the jackknife module version in the linear model. It is calculated as the difference between the observed and expected (fitted) total tree length/species-to-other branch length.}
#' }
#' It also contains a new summary row with the aggregate values across all jackknife values of the module.
#' @export
#'
#' @examples POU5F1_target_conservation <- findConservedDivergedTargets("POU5F1", tree_stats_jk, lm_overall)
findConservedDivergedTargets <- function(module_name, tree_stats_jk, lm_tree_stats, conf_level = 0.95) {

  if (!is.data.frame(tree_stats_jk))
    stop("The argument \"tree_stats_jk\" should be a data frame.")

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

  if (any(!c("regulator", "module_size", "type", "id", "gene_removed") %in% colnames(tree_stats_jk)))
    stop("The argument \"tree_stats_jk\" should contain the columns \"regulator\", \"module_size\", \"type\", \"id\" and \"gene_removed\".")

  if (any(!(c(x_var, y_var) %in% colnames(tree_stats_jk))))
    stop(paste0("The argument \"tree_stats_jk\" should contain the columns corresponding to the dependent and independent variables in \"lm_tree_stats\" (\"", x_var, "\" and \"", y_var, "\"). Please include these in the data frame or rerun the function \"fitTreeStatsLm()\" with the variables of interest."))

  if (!inherits(module_name, "character") || length(module_name) != 1)
    stop("The argument \"module_name\" should be a string specifying a single module in which the conserved and diverged target genes should be identified.")

  if (!module_name %in% unique(tree_stats_jk$regulator))
    stop("Tree statistics for the specified module cannot be found in \"tree_stats_jk\". Please make sure that the argument \"module_name\" is present in the column \"regulator\" of \"tree_stats_jk\".")

  if (!inherits(conf_level, "numeric") || length(conf_level) != 1 || conf_level < 0 || conf_level >= 1)
    stop("The argument \"conf_level\" should be a numeric value between 0 and 1.")

  # focus of interest
  focus <- lm_tree_stats$focus

  # variable for which the lm was fit
  y_var <- lm_tree_stats$terms[[2]]
  x_var <- lm_tree_stats$terms[[3]]

  # subset module
  tree_stats_jk_mod <- tree_stats_jk %>%
    dplyr::filter(.data[["regulator"]] == module_name) %>%
    dplyr::select(.data[["regulator"]], .data[["module_size"]], .data[["type"]], .data[["id"]], .data[["gene_removed"]], .data[[x_var]], .data[[y_var]]) %>%
    dplyr::bind_rows(lm_tree_stats$model %>%
                       dplyr::filter(.data[["regulator"]] == module_name) %>%
                       dplyr::transmute(.data[["regulator"]],
                                        type = "summary",
                                        id = paste0(module_name, "_summary"),
                                        .data[["module_size"]],
                                        gene_removed = NA,
                                        .data[[x_var]],
                                        .data[[y_var]])) %>%
    dplyr::mutate(focus = focus, .before = 1)

  # slope and intercept
  intercept <- lm_tree_stats$coefficients[1]
  slope <- lm_tree_stats$coefficients[2]

  # weighted fit
  if (!is.null(lm_tree_stats$weight)) {

    # get weight for original module
    weight <- lm_tree_stats$weights[lm_tree_stats$model$regulator == module_name]

    # 95%  prediction interval of the fit
    pred_interval <- stats::predict(lm_tree_stats,
                                    interval = "prediction",
                                    level = conf_level,
                                    se.fit = TRUE,
                                    newdata = tree_stats_jk_mod,
                                    weights = weight)

    # rank targets based on the effect of removing them (inverse logic: if removing a gene shifts the module towards the conserved section of the plot, that target is likely to be diverged; analogously for the rest of the categories)
    target_conservation <- cbind(tree_stats_jk_mod,
                                 pred_interval %>%
                                  as.data.frame() %>%
                                   dplyr::select(fit = .data[["fit.fit"]],
                                                 lwr_fit = .data[["fit.lwr"]],
                                                 upr_fit = .data[["fit.upr"]])) %>%
      dplyr::mutate(residual = .data[[y_var]] - .data[["fit"]])

  } else {

    # 95%  prediction interval of the fit
    pred_interval <- stats::predict(lm_tree_stats,
                                    interval = "prediction",
                                    level = conf_level,
                                    se.fit = TRUE,
                                    newdata = tree_stats_jk_mod)

    # rank targets based on the effect of removing them (inverse logic: if removing a gene shifts the module towards the conserved section of the plot, that target is likely to be diverged; analogously for the rest of the categories)
    target_conservation <- cbind(tree_stats_jk_mod,
                                 pred_interval %>%
                                   as.data.frame() %>%
                                   dplyr::select(fit = .data[["fit.fit"]],
                                                 lwr_fit = .data[["fit.lwr"]],
                                                 upr_fit = .data[["fit.upr"]])) %>%
      dplyr::mutate(residual = .data[[y_var]] - .data[["fit"]])

  }

  target_conservation

}
