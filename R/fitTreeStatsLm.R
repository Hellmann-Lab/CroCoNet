#' Fit (weighted) linear model between tree-based statistics
#'
#' Fits a (weighted) linear model between the total tree length and within-species diversity (in case the focus of interest is conservation and overall divergence) or between the subtree length and the diversity of a species (in case the focus of interest is the species-specific divergence).
#'
#' The linear models output by this function can be used to identify conserved and diverged modules, and to identify target genes within these modules that contribute the most to the conservation/divergence. For details, please see \code{\link{findConservedDivergedModules}} and \code{\link{findConservedDivergedTargets}}.
#'
#' The focus of interest can be specified using the parameter \code{focus}. If \code{focus} is set to "overall" (default), the linear model will be fit between the total tree length and within-species diversity, and subsequent analysis using \code{\link{findConservedDivergedModules}} and \code{\link{findConservedDivergedTargets}} can identify modules and target genes that are conserved or diverged across all species. If \code{focus} is set to the name of a species in the dataset, the linear model will be fit between the subtree length and the diversity of that species, and subsequent analysis using \code{\link{findConservedDivergedModules}} and \code{\link{findConservedDivergedTargets}} can identify modules and target genes that are diverged between the species and all others. Please note that if the aim is to find conserved modules, \code{focus} should always be set to "overall".
#'
#' The function fits the linear model corresponding to the focus of interest by calling \code{\link{lm}}. If a weighted model is desired (default), the weights are defined to be inversely proportional to the variance of the dependent variable (total tree length or the subtree length of a species). If no jackknifing was performed and thus the variance is unknown, please set \code{weighted_lm} to FALSE.
#' @param tree_stats Data frame of the tree-based statistics for the pruned modules.
#'
#' Columns required in case the focus of interest is conservation and overall divergence:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{total_tree_length}{Numeric, total tree length per module (typically the median across all jackknife versions of the module).}
#' \item{var_total_tree_length}{Numeric, the variance of total tree lengths across all jackknifed versions of the module (optional, only needed for weighted regression).}
#' \item{within_species_diversity}{Numeric, within-species diveristy per module (typically the median across all jackknife versions of the module).}
#' }
#' Columns required in case the focus of interest is species-specific divergence:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{\{\{species\}\}_subtree_length}{Numeric, the sum of the branch lengths in the subtree that is defined by the replicates of the species and includes the internal branch connecting these replicates to the rest of the tree (typically the median across all jackknife versions of the module).}
#' \item{var_\{\{species\}\}_subtree_length}{Numeric, the variance of the subtree length of a species across all jackknifed versions of the module (optional, only needed for weighted regression).}
#' \item{\{\{species\}\}_diversity}{Numeric, within-species diveristy per module (typically the median across all jackknife versions of the module).}
#' }
#' @param focus Character, the focus of interest in terms of cross-species conservation, either "overall" if the focus of interest is conservation and overall divergence, or the name of a species if the focus of interest is the divergence between that particular species and all others.
#' @param weighted_lm Logical indicating whether the linear regression should be weighted or not (default: TRUE). If TRUE, \code{tree_stats} is expected to contain the column \code{var_total_tree_length} and the weights will be inversely proportional to these variances. If no jackknifing was performed and thus these variances were not calculated, please set this parameter to FALSE.
#'
#' @return An object of class \code{\link{lm}}.
#' @export
#'
#' @examples lm_overall <- fitTreeStatsLm(tree_stats, focus = "overall")
fitTreeStatsLm <- function(tree_stats, focus = "overall", weighted_lm = TRUE) {

  if (!is.data.frame(tree_stats))
    stop("The argument \"tree_stats\" should be a data frame.")

  if (any(!c("regulator", "module_size") %in% colnames(tree_stats)))
    stop("The argument \"tree_stats\" should contain the columns \"regulator\" and \"module_size\".")

  species <- lapply(strsplit(colnames(tree_stats)[grepl("_subtree_length", colnames(tree_stats))], "_|_subtree_length"),
                    function(x) {x[length(x)]}) %>%
    unlist() %>%
    unique()

  if (!inherits(focus, "character") || length(focus) != 1 || !focus %in% c("overall", species))
    stop("The argument \"focus\" should be \"overall\" if the focus of interest is conservation and overall divergence or the name of a species if the focus of interest is species-specific divergence. If you set it to a species, please make sure that the columns \"{{species}}_diversity\" and \"{{species}}_subtree_length\" are present in \"tree_stats\".")

  if (!inherits(weighted_lm, "logical") || length(weighted_lm) != 1)
    stop("The argument \"weighted_lm\" should be a logical value.")

  if (focus == "overall") {

    if (any(!(c("within_species_diversity", "total_tree_length") %in% colnames(tree_stats))))
      stop("The argument \"tree_stats\" should contain the columns corresponding to \"focus\". If \"focus\" is set to \"overall\", the required columns are \"within_species_diversity\" and \"total_tree_length\".")

    if (weighted_lm && !"var_total_tree_length" %in% colnames(tree_stats))
      stop("If \"focus\" is set to \"overall\" and \"weighted_lm\" is set to TRUE, \"tree_stats\" is expected to contain the column \"var_total_tree_length\". Please add this column to the data frame or change the value of \"weighted_lm\" or \"focus\".")

  } else {

    if (any(!(c(paste0(focus, "_diversity"), paste0(focus, "_subtree_length")) %in% colnames(tree_stats))))
      stop(paste0("The argument \"tree_stats\" should contain the columns corresponding to \"focus\". If \"focus\" is set to \"", focus, "\", the required columns are \"", focus, "_diversity\" and \"", focus, "_subtree_length\"."))

    if (weighted_lm && !paste0("var_", focus, "_subtree_length") %in% colnames(tree_stats))
      stop(paste0("If \"focus\" is set to \"", focus, "\" and \"weighted_lm\" is set to TRUE, \"tree_stats\" is expected to contain the column \"var_", focus, "_subtree_length\". Please add this column to the data frame or change the value of \"weighted_lm\" or \"focus\"."))

  }

  weight = regulator = NULL

  if (focus == "overall") {

    # weighted fit
    if (weighted_lm) {

      # get weights
      tree_stats <- tree_stats %>%
        tidyr::drop_na(.data[["var_total_tree_length"]]) %>%
        dplyr::mutate(weight = 1 / .data[["var_total_tree_length"]],
                      weight = .data[["weight"]] / sum(.data[["weight"]]))

      # fit weighted lm
      lm_tree_stats <- stats::lm(total_tree_length ~ within_species_diversity, tree_stats, weights = weight)

    } else {

      # fit weighted lm
      lm_tree_stats <- stats::lm(total_tree_length ~ within_species_diversity, tree_stats)

    }

  } else {

    # weighted fit
    if (weighted_lm) {

      # get weights
      tree_stats <- tree_stats %>%
        tidyr::drop_na(.data[[paste0("var_", focus, "_subtree_length")]]) %>%
        dplyr::mutate(weight = 1 / .data[[paste0("var_", focus, "_subtree_length")]],
                      weight = .data[["weight"]] / sum(.data[["weight"]]))

      # fit weighted lm
      lm_tree_stats <- stats::lm(stats::as.formula(paste0(focus, "_subtree_length ~ ", focus, "_diversity")), tree_stats, weights = weight)

    } else {

      # fit lm
      lm_tree_stats <- stats::lm(stats::as.formula(paste0(focus, "_subtree_length ~ ", focus, "_diversity")), tree_stats)

    }

  }

  lm_tree_stats$focus <- focus
  lm_tree_stats$model$regulator <- tree_stats$regulator
  lm_tree_stats$model$module_size <- tree_stats$module_size

  lm_tree_stats

}
