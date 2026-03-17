#' Plot module conservation
#'
#' Plots the tree-based statistics that characterize the cross-species conservation of the modules and marks the modules that were found to be conserved/diverged.
#'
#' To determine whether a module as a whole is conserved, diverged overall or diverged on a specific lineage, the CroCoNet approach relies on module trees reconstructed from pairwise preservation scores between replicates and statistics calculated based on these trees (total tree length, subtree lengths of the species, overall within-species and species-specific diversity, for details please see \code{\link{calculatePresStats}}, \code{\link{reconstructTrees}}, \code{\link{calculateTreeStats}} and \code{\link{findConservedDivergedModules}}).
#'
#' The function plots each input module along the regression line that captures the relationship between total tree length and within-species diversity (in case the focus of interest is conservation and overall divergence) or between the subtree length and diversity of a species (in case the focus of interest is divergence between this particular species and all others). A module is considered diverged overall/in a species-specific manner if it falls above the prediction interval of the regression line, while a module is considered conserved if it falls below the prediction interval of the regression line.
#'
#' The regression line is drawn in dark grey, and the 95\% prediction interval of the line is shown as a light grey area. If the focus of interest is conservation and overall divergence, the conserved modules are colored green, the diverged modules are colored red, and the top \code{N} conserved and diverged modules based on the variable specified in \code{rank_by} are labelled. If \code{robust_only = TRUE}, only modules that pass the robustness filter (\code{robust == TRUE}) are considered when selecting the top modules to label. If the focus of interest is divergence between a particular species and all others, only the diverged modules are colored and labelled. Modules below the prediction interval are not meaningful in this context and are therefore not highlighted.
#'
#' @param module_conservation Data frame of cross-species conservation measures per module.
#'
#' Columns required in case the focus of interest is conservation and overall divergence:
#' \describe{
#' \item{focus}{Character, the focus of interest in terms of cross-species conservation, "overall" if the focus of interest is conservation and overall divergence.}
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{total_tree_length}{Numeric, total tree length per module (typically the median across all jackknife versions of the module).}
#' \item{lwr_total_tree_length}{Numeric, the lower bound of the confidence interval of the total tree length calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{upr_total_tree_length}{Numeric, the upper bound of the confidence interval of the total tree length calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{within_species_diversity}{Numeric, within-species diversity per module (typically the median across all jackknife versions of the module).}
#' \item{lwr_within_species_diversity}{Numeric, the lower bound of the confidence interval of the within-species diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{upr_within_species_diversity}{Numeric, the upper bound of the confidence interval of the within-species diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{fit}{Numeric, the fitted total tree length at the within-species diversity value of the module.}
#' \item{lwr_fit}{Numeric, the lower bound of the prediction interval of the fit.}
#' \item{upr_fit}{Numeric, the upper bound of the prediction interval of the fit.}
#' \item{residual}{Numeric, the residual of the module in the linear model. It is calculated as the difference between the observed and expected (fitted) total tree lengths.}
#' \item{studentized_residual}{Numeric, the externally studentized residual of the module. This statistic measures how strongly a module deviates from the fitted regression line relative to the expected variability.}
#' \item{category}{Character, classification of the module based on the prediction interval: "within_expectation" if the module falls inside the prediction interval of the fit, "diverged" if a module has a higher total tree length than the upper boundary of the prediction interval, and "conserved" if a module has a lower total tree length than the lower boundary of the prediction interval (only for the overall analysis).}
#' \item{robust}{Logical, indicates whether the module passes an additional robustness filter based on the false discovery rate (FDR). Modules with an FDR below \code{fdr_cutoff} are marked as \code{TRUE}.}
#' }
#' Columns required in case the focus of interest is species-specific divergence:
#' \describe{
#' \item{focus}{Character, the focus of interest in terms of cross-species conservation, the name of a species if the focus of interest is species-specific divergence}
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{\{\{species\}\}_subtree_length}{Numeric, the subtree length of the species per module (typically the median across all jackknife versions of the module).}
#' \item{lwr_\{\{species\}\}_subtree_length}{Numeric, the lower bound of the confidence interval of the species subtree length calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{upr_\{\{species\}\}_subtree_length}{Numeric, the upper bound of the confidence interval of the species subtree length calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{\{\{species\}\}_diversity}{Numeric, the diversity of the species of interest per module (typically the median across all jackknife versions of the module).}
#' \item{lwr_\{\{species\}\}_diversity}{Numeric, the lower bound of the confidence interval of the species-specific diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{upr_\{\{species\}\}_diversity}{Numeric, the upper bound of the confidence interval of the species-specific diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{fit}{Numeric, the fitted subtree length of a species.}
#' \item{lwr_fit}{Numeric, the lower bound of the prediction interval of the fit.}
#' \item{upr_fit}{Numeric, the upper bound of the prediction interval of the fit.}
#' \item{residual}{Numeric, the residual of the module in the linear model. It is calculated as the difference between the observed and expected (fitted) subtree length of a species.}
#' \item{studentized_residual}{Numeric, the externally studentized residual of the module. This statistic measures how strongly a module deviates from the fitted regression line relative to the expected variability.}
#' \item{category}{Character, classification of the module based on the prediction interval: "diverged" if a module has a higher subtree length for the species of interest than the upper boundary of the prediction interval, "within_expectation" otherwise.}
#' \item{robust}{Logical, indicates whether the module passes an additional robustness filter based on the false discovery rate (FDR). Modules with an FDR below \code{fdr_cutoff} are marked as \code{TRUE}.}
#' }
#' @param N Integer, the number of top conserved and diverged modules to label.
#' @param rank_by Character, the name of the variable used to rank modules when selecting the top N conserved and diverged modules to label. Options are "residual" (default) and "studentized_residual".
#' @param robust_only Logical, if \code{TRUE}, the top conserved and diverged modules are selected only among modules that pass the robustness filter (\code{robust == TRUE}). If \code{FALSE} (default), the top modules are selected among all modules regardless of robustness.
#' @param colors (Named) character vector of length 2, the colors for the diverged and conserved modules.
#' @param font_size Numeric, font size (default: 14).
#' @param label_size Numeric, the size of the labels for the most conserved and diverged modules (default: 3.5).
#'
#' @return A \code{\link{ggplot}} object.
#' @export
#'
#' @examples plotConservedDivergedModules(module_conservation_overall)
plotConservedDivergedModules <- function(module_conservation, N = 5L, rank_by = "residual", robust_only = FALSE, colors = NULL, font_size = 14, label_size = 3.5) {

  if (!is.data.frame(module_conservation))
    stop("The argument \"module_conservation\" should be a data frame.")

  if (any(!c("focus", "regulator", "fit", "lwr_fit", "upr_fit", "residual", "studentized_residual", "category", "robust") %in% colnames(module_conservation)))
    stop(paste0("The argument \"module_conservation\" should contain the columns \"focus\", \"regulator\", \"fit\", \"lwr_fit\", \"upr_fit\", \"residual\", \"studentized_residual\",  \"category\" and \"robust\"."))

  if (length(N) != 1 || (!inherits(N, "integer") && !(inherits(N, "numeric") && N == round(N))) || N < 0)
    stop("The argument \"N\" should be a positive integer or 0.")

  if (is.null(rank_by) || !rank_by %in% c("residual", "studentized_residual"))
    stop("The argument \"rank_by\" should be one of \"residual\", \"studentized_residual\".")

  if (!inherits(robust_only, "logical") || length(robust_only) != 1)
    stop("The argument \"robust_only\" should be TRUE or FALSE.")

  if (!is.null(colors) && (!inherits(colors, "character") || any(!areColors(colors)) || length(colors) != 2))
    stop("The argument \"colors\" should be a character vector of valid color representations.")

  if (!is.null(names(colors)) && any(sort(names(colors)) != c("conserved", "diverged")))
    stop("The names of \"colors\" should be \"conserved\" and \"diverged\".")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  focus <- unique(module_conservation$focus)

  if (is.null(focus))
    stop("Please specify the focus of analysis by adding the column \"focus\" to \"module_conservation\". This should be \"overall\" if the focus of interest is conservation and overall divergence and the name of a species if the focus of interest is species-specific divergence.")

  if (!inherits(focus, "character") || length(focus) != 1)
    stop("The column \"focus\" of \"module_conservation\" should be \"overall\" if the focus of interest is conservation and overall divergence and the name of a species if the focus of interest is species-specific divergence.")

  if (focus == "overall") {

    if (any(!(c("within_species_diversity", "total_tree_length") %in% colnames(module_conservation))))
    stop("The argument \"module_conservation\" should contain the columns corresponding to the value of the column \"focus\". If \"focus\" is \"overall\", the required columns are \"within_species_diversity\" and \"total_tree_length\".")

  } else {

    if (any(!(c(paste0(focus, "_diversity"), paste0(focus, "_subtree_length")) %in% colnames(module_conservation))))
      stop(paste0("The argument \"module_conservation\" should contain the columns corresponding to the value of the column \"focus\". If \"focus\" is \"", focus, "\", the required columns are \"", focus, "_diversity\" and \"", focus, "_subtree_length\"."))

  }

  # Determine x/y variables
  if (focus == "overall") {
    x_var <- "within_species_diversity"
    y_var <- "total_tree_length"
    xlab <- "within-species diversity"
    ylab <- "total tree length"
  } else {
    x_var <- paste0(focus, "_diversity")
    y_var <- paste0(focus, "_subtree_length")
    xlab <- paste0(focus, " diversity")
    ylab <- paste0(focus, " subtree length")
  }
  x_lwr <- paste0("lwr_", x_var)
  x_upr <- paste0("upr_", x_var)
  y_lwr <- paste0("lwr_", y_var)
  y_upr <- paste0("upr_", y_var)

  if (is.null(colors)) {

    if (focus == "overall") {

      colors <- c(diverged = "#AA4139", conserved = "#2B823A", within_expectation = "black")

    } else {

      colors <- c(diverged = "salmon", within_expectation = "black")

    }

  } else {

    if (is.null(names(colors))) names(colors) <- c("diverged", "conserved")
    colors <- c(colors, within_expectation = "black")

  }

  # Dataset used for labeling
  label_data <- module_conservation

  if (robust_only)
    label_data <- dplyr::filter(label_data, .data$robust)

  # Select modules to label
  if (focus == "overall") {

    label_data <- dplyr::bind_rows(
      label_data %>%
        dplyr::filter(.data$category=="diverged") %>%
        dplyr::slice_max(order_by=.data[[rank_by]], n=N),
      label_data %>%
        dplyr::filter(.data$category=="conserved") %>%
        dplyr::slice_min(order_by=.data[[rank_by]], n=N)
    )

  } else {

    label_data <- label_data %>%
      dplyr::filter(.data$category=="diverged") %>%
      dplyr::slice_max(order_by=.data[[rank_by]], n=N)

  }

  # Base plot
  p <- module_conservation %>%
    ggplot2::ggplot(ggplot2::aes(x=.data[[x_var]], y=.data[[y_var]])) +
    ggplot2::geom_line(ggplot2::aes(y=.data$fit),
                       color="grey30", linewidth=0.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lwr_fit, ymax=.data$upr_fit),
                         fill="grey80", alpha=0.5) +
    ggplot2::geom_point(ggplot2::aes(color=.data$category,
                                     size=.data$category)) +
    ggplot2::theme_bw(base_size=font_size) +
    ggplot2::scale_color_manual(values=colors,
                                breaks=intersect(names(colors),
                                                 c("diverged","conserved"))) +
    ggplot2::scale_size_manual(values=c(diverged=0.5,
                                        conserved=0.5,
                                        within_expectation=0.05),
                               guide="none") +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::theme(legend.title=ggplot2::element_blank())

  if (all(c(x_lwr, x_upr, y_lwr, y_upr) %in% colnames(module_conservation))) {

    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data[[x_lwr]], xmax = .data[[x_upr]], color = .data$category),
        linewidth = 0.2
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[[y_lwr]], ymax = .data[[y_upr]], color = .data$category),
        linewidth = 0.2
      )

  }

  p +
    ggrepel::geom_label_repel(
      data = label_data,
      ggplot2::aes(label=.data$regulator, color=.data$category),
      fill="white",
      size=label_size,
      label.size=0.08,
      segment.size=0.08,
      box.padding=0.05,
      label.padding=0.05,
      force=2,
      max.overlaps=20,
      show.legend=FALSE
    )

}
