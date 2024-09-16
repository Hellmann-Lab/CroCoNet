#' Plot module conservation
#'
#' Plots the tree-bases statistics that characterize the cross-species conservation of the modules and marks the modules that were found to be conserved/diverged.
#'
#' To determine whether a module as whole is conserved, diverged overall or diverged on a specific lineage, the CroCoNet approach relies on module trees reconstructed from pairwise preservation scores between clones and statistics calculated based on these trees (total tree length, species-to-other branch length, species-specific and overall diversity, for details please see \code{\link{calculatePresStats}}, \code{\link{reconstructTrees}}, \code{\link{calculateTreeStats}} and \code{\link{findConservedDivergedModules}}).
#'
#' The function plots each input module along the regression line that captures the relationship between total tree length and within-species diversity (in case the focus of interest is conservation and overall divergence) or between the species-to-other branch length and diversity of a species (in case the focus of interest is divergence between this particular species and all others). A module is considered diverged overall/in a species-specific manner if it falls above the prediction interval of the regression line, while a module is considered conserved if it falls below the prediction interval of the regression line.
#'
#' The regression line is drawn in dark grey, and the 95\% prediction interval of the line is shown as a light grey area. If the focus of interest is conservation and overall divergence, the conserved modules are colored green, the diverged modules are colored red, and the top 'N' conserved and diverged modules based on the variable specified in \code{rank_by} are labelled. If the focus of interest is divergence between a particular species and all others, only the diverged modules are colored and labelled, the modules below the prediction interval are not meaningful and thus not highlighted. If confidence intervals are provided for the variables, each module is plotted with error bars.
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
#' \item{within_species_diversity}{Numeric, within-species diveristy per module (typically the median across all jackknife versions of the module).}
#' \item{lwr_within_species_diversity}{Numeric, the lower bound of the confidence interval of the within-species diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{upr_within_species_diversity}{Numeric, the upper bound of the confidence interval of the within-species diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{fit}{Numeric, the fitted total tree length at the within-species diversity value of the module.}
#' \item{lwr_fit}{Numeric, the lower bound of the prediction interval of the fit.}
#' \item{upr_fit}{Numeric, the upper bound of the prediction interval of the fit.}
#' \item{residual}{Numeric, the residual of the module in the linear model. It is calculated as the difference between the observed and expected (fitted) total tree lengths.}
#' \item{t_score}{Numeric, the t-score of the module. It is calculated as the residual normalized by the standard error of the total tree length prediction at the given within-species diversity value}
#' \item{conservation}{Character, "not_significant" if the module falls inside the prediction interval of the fit, "diverged" if a module has a higher total tree length than the upper boundary of the prediction interval, and "conserved" if a module has a lower total tree length than the lower boundary of the prediction interval}
#' }
#' Columns required in case the focus of interest is species-specific divergence:
#' \describe{
#' \item{focus}{Character, the focus of interest in terms of cross-species conservation, the name of a species if the focus of interest is species-specific divergence}
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{\{\{species\}\}_to_other_branch_length}{Numeric, the species-to-other branch length per module (typically the median across all jackknife versions of the module).}
#' \item{lwr_\{\{species\}\}_to_other_branch_length}{Numeric, the lower bound of the confidence interval of the species-to-other branch length calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{upr_\{\{species\}\}_to_other_branch_length}{Numeric, the upper bound of the confidence interval of the species-to-other branch length calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{\{\{species\}\}_diversity}{Numeric, the diveristy of the species of interest per module (typically the median across all jackknife versions of the module).}
#' \item{lwr_\{\{species\}\}_diversity}{Numeric, the lower bound of the confidence interval of the species-specific diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{upr_\{\{species\}\}_diversity}{Numeric, the upper bound of the confidence interval of the species-specific diversity calculated based on the jackknifed versions of the module (optional, only needed for plotting error bars later on).}
#' \item{fit}{Numeric, the fiited species-to_other branch length at the species-specific diversity value of the module.}
#' \item{lwr_fit}{Numeric, the lower bound of the prediction interval of the fit.}
#' \item{upr_fit}{Numeric, the upper bound of the prediction interval of the fit.}
#' \item{residual}{Numeric, the residual of the module in the linear model. It is calculated as the difference between the observed and expected (fitted) species-to_other branch lengths.}
#' \item{t_score}{Numeric, the t-score of the module. It is calculated as the residual normalized by the standard error of the species-to_other branch length prediction at the given within-species diversity value.}
#' \item{conservation}{Character, 'not_significant' if the module falls inside the prediction interval of the fit, 'diverged' if a module has a higher species-to_other branch length than the upper boundary of the prediction interval, and 'conserved' if a module has a lower species-to_other branch length than the lower boundary of the prediction interval.}
#' }
#' @param N Integer, the number of top conserved and diverged modules to label.
#' @param rank_by Character, one of "residual" and "t_score". The name of the variable to rank the by when selecting the top N conserved and diverged modules to label.
#' @param colors (Named) character vector of length 2, the colors for the diverged and conserved modules.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A \code{\link{ggplot}} object.
#' @export
#'
#' @examples plotConservedDivergedModules(module_conservation_overall)
plotConservedDivergedModules <- function(module_conservation, N = 5L, rank_by = "residual", colors = NULL, font_size = 14) {

  if (!is.data.frame(module_conservation))
    stop("The argument \"module_conservation\" should be a data frame.")

  if (any(!c("focus", "regulator", "fit", "lwr_fit", "upr_fit", "residual", "t_score", "conservation") %in% colnames(module_conservation)))
    stop(paste0("The argument \"module_conservation\" should contain the columns \"focus\", \"regulator\", \"fit\", \"lwr_fit\", \"upr_fit\", \"residual\", \"t_score\" and \"conservation\"."))

  if (length(N) != 1 || (!inherits(N, "integer") && !(inherits(N, "numeric") && N == round(N))) || N < 0)
    stop("The argument \"N\" should be a positive integer or 0.")

  if (is.null(rank_by) || !rank_by %in% c("residual", "t_score"))
    stop("The argument \"rank_by\" should be one of \"residual\", \"t_score\".")

  if (!is.null(colors) && (!inherits(colors, "character") || any(!areColors(colors)) || length(colors) == 2))
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

    if (any(!(c(paste0(focus, "_diversity"), paste0(focus, "_to_other_branch_length")) %in% colnames(module_conservation))))
      stop(paste0("The argument \"module_conservation\" should contain the columns corresponding to the value of the column \"focus\". If \"focus\" is \"", focus, "\", the required columns are \"", focus, "_diversity\" and \"", focus, "_to_other_branch_length\"."))

  }

  if (is.null(colors)) {

    colors <- c(diverged = "#AA4139", conserved = "#2B823A", not_significant = "black")

  } else {

    if (is.null(names(colors))) names(colors) <- c("diverged", "conserved")
    colors <- c(colors, not_significant = "black")

  }

  if (focus == "overall") {

    p <- module_conservation %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[["within_species_diversity"]], y = .data[["total_tree_length"]])) +
      ggplot2::geom_line(ggplot2::aes(y = .data[["fit"]]), color = "grey30", linewidth = 0.5) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["lwr_fit"]], ymax = .data[["upr_fit"]]), fill = "grey80", alpha = 0.5) +
      ggplot2::geom_point(ggplot2::aes(color = .data[["conservation"]], size = .data[["conservation"]])) +
      ggplot2::theme_bw(base_size = font_size) +
      ggplot2::scale_color_manual(values = colors, breaks = c("diverged", "conserved")) +
      ggplot2::scale_size_manual(values = c(diverged = 0.5, conserved = 0.5, not_significant = 0.05), guide = "none") +
      ggplot2::xlab("within-species diversity") +
      ggplot2::ylab("total tree length") +
      ggplot2::theme(legend.title = ggplot2::element_blank())

    if (all(c("lwr_total_tree_length", "upr_total_tree_length", "lwr_within_species_diversity", "upr_within_species_diversity") %in% colnames(module_conservation))) {

      p <- p +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = .data[["lwr_within_species_diversity"]], xmax = .data[["upr_within_species_diversity"]], color = .data[["conservation"]]), linewidth = 0.2) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[["lwr_total_tree_length"]], ymax = .data[["upr_total_tree_length"]], color = .data[["conservation"]]), linewidth = 0.2)

    }

    p <- p +
      ggrepel::geom_label_repel(data = dplyr::bind_rows(module_conservation %>%
                                                          dplyr::filter(.data[["conservation"]] == "diverged") %>%
                                                          dplyr::slice_max(order_by = .data[[rank_by]], n = N),
                                                        module_conservation %>%
                                                          dplyr::filter(.data[["conservation"]] == "conserved") %>%
                                                          dplyr::slice_min(order_by = .data[[rank_by]], n = N)),
                                ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                                fill = "white", size = 2.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE)

  } else {

    p <- module_conservation %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[[paste0(focus, "_diversity")]], y = .data[[paste0(focus, "_to_other_branch_length")]])) +
      ggplot2::geom_line(ggplot2::aes(y = .data[["fit"]]), color = "grey30", linewidth = 0.5) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["lwr_fit"]], ymax = .data[["upr_fit"]]), fill = "grey80", alpha = 0.5) +
      ggplot2::geom_point(ggplot2::aes(color = .data[["conservation"]], size = .data[["conservation"]])) +
      ggplot2::theme_bw(base_size = font_size) +
      ggplot2::scale_color_manual(values = colors, breaks = "diverged", labels = paste0("diverged between\n", focus, " and other")) +
      ggplot2::scale_size_manual(values = c(diverged = 0.5, not_significant = 0.05), guide = "none") +
      ggplot2::xlab(paste0(focus, " diversity")) +
      ggplot2::ylab(paste0(focus, "-to-other branch length")) +
      ggplot2::theme(legend.title = ggplot2::element_blank())

    if (sum(c(paste0("lwr_", focus, "_to_other_branch_length"), paste0("upr_", focus, "_to_other_branch_length"), paste0("lwr_", focus, "_diversity"), paste0("upr_", focus, "_diversity")) %in% colnames(module_conservation)) == 4) {

      p <- p +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = .data[[paste0("lwr_", focus, "_diversity")]], xmax = .data[[paste0("upr_", focus, "_diversity")]], color = .data[["conservation"]]), linewidth = 0.2) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = .data[[paste0("lwr_", focus, "_to_other_branch_length")]], ymax = .data[[paste0("upr_", focus, "_to_other_branch_length")]], color = .data[["conservation"]]), linewidth = 0.2)

    }

    p <- p +
      ggrepel::geom_label_repel(data = module_conservation %>%
                                  dplyr::filter(.data[["conservation"]] == "diverged") %>%
                                  dplyr::slice_max(order_by = .data[[rank_by]], n = N),
                                ggplot2::aes(label = .data[["regulator"]], color = .data[["conservation"]]),
                                fill = "white", size = 2.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20, show.legend = FALSE)

  }

  p

}
