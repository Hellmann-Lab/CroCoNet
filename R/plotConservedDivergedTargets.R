#' Plot target conservation
#'
#' Plots the cross-species conservation of target genes within a conserved/diverged module and marks the targets that are particularly responsible for the conservation/divergence.
#'
#' To determine whether a module as whole is conserved, diverged overall or diverged on a specific lineage, the CroConet approach relies on module trees reconstrcuted from pairwise preservation scores between replicates and statistics calculated based on these trees (total tree length, overall within-species diversity, subtree length and diversity of each species). To identify individual genes that contribute the most to conservation/divergence, the same statistics can be used in combination with jackknifing. For more details on the approach and the generation of the \code{target_conservation} object please see \code{\link{findConservedDivergedTargets}}.
#'
#' The function plots each jackknife version of a module along the regression line that captures the relationship between total tree length and within-species diversity (in case the focus of interest is conservation and overall divergence) or between the subtree length and diversity of a species (in case the focus of interest is species-specific divergence). The residuals compared to this regression line characterize the cross-species conservation of the different jackknife module versions and the corresponding absent target genes: the jackknife version that has the highest residual is the most diverged, therefore the corresponding target gene is the most conserved, while the jackknife version that has the lowest residual is the most conserved, therefore the corresponding target gene is the most diverged.
#'
#' The regression line is drawn in dark grey, and the 95\% prediction interval of the line is shown as a light grey area. The jackknife module versions are colored by their residuals: the red to yellow part of the color spectrum denotes positive residuals, while the yellow to green part of the color spectrum denotes negative residuals (the entire spectrum always spans residuals from \code{-max(abs(residual))} to \code{max(abs(residual))} across all jackknife module versions of the module of interest).
#'
#' In case of a conserved module, the lightest green data points are typically the most interesting: these correspond to the smallest negative residuals, i.e. the least conserved jackknife module versions with the most conserved target genes removed. Analogously, in case of a diverged module, the lightest red data points are typically the most interesting: these correspond to the smallest positive residuals, i.e. the least diverged jackknife module versions with the most diverged target genes removed. In an extreme scenario, the least conserved jackknife module versions of a conserved module can reach positive residuals (red color) and the least diverged jackknife module versions of a diverged module can reach a negative residuals (green color). However, this should be rare if we assume that the divergence/conservation of a module cannot be contributed to a single target gene, but it is rather the combined signal of all targets together.
#'
#' As a default, in case of a conserved module, the top N least conserved jackknife module versions (most conserved targets) are labelled, while in case of a diverged module, the top N least diverged jackknife module versions (most diverged targets) are labelled. Using the parameter \code{label}, it is possible to manually choose whether the most conserved targets ("conserved"), the most diverged targets ("diverged"), or both the most conserved and most diverged targets ("both") should be labelled. The gene name in the label always refers to the target gene removed in the given jackknife version (format: "w/o \{\{geneName\}\}").
#'
#' In addition to the jackknife module versions corresponding to the most conserved/diverged target genes, the function always labels the original module version as well where no target gene was removed. The original module is always labelled as "original".
#'
#' @param target_conservation Data frame of cross-species conservation measures for each jackknife version of a module.
#'
#' Columns required in case the focus of interest is conservation and overall divergence:
#' \describe{
#' \item{focus}{Character, the focus of interest in terms of cross-species conservation, "overall" if the focus of interest is conservation and overall divergence.}
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{type}{Character, module type (orig = original, jk = jackknifed or summary = the summary across all jackknife module versions).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk', nameOfRegulator_orig in case of module type 'orig' and nameOfRegulator_summary in case of module type 'summary').}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module types 'orig' and 'summary').}
#' \item{total_tree_length}{Numeric, total tree length per jackknife module version.}
#' \item{within_species_diversity}{Numeric, within-species diversity per jackknife module version.}
#' \item{fit}{Numeric, the fitted total tree length at the within-species diversity value of a jackknife module version.}
#' \item{lwr_fit}{Numeric, the lower bound of the prediction interval of the fit.}
#' \item{upr_fit}{Numeric, the upper bound of the prediction interval of the fit.}
#' \item{residual}{Numeric, the residual of the jackknife module version in the linear model. It is calculated as the difference between the observed and expected (fitted) total tree lengths.}
#' }
#'
#' Columns required in case the focus of interest is species-specific divergence:
#' \describe{
#' \item{focus}{Character, the focus of interest in terms of cross-species conservation, the name of a species if the focus of interest is species-specific divergence.}
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{type}{Character, module type (orig = original, jk = jackknifed or summary = the summary across all jackknife module versions).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk', nameOfRegulator_orig in case of module type 'orig' and nameOfRegulator_summary in case of module type 'summary').}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module types 'orig' and 'summary').}
#' \item{\{\{species\}\}_subtree_length}{Numeric, the subtree length of a species per jackknife module version.}
#' \item{\{\{species\}\}_diversity}{Numeric, diversity within the species of interest per jackknife module version.}
#' \item{fit}{Numeric, the fitted subtree length of the species at the diversity value of a jackknife module version.}
#' \item{lwr_fit}{Numeric, the lower bound of the prediction interval of the fit.}
#' \item{upr_fit}{Numeric, the upper bound of the prediction interval of the fit.}
#' \item{residual}{Numeric, the residual of the jackknife module version in the linear model. It is calculated as the difference between the observed and expected (fitted) subtree length of the species.}
#' }
#' @param N Integer, the number of top conserved/diverged target genes to label.
#' @param label Character, specifies which data points should be labelled. If "auto" (default), in case of a conserved module the top N least conserved jackknife module versions (most conserved targets) are labelled, while in case of a diverged module the top N least diverged jackknife module versions (most diverged targets) are labelled. If "conserved", the top N most conserved targets are labelled, and if "diverged", the top N most diverged targets are labelled, regardless of the degree of conservation of the module as a whole. If "both", both the top N most conserved and the top N most diverged targets are labelled.
#' @param colors Character vector, the colors to visualize the residuals. The vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A \code{\link{ggplot}} object.
#' @export
#'
#' @examples plotConservedDivergedTargets(POU5F1_target_conservation)
plotConservedDivergedTargets <- function(target_conservation, N = 5, label = "auto", colors = NULL, font_size = 14) {

  if (!is.data.frame(target_conservation))
    stop("The argument \"target_conservation\" should be a data frame.")

  if (any(!c("focus", "regulator", "type", "id", "gene_removed", "fit", "lwr_fit", "upr_fit", "residual") %in% colnames(target_conservation)))
    stop("The argument \"target_conservation\" should contain the columns \"focus\", \"regulator\", \"type\", \"id\", \"gene_removed\", \"fit\", \"lwr_fit\", \"upr_fit\" and \"residual\".")

  if (length(N) != 1 || (!inherits(N, "integer") && !(inherits(N, "numeric") && N == round(N))) || N < 0)
    stop("The argument \"N\" should be a positive integer or 0.")

  if (is.null(label) || !all(label %in% c("auto", "diverged", "conserved", "both")))
    stop("The argument \"label\" should be one or more of \"auto\", \"diverged\", \"conserved\", \"both\".")

  if (label == "auto" && all(target_conservation$type != "summary"))
    stop("If the argument \"label\" is set to \"auto\", a row of type \"summary\" is expected to be present in \"target_conservation\". This row should contain the median/mean tree statistics across all jackknife version of the module which is needed to determin whether the module as a whole is conserved or diverged. Please add this information to the \"target_conservation\" object or set \"label\" to a different value.")

  if (!is.null(colors) && (!inherits(colors, "character") || any(!areColors(colors))))
    stop("The argument \"colors\" should be a character vector of valid color representations.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  focus <- unique(target_conservation$focus)

  if (is.null(focus))
    stop("Please specify the focus of analysis by adding the column \"focus\" to \"target_conservation\". This should be \"overall\" if the focus of interest is conservation and overall divergence and the name of a species if the focus of interest is species-specific divergence.")

  if (!inherits(focus, "character") || length(focus) != 1)
    stop("The column \"focus\" of \"target_conservation\" should be \"overall\" if the focus of interest is conservation and overall divergence and the name of a species if the focus of interest is species-specific divergence.")

  if (focus == "overall") {

    if (any(!(c("within_species_diversity", "total_tree_length") %in% colnames(target_conservation))))
      stop("The argument \"target_conservation\" should contain the columns corresponding to the value of the column \"focus\". If \"focus\" is \"overall\", the required columns are \"within_species_diversity\" and \"total_tree_length\".")

  } else {

    if (any(!(c(paste0(focus, "_diversity"), paste0(focus, "_subtree_length")) %in% colnames(target_conservation))))
      stop(paste0("The argument \"target_conservation\" should contain the columns corresponding to the value of the column \"focus\". If \"focus\" is \"", focus, "\", the required columns are \"", focus, "_diversity\" and \"", focus, "_subtree_length\"."))

  }

  if (is.null(colors))
    colors <- residual_colors

  if (label == "auto") {

    residual_summary <- target_conservation %>%
      dplyr::filter(.data[["type"]] == "summary") %>%
      dplyr::pull(.data[["residual"]])

    if (residual_summary > 0) label = "diverged" else label = "conserved"

  }

  target_conservation <- target_conservation %>%
    dplyr::filter(.data[["type"]] != "summary")

  if (label == "conserved") {

    target_conservation <- target_conservation %>%
      dplyr::mutate(to_label = ifelse(.data[["residual"]] %in% sort(.data[["residual"]], decreasing = TRUE)[1:N] | .data[["type"]] == "orig",
                                      TRUE, FALSE))

  } else if (label == "diverged") {

    target_conservation <- target_conservation %>%
      dplyr::mutate(to_label = ifelse(.data[["residual"]] %in% sort(.data[["residual"]])[1:N] | .data[["type"]] == "orig",
                                      TRUE, FALSE))

  } else {

    target_conservation <- target_conservation %>%
      dplyr::mutate(to_label = ifelse(.data[["residual"]] %in% sort(.data[["residual"]])[1:N] |
                                        .data[["residual"]] %in% sort(.data[["residual"]], decreasing = TRUE)[1:N] |
                                        .data[["type"]] == "orig",
                                      TRUE, FALSE))

  }

  max_residual <- max(abs(target_conservation$residual))

  if (focus == "overall") {

    target_conservation %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[["within_species_diversity"]], y = .data[["total_tree_length"]])) +
      ggplot2::geom_line(ggplot2::aes(y = .data[["fit"]]), color = "grey30", linewidth = 0.5) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["lwr_fit"]], ymax = .data[["upr_fit"]]), fill = "grey80", alpha = 0.5) +
      ggplot2::geom_point(ggplot2::aes(color = .data[["residual"]], size = .data[["to_label"]])) +
      ggplot2::theme_bw(base_size = font_size) +
      ggplot2:: scale_color_gradientn(colors = colors,
                                      limits = c(-max_residual, max_residual)) +
      ggplot2::scale_size_manual(values = c("TRUE" = 0.7, "FALSE" = 0.3), guide = "none") +
      ggplot2::ylab("total tree length") +
      ggplot2::xlab("within-species diversity") +
      ggrepel::geom_label_repel(data = target_conservation %>%
                                   dplyr::filter(.data[["to_label"]]) %>%
                                   dplyr::mutate(label = ifelse(.data[["type"]] == "orig", "original", paste0("w/o ", .data[["gene_removed"]]))),
                                 ggplot2::aes(label = .data[["label"]], color = .data[["residual"]]),
                                 fill = "white", size = 3, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20)


  } else {

    target_conservation %>%
      ggplot2::ggplot(ggplot2::aes(x = .data[[paste0(focus, "_diversity")]], y = .data[[paste0(focus, "_subtree_length")]])) +
      ggplot2::geom_line(ggplot2::aes(y = .data[["fit"]]), color = "grey30", linewidth = 0.5) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[["lwr_fit"]], ymax = .data[["upr_fit"]]), fill = "grey80", alpha = 0.5) +
      ggplot2::geom_point(ggplot2::aes(color = .data[["residual"]], size = .data[["to_label"]])) +
      ggplot2::theme_bw(base_size = font_size) +
      ggplot2:: scale_color_gradientn(colors = colors,
                                      limits = c(-max_residual, max_residual)) +
      ggplot2::scale_size_manual(values = c("TRUE" = 0.7, "FALSE" = 0.3), guide = "none") +
      ggplot2::xlab(paste0(focus, " diversity")) +
      ggplot2::ylab(paste0(focus, " subtree length")) +
      ggrepel::geom_label_repel(data = target_conservation %>%
                                  dplyr::filter(.data[["to_label"]]) %>%
                                  dplyr::mutate(label = ifelse(.data[["type"]] == "orig", "original", paste0("w/o ", .data[["gene_removed"]]))),
                                ggplot2::aes(label = .data[["label"]], color = .data[["residual"]]),
                                fill = "white", size = 3, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05, force = 2, max.overlaps = 20)

  }

}
