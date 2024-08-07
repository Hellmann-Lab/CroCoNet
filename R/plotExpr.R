#' Plot expression profiles as a heatmap
#'
#' Plots the expression profiles of one or more genes as a heatmap, with the cells ordered by pseudotime or cell type.
#'
#' The species, pseudotime and cell type information are taken from the metadata slot of the input \code{sce} object, and the expression data are taken from the logcounts assay of the input \code{sce} object. The function a heatmap, with columns corresponding to cells, rows corresponding to genes and colors corresponding to scaled and centered expression levels. The colors to represent the expression levels can be controlled by the argument \code{heatmap_colors}.
#'
#' Expression levels are clipped to the range of mean \eqn{\pm} 3\eqn{\sigma} per gene, this can be changed via the argument \code{clip}. Clipping aids visualization by preventing the outlier data points from squishing the rest of the data into a small color range.
#'
#' The cells can be order by pseudotime, cell type or any other variable that makes sense in the given context, the name of the column containing the desired cell property can be specified by the argument \code{order_by} (default: "pseudotime"). The colors to represent this cell property can be controlled by the argument \code{annotation_colors}.
#'
#' @param genes Character vector, the names of the genes for which the expression profiles should be plotted.
#' @param sce \code{\link{SingleCellExperiment-class}} object containing the expression data (raw counts, logcounts and metadata) for all network genes. Required metadata columns:
#'\describe{
#' \item{species}{Character, the name of the species.}
#' \item{\{\{order_by\}\}}{Numeric or character, the column by which the cells should be ordered on the plot. Typically contains the inferred pseudotime, or the cell type labels.}
#' }
#' @param order_by Character specifying the column in the metadata of \code{sce} by which the cells should be ordered on the plot. This column typically contains the inferred pseudotime (default: "pseudotime"), or the cell type labels.
#' @param heatmap_colors Character vector, the heatmap colors for the expression levels. The vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}.
#' @param annotation_colors Character vector, the colors for the variable specified by \code{order_by}. If the variable is discrete, the vector should contain as many colors as there are unique values of the variable, if the variable is continuous, the vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}.
#' @param clip Numeric specifying the degree of clipping. The expression levels of each gene that are more standard deviations away from the mean than \code{clip} are treated as NA. The Default is 3 meaning that expression levels are clipped to the range of mean \eqn{\pm} 3\eqn{\sigma} per gene.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A heatmap as a \code{\link{ggplot}} object showing the expression profiles of the input genes across all cells.
#' @export
#'
#' @examples plotExprHeatmap(c("POU5F1","DNMT3B","TERF1","TDGF1","L1TD1","VIM","MAP1B","MARCKS","PTN","CDH2"),
#'                           sce)
plotExprHeatmap <- function(genes, sce, order_by = "pseudotime", heatmap_colors = NULL, annotation_colors = NULL, clip = 3, font_size = 14) {

  if (!inherits(genes, "character"))
    stop("The argument \"genes\" should be a character vector.")

  if (!inherits(sce, "SingleCellExperiment"))
    stop("The argument \"sce\" should be of class \"SingleCellExperiment\".")

  if (!("logcounts" %in% names(SummarizedExperiment::assays(sce))))
    stop("The argument \"sce\" should contain the assay \"logcounts\".")

  if (is.null(sce$species))
    stop("The argument \"sce\" should contain the metadata column \"species\".")

  if (!inherits(order_by, "character") || length(order_by) != 1)
    stop("The argument \"order_by\" should be a string specifying the metadata column of \"sce\" by which the cells should be ordered on the plot.")

  if (is.null(sce[[order_by]]))
    stop("The argument \"sce\" should contain the metadata column specified by \"order_by\".")

  genes_not_found <- genes[!genes %in% rownames(sce)]

  if (length(genes_not_found) > 0)
    stop(paste0("The following input genes cannot be found in \"sce\": ", paste(genes_not_found, collapse = ", "), "."))

  if (!is.null(heatmap_colors) & (!inherits(heatmap_colors, "character") || any(!areColors(heatmap_colors))))
    stop("The argument \"heatmap_colors\" should be a character vector of valid color representations.")

  if (!is.null(annotation_colors) & (!inherits(annotation_colors, "character") || any(!areColors(annotation_colors))))
    stop("The argument \"annotation_colors\" should be a character vector of valid color representations.")

  if (!is.null(annotation_colors) && !inherits(sce[[order_by]], "numeric") && length(annotation_colors) != length(unique(sce[[order_by]])))
    stop("The argument \"annotation_colors\" should contain as many colors as there are unique values in the metadata column of \"sce\" specified by \"order_by\".")

  if (!inherits(clip, "numeric") || length(clip) != 1 || clip <= 0)
    stop("The argument \"clip\" should be a positive numeric value.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  n_genes <- length(genes)

  gene_name = NULL

  # extract metadata and logcounts of the genes of interest from the SCE object
  expr <- foreach::foreach(gene_name = genes,
                           .combine = dplyr::bind_rows,
                           .multicombine = TRUE) %do% {

                             data.frame(cell = colnames(sce),
                                        order_by = sce[[order_by]],
                                        gene = gene_name,
                                        expr = scale(SingleCellExperiment::logcounts(sce)[gene_name,]))

                           } %>%
    dplyr::mutate(gene = factor(.data[["gene"]], levels = rev(genes)))

  max_abs_expr <- min(max(abs(expr$expr)), clip)

  if (is.null(heatmap_colors))
    heatmap_colors <- eigengene_colors

  p <- expr %>%
    ggplot2::ggplot(ggplot2::aes(x = stats::reorder(.data[["cell"]], as.numeric(.data[["order_by"]])), y = .data[["gene"]], fill = .data[["expr"]])) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::scale_fill_gradientn(colors = heatmap_colors, limits = c(-max_abs_expr, max_abs_expr), guide = ggplot2::guide_colorbar(order = 1), name = "expression") +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(2, 2, 2, 2)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.85, 0))) +
    ggplot2::geom_rug(ggplot2::aes(color = .data[["order_by"]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(1 / (3 * n_genes), "npc")) +
    ggplot2::xlab("cells") +
    ggplot2::ylab("genes")

  if (is.numeric(expr$order_by)) {

    if (is.null(annotation_colors))
      annotation_colors <- pseudotime_colors

    return(
      p + ggplot2::scale_color_gradientn(colors = annotation_colors, name = order_by)
    )

  } else {

    if (is.null(annotation_colors))
      annotation_colors <- cell_type_color_ramp(seq(0, 1, length = length(unique(expr[[order_by]]))))

    p + ggplot2::scale_color_manual(values = annotation_colors, name = order_by) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 8.5))) +
      ggplot2::theme(legend.key.width = ggplot2::unit(0.035, "npc"))

  }

}


#' Plot module eigengenes as a heatmap
#'
#' Plots the eigengenes (or other types of summarized module expression profiles) as a heatmap, with the cells ordered by pseudotime or cell type.
#'
#' A concept adapted from WGCNA, the eigengene summarizes the expression profile of an entire module, and it is calculated as the first principal component of the module expression data (see also \code{\link{calculateEigengenes}}). Other possible ways of representing the expression profile of a module include the mean expression and the regulator expression.
#'
#' The function takes a data frame containing any of these summarized module expression profiles as input (normally the output of \code{\link{calculateEigengenes}}) and plots it as a heatmap, with columns corresponding to cells, rows corresponding to modules and colors corresponding to expression levels.
#'
#' The column of the data frame containing the module expression profiles can be specified by the argument \code{expr_column} (default: "eigengene"). The colors to represent the expression levels can be controlled by the argument \code{heatmap_colors}.
#'
#' The cells can be order by pseudotime, cell type or any other variable that makes sense in the given context, the name of the column containing the desired cell property can be specified by the argument \code{order_by} (default: "pseudotime"). The colors to represent this cell property can be controlled by the argument \code{annotation_colors}.
#'
#' If the eigengene (or mean expression) has been calculated for positively and negatively regulated targets separately, then these will appear as separate rows in the heatmap.
#'
#' @param eigengenes Data frame of eigengenes, required columns:
#'\describe{
#' \item{cell}{Character, the cell barcode.}
#' \item{\{\{order_by\}\}}{Numeric or character, the column by which the cells should be ordered on the plot. Typically contains the inferred pseudotime, or the cell type labels.}
#' \item{module}{Character, transcriptional regulator and in case the eigengene was calculated for the positively or negatively regulated targets only, the direction of regulation (format: nameOfRegulator(+) or nameOfRegulator(-)).}
#' \item{\{\{expr_column\}\}}{Numeric, summarized module expression profiles (typically the eigengene, the mean expression of the module, or the expression of the regulator).}
#' }
#' @param expr_column Character specifying the column of \code{eigengenes} by which the heatmap should be colored. This column is expected to contain summarized module expression profiles, typically the eigengene (default: "eigengene"), the mean expression of the module, or the expression of the regulator.
#' @param order_by Character specifying which column of \code{eigengenes} by which the cells should be ordered on the plot. This column typically contains the inferred pseudotime (default: "pseudotime"), or the cell type labels.
#' @param heatmap_colors Character vector, the heatmap colors for the variable specified by \code{expr_column}. The vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}.
#' @param annotation_colors Character vector, the colors for the variable specified by \code{order_by}. If the variable is discrete, the vector should contain as many colors as there are unique values of the variable, if the variable is continuous, the vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{scale_color_gradientn}.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A heatmap as a \code{\link{ggplot}} object showing the summarized expression profiles of the modules across all cells.
#' @export
#'
#' @examples plotEigengeneHeatmap(eigengenes)
plotEigengeneHeatmap <- function(eigengenes, expr_column = "eigengene", order_by = "pseudotime", heatmap_colors = NULL, annotation_colors = NULL, font_size = 14) {

  if (!is.data.frame(eigengenes))
    stop("The argument \"eigengenes\" should be a data frame.")

  if (any(!c("cell", "module") %in% colnames(eigengenes)))
    stop("The argument \"eigengenes\" should contain the columns \"cell\" and \"module\".")

  if (!inherits(expr_column, "character") || length(expr_column) != 1 || !expr_column %in% colnames(eigengenes))
    stop("The argument \"expr_column\" should be a string specifying the column of \"eigengenes\" that contains the summarized module expression profiles to be plotted.")

  if (!inherits(order_by, "character") || length(order_by) != 1 || !order_by %in% colnames(eigengenes))
    stop("The argument \"order_by\" should be a string specifying the column of \"eigengenes\" by which the cells should be ordered on the plot.")

  if (!is.null(heatmap_colors) && (!inherits(heatmap_colors, "character") || any(!areColors(heatmap_colors))))
    stop("The argument \"heatmap_colors\" should be a character vector of valid color representations.")

  if (!is.null(annotation_colors) && (!inherits(annotation_colors, "character") || any(!areColors(annotation_colors))))
    stop("The argument \"annotation_colors\" should be a character vector of valid color representations.")

  if (!is.null(annotation_colors) && !inherits(eigengenes[[order_by]], "numeric") && length(annotation_colors) != length(unique(eigengenes[[order_by]])))
    stop("The argument \"annotation_colors\" should contain as many colors as there are unique values in the column of \"eigengenes\" specified by \"order_by\".")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  n_modules <- length(unique(eigengenes$module))

  if (is.null(heatmap_colors))
    heatmap_colors <- eigengene_colors

  p <- eigengenes %>%
    ggplot2::ggplot(ggplot2::aes(x = stats::reorder(.data[["cell"]], as.numeric(.data[[order_by]])), y = .data[["module"]], fill = .data[[expr_column]])) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::scale_fill_gradientn(colors = heatmap_colors, guide = ggplot2::guide_colorbar(order = 1)) +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(2, 2, 2, 2)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.85, 0)), limits = rev) +
    ggplot2::geom_rug(ggplot2::aes(color = .data[[order_by]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(1 / (3 * n_modules), "npc")) +
    ggplot2::xlab("cells") +
    ggplot2::ylab("modules")

  if (is.numeric(eigengenes[[order_by]])) {

    if (is.null(annotation_colors))
      annotation_colors <- pseudotime_colors

    return(
      p + ggplot2::scale_color_gradientn(colors = annotation_colors)
    )

  } else {

    if (is.null(annotation_colors))
      annotation_colors <- cell_type_color_ramp(seq(0, 1, length = length(unique(eigengenes[[order_by]]))))

    p + ggplot2::scale_color_manual(values = annotation_colors) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 8.5))) +
      ggplot2::theme(legend.key.width = ggplot2::unit(0.035, "npc"))

  }

}


#' Plot expression profiles along a pseudotime trajectory
#'
#' Plots the expression profiles of one or more genes along a pseudotime trajectory per species, and thus allows the expression patterns to be visually compared across species.
#'
#' The species, pseudotime and cell type information are taken from the metadata slot of the input \code{sce} object, and the expression data are taken from the logcounts assay of the input \code{sce} object. The smoothed expression profiles are fitted per species and gene using \code{\link{loess}} with the formula = "expression ~ pseusotime" and plotted as continuous lines colored by species and faceted by gene. The 95\% confidence intervals of the fitted lines are calculated using \code{\link{predict}} and shown as lightly colored areas around the lines. If a cell type column is specified by the argument \code{cell_type_column}, the cell are shown as a rug plot along the pseudotime axis colored by cell types.
#'
#' @param genes Character vector, the names of the genes for which the expression profiles should be plotted.
#' @param sce \code{\link{SingleCellExperiment-class}} object containing the expression data (raw counts, logcounts and metadata) for all network genes. Required metadata columns:
#'\describe{
#' \item{species}{Character, the name of the species.}
#' \item{\{\{pseudotime_column\}\}}{Numeric, inferred pseudotime.}
#' \item{\{\{cell_type_column\}\}}{Character, cell type annotation (optional).}
#' }
#' @param pseudotime_column Character, the name of the pseudotime column in the metadata of \code{sce} (default: "pseudotime").
#' @param cell_type_column Character, the name of the cell type annotation column in the metadata of \code{sce} (default: "cell_type", if there is no cell type annotation available or the user wants to omit the cell type rug plot, this argument should be set to NULL).
#' @param species_colors Character vector, colors per species.
#' @param cell_type_colors Character vector, colors per cell type.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A \code{\link{ggplot}} object.
#' @export
#'
#' @examples plotExprAlongPseudotime(regulators, sce)
plotExprAlongPseudotime <- function(genes, sce, pseudotime_column = "pseudotime", cell_type_column = "cell_type", species_colors = NULL, cell_type_colors = NULL, font_size = 14) {

  if (!inherits(genes, "character"))
    stop("The argument \"genes\" should be a character vector.")

  if (!inherits(sce, "SingleCellExperiment"))
    stop("The argument \"sce\" should be of class \"SingleCellExperiment\".")

  if (!("logcounts" %in% names(SummarizedExperiment::assays(sce))))
    stop("The argument \"sce\" should contain the assay \"logcounts\".")

  if (is.null(sce$species))
    stop("The argument \"sce\" should contain the metadata column \"species\".")

  if (!inherits(pseudotime_column, "character") || length(pseudotime_column) != 1 || is.null(sce[[pseudotime_column]]))
    stop("The argument \"pseudotime_column\" should be a string specifying the metadata column of \"sce\" that contains the pseudotime values.")

  if (!is.null(cell_type_column) && (!inherits(cell_type_column, "character") || length(cell_type_column) != 1 || is.null(sce[[cell_type_column]])))
    stop("The argument \"cell_type_column\" should be a string specifying the metadata column of \"sce\" that contains the cell type labels.")

  genes_not_found <- genes[!genes %in% rownames(sce)]

  if (length(genes_not_found) > 0)
    stop(paste0("The following input genes cannot be found in \"sce\": ", paste(genes_not_found, collapse = ", "), "."))

  if (!is.null(species_colors) && (!inherits(species_colors, "character") || any(!areColors(species_colors))))
    stop("The argument \"species_colors\" should be a character vector of valid color representations.")

  if (!is.null(cell_type_colors) && (!inherits(cell_type_colors, "character") || any(!areColors(cell_type_colors))))
    stop("The argument \"cell_type_colors\" should be a character vector of valid color representations.")

  if (!is.null(species_colors) && length(species_colors) != length(unique(sce$species)))
    stop("The argument \"species_colors\" should contain as many colors as there are unique species in the \"species\" metadata column of \"sce\".")

  if (!is.null(cell_type_column) && !is.null(cell_type_colors) && length(cell_type_colors) != length(unique(sce[[cell_type_column]])))
    stop("The argument \"cell_type_colors\" should contain as many colors as there are unique cell types in the metadata column of \"sce\" specified by \"cell_type_column\".")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  gene_name = species_name = NULL

  # grab default colors if the user did not input any custom colors
  species_names <- unique(sce$species)
  if (is.null(species_colors))
    species_colors <- species_color_ramp(seq(0, 1, len = length(species_names)))

  # extract metadata and logcounts of the genes of interest from the SCE object
  expr <- foreach::foreach(gene_name = genes,
                           .combine = dplyr::bind_rows,
                           .multicombine = TRUE) %do% {

                             expr_gene <- data.frame(species = sce$species,
                                                     pseudotime = sce[[pseudotime_column]],
                                                     gene = gene_name,
                                                     expr = SingleCellExperiment::logcounts(sce)[gene_name,])

                             if (!is.null(cell_type_column)) {

                               expr_gene$cell_type <- sce[[cell_type_column]]

                             }

                             expr_gene

                           } %>%
    dplyr::mutate(gene = factor(.data[["gene"]], levels = genes))

  # loess fit per gene and species
  loess_fit <- foreach::foreach(gene_name = genes,
                                .combine = dplyr::bind_rows,
                                .multicombine = TRUE) %:%
               foreach::foreach(species_name = species_names,
                                .combine = dplyr::bind_rows,
                                .multicombine = TRUE)  %do% {

                       expr_gene_spec <- expr %>%
                         dplyr::filter(.data[["gene"]] == gene_name & .data[["species"]] == species_name)
                       fit  = stats::loess(expr ~ pseudotime, data = expr_gene_spec)
                       newx = data.frame(pseudotime = seq(min(expr_gene_spec$pseudotime), max(expr_gene_spec$pseudotime), len = 80))
                       pred = stats::predict(fit, newdata = newx, se = TRUE)
                       cbind(pseudotime = newx, expr = pred$fit, se = pred$se.fit) %>%
                         dplyr::mutate(moe = stats::qt(p = 0.975, df = nrow(expr_gene_spec) - 2)*.data[["se"]],
                                       lwr = .data[["expr"]] - .data[["moe"]],
                                       upr = .data[["expr"]] + .data[["moe"]],
                                       gene = gene_name,
                                       species = species_name)

                     } %>%
    dplyr::mutate(gene = factor(.data[["gene"]], levels = genes))

  # rug positions and colors if cell type information is present
  if (!is.null(cell_type_column)) {

    if (is.null(cell_type_colors))
      cell_type_colors <- cell_type_color_ramp(seq(0, 1, len = length(unique(sce[[cell_type_column]]))))

    rug_positions <- loess_fit %>%
      dplyr::group_by(.data[["gene"]]) %>%
      dplyr::summarize(pos = min(.data[["lwr"]]) - 0.1*(max(.data[["upr"]]) - min(.data[["lwr"]]))) %>%
      tibble::deframe()

  }

  # actual plot
  p <- expr %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["pseudotime"]], y = .data[["expr"]])) +
    ggplot2::xlab("pseudotime") +
    ggplot2::ylab("expression (logcounts)") +
    ggplot2::theme_bw(base_size = font_size)

  if (!is.null(cell_type_column)) {

    p <- p +
      ggplot2::geom_rug(data = expr %>% dplyr::mutate(expr = rug_positions[.data[["gene"]]]),
                        ggplot2::aes(color = .data[["cell_type"]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(0.1, "npc"), linewidth = 0.15) +
      ggplot2::scale_color_manual(values = cell_type_colors,  name = "cell type") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
      ggnewscale::new_scale_color()

  }

  p +
    ggplot2::geom_line(data = loess_fit, ggplot2::aes(color = .data[["species"]], group = .data[["species"]]), linewidth = 1) +
    ggplot2::geom_ribbon(data = loess_fit, ggplot2::aes(ymin = .data[["lwr"]], ymax = .data[["upr"]], fill = .data[["species"]]), alpha = 0.1) +
    ggplot2::scale_color_manual(values = species_colors) +
    ggplot2::scale_fill_manual(values = species_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
    ggplot2::facet_wrap(~.data[["gene"]], scales = "free_y", ncol = 1, strip.position = "left") +
    ggplot2::scale_y_continuous(position = "right")

}


#' Plot module eigengenes along a pseudotime trajectory
#'
#' Plots the eigengenes (or other types of summarized module expression profiles) along a pseudotime trajectory per species, and thus allows the expression patterns to be visually compared across species.
#'
#' A concept adapted from WGCNA, the eigengene summarizes the expression profile of an entire module, and it is calculated as the first principal component of the module expression data (see also \code{\link{calculateEigengenes}}). Other possible ways of representing the expression profile of a module include the mean expression and the regulator expression.
#'
#' The function takes a data frame containing any of these summarized module expression profiles as input (normally the output of \code{\link{calculateEigengenes}}). The smoothed expression profiles are fitted per species and module using \code{\link{loess}} with the formula = "expression ~ pseusotime" and plotted as continuous lines colored by species and faceted by module. The 95\% confidence intervals of the fitted lines are calculated using \code{\link{predict}} and shown as lightly colored areas around the lines. If a cell type column is specified by the argument \code{cell_type_column}, the cell are shown as a rug plot along the pseudotime axis colored by cell types.
#'
#' If the eigengene (or mean expression) has been calculated for positively and negatively regulated targets separately, then these will appear as separate facets in the plot.
#'
#' @param eigengenes Data frame of eigengenes, required columns:
#'\describe{
#' \item{cell}{Character, the cell barcode.}
#' \item{species}{Character, the name of the species.}
#' \item{\{\{pseudotime_column\}\}}{Numeric, inferred pseudotime.}
#' \item{\{\{cell_type_column\}\}}{Character, cell type annotation (optional).}
#' \item{module}{Character, transcriptional regulator and in case the eigengene was calculated for the positively or negatively regulated targets only, the direction of regulation (format: nameOfRegulator(+) or nameOfRegulator(-)).}
#' \item{\{\{expr_column\}\}}{Numeric, summarized module expression profiles (typically the eigengene, the mean expression of the module, or the expression of the regulator).}
#' }
#' @param expr_column Character specifying which column of \code{eigengenes} the heatmap should be colored by. This column is expected to contain summarized module expression profiles, typically the eigengene (default: "eigengene"), the mean expression of the module, or the expression of the regulator.
#' @param pseudotime_column Character, the name of the pseudotime column in \code{eigengenes} (default: "pseudotime").
#' @param cell_type_column Character, the name of the cell type annotation column in \code{eigengenes} (default: "cell_type", if there is no cell type annotation available or the user wants to omit the cell type rug plot, this argument should be set to NULL).
#' @param species_colors Character vector, colors per species.
#' @param cell_type_colors Character vector, colors per cell type.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A \code{\link{ggplot}} object.
#' @export
#'
#' @examples
#' plotEigengenesAlongPseudotime(eigengenes_per_species)
#' plotEigengenesAlongPseudotime(eigengenes_per_species, expr_column = "mean_expr")
#' plotEigengenesAlongPseudotime(eigengenes_per_species, cell_type_column = NULL)
plotEigengenesAlongPseudotime <- function(eigengenes, expr_column = "eigengene",  pseudotime_column = "pseudotime", cell_type_column = "cell_type", species_colors = NULL, cell_type_colors = NULL, font_size = 14) {

  if (!is.data.frame(eigengenes))
    stop("The argument \"eigengenes\" should be a data frame.")

  if (any(!c("cell", "species", "module") %in% colnames(eigengenes)))
    stop("The argument \"eigengenes\" should contain the columns \"cell\", \"species\" and \"module\".")

  if (!inherits(expr_column, "character") || length(expr_column) != 1 || !expr_column %in% colnames(eigengenes))
    stop("The argument \"expr_column\" should be a string specifying the column of \"eigengenes\" that contains the summarized module expression profiles to be plotted.")

  if (!inherits(pseudotime_column, "character") || length(pseudotime_column) != 1 || !pseudotime_column %in% colnames(eigengenes))
    stop("The argument \"pseudotime_column\" should be a string specifying the column of \"eigengenes\" that contains the pseudotime values.")

  if (!is.null(cell_type_column) && (!inherits(cell_type_column, "character") || length(cell_type_column) != 1 || !cell_type_column %in% colnames(eigengenes)))
    stop("The argument \"cell_type_column\" should be a string specifying the metadata column of \"eigengenes\" that contains the cell type labels.")

  if (!is.null(species_colors) && (!inherits(species_colors, "character") || any(!areColors(species_colors))))
    stop("The argument \"species_colors\" should be a character vector of valid color representations.")

  if (!is.null(cell_type_colors) && (!inherits(cell_type_colors, "character") || any(!areColors(cell_type_colors))))
    stop("The argument \"cell_type_colors\" should be a character vector of valid color representations.")

  if (!is.null(species_colors) && length(species_colors) != length(unique(eigengenes$species)))
    stop("The argument \"species_colors\" should contain as many colors as there are unique species in the \"species\" column of \"eigengenes\".")

  if (!is.null(cell_type_column) && !is.null(cell_type_colors) && length(cell_type_colors) != length(unique(eigengenes[[cell_type_column]])))
    stop("The argument \"cell_type_colors\" should contain as many colors as there are unique cell types in the column of \"eigengenes\" specified by \"cell_type_column\".")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  module_name = species_name = NULL

  # get all modules and species
  eigengenes <- eigengenes %>%
    dplyr::rename(expr = expr_column,
                  pseudotime = pseudotime_column)
  module_names <- unique(eigengenes$module)
  species_names <- unique(eigengenes$species)

  # grab default colors if the user did not input any custom colors
  if (is.null(species_colors))
    species_colors <- species_color_ramp(seq(0, 1, len = length(species_names)))

  # loess fit per module and species
  loess_fit <- foreach::foreach(module_name = module_names,
                                .combine = dplyr::bind_rows,
                                .multicombine = TRUE) %:%
               foreach::foreach(species_name = species_names,
                                .combine = dplyr::bind_rows,
                                .multicombine = TRUE)  %do% {

                       eigengene_mod_spec <- eigengenes %>%
                         dplyr::filter(.data[["module"]] == module_name & .data[["species"]] == species_name)
                       fit  = stats::loess(expr ~ pseudotime, data = eigengene_mod_spec)
                       newx = data.frame(pseudotime = seq(min(eigengene_mod_spec$pseudotime), max(eigengene_mod_spec$pseudotime), len = 80))
                       pred = stats::predict(fit, newdata = newx, se = TRUE)
                       cbind(pseudotime = newx, expr = pred$fit, se = pred$se.fit) %>%
                         dplyr::mutate(moe = stats::qt(p = 0.975, df = nrow(eigengene_mod_spec) - 2)*.data[["se"]],
                                       lwr = .data[["expr"]] - .data[["moe"]],
                                       upr = .data[["expr"]] + .data[["moe"]],
                                       module = module_name,
                                       species = species_name)

                     } %>%
    dplyr::mutate(module = factor(.data[["module"]], levels = module_names))

  # rug positions and colors if cell type information is present
  if (!is.null(cell_type_column)) {

    if (is.null(cell_type_colors))
      cell_type_colors <- cell_type_color_ramp(seq(0, 1, len = length(unique(eigengenes[[cell_type_column]]))))

    rug_positions <- loess_fit %>%
      dplyr::group_by(.data[["module"]]) %>%
      dplyr::summarize(pos = min(.data[["lwr"]]) - 0.1*(max(.data[["upr"]]) - min(.data[["lwr"]]))) %>%
      tibble::deframe()

  }

  # actual plot
  p <- eigengenes %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["pseudotime"]], y = .data[["expr"]])) +
    ggplot2::xlab("pseudotime") +
    ggplot2::ylab(expr_column) +
    ggplot2::theme_bw(base_size = font_size)

  if (!is.null(cell_type_column)) {

    p <- p +
      ggplot2::geom_rug(data = eigengenes %>% dplyr::mutate(expr = rug_positions[.data[["module"]]]),
                        ggplot2::aes(color = .data[[cell_type_column]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(0.1, "npc"), linewidth = 0.15) +
      ggplot2::scale_color_manual(values = cell_type_colors,  name = "cell type") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
      ggnewscale::new_scale_color()

  }

  p +
    ggplot2::geom_line(data = loess_fit, ggplot2::aes(color = .data[["species"]], group = .data[["species"]]), linewidth = 1) +
    ggplot2::geom_ribbon(data = loess_fit, ggplot2::aes(ymin = .data[["lwr"]], ymax = .data[["upr"]], fill = .data[["species"]]), alpha = 0.1) +
    ggplot2::scale_color_manual(values = species_colors) +
    ggplot2::scale_fill_manual(values = species_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
    ggplot2::facet_wrap(~.data[["module"]], scales = "free_y", ncol = 1, strip.position = "left") +
    ggplot2::scale_y_continuous(position = "right")

}


#' Plot expression levels per cell type
#'
#' Plots the expression profiles of one or more genes per cell type and species, and thus allows the expression patterns to be visually compared across species.
#'
#' The species and cell type information are taken from the metadata slot of the input \code{sce} object, and the expression data are taken from the logcounts assay of the input \code{sce} object.
#'
#' @param genes Character vector, the names of the genes for which the expression profiles should be plotted.
#' @param sce \code{\link{SingleCellExperiment-class}} object containing the expression data (raw counts, logcounts and metadata) for all network genes. Required metadata columns:
#'\describe{
#' \item{species}{Character, the name of the species.}
#' \item{\{\{cell_type_column\}\}}{Character, cell type annotation.}
#' }
#' @param cell_type_column Character, the name of the cell type annotation column in the metadata of \code{sce}.
#' @param species_colors Character vector, colors per species.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A \code{\link{ggplot}} object.
#' @export
#'
#' @examples plotExprPerCellType(regulators, sce)
plotExprPerCellType <- function(genes, sce, cell_type_column = "cell_type", species_colors = NULL, font_size = 14) {

  if (!inherits(genes, "character"))
    stop("The argument \"genes\" should be a character vector.")

  if (!inherits(sce, "SingleCellExperiment"))
    stop("The argument \"sce\" should be of class \"SingleCellExperiment\".")

  if (!("logcounts" %in% names(SummarizedExperiment::assays(sce))))
    stop("The argument \"sce\" should contain the assay \"logcounts\".")

  if (is.null(sce$species))
    stop("The argument \"sce\" should contain the metadata column \"species\".")

  if (!inherits(cell_type_column, "character") || length(cell_type_column) != 1)
    stop("The argument \"cell_type_column\" should be a string specifying the metadata column of \"sce\" that contains the cell type labels.")

  if (is.null(sce[[cell_type_column]]))
    stop("The argument \"sce\" should contain the metadata column specified by \"cell_type_column\".")

  genes_not_found <- genes[!genes %in% rownames(sce)]

  if (length(genes_not_found) > 0)
    stop(paste0("The following input genes cannot be found in \"sce\": ", paste(genes_not_found, collapse = ", "), "."))

  if (!is.null(species_colors) && (!inherits(species_colors, "character") || any(!areColors(species_colors))))
    stop("The argument \"species_colors\" should be a character vector of valid color representations.")

  if (!is.null(species_colors) && length(species_colors) != length(unique(sce$species)))
    stop("The argument \"species_colors\" should contain as many colors as there are unique species in the \"species\" metadata column of \"sce\".")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  gene = NULL

  if (is.null(species_colors))
    species_colors <- species_color_ramp(seq(0, 1, length = dplyr::n_distinct(sce$species)))

  expr <- foreach::foreach(gene = genes,
                           .combine = dplyr::bind_rows,
                           .multicombine = TRUE) %do% {

                             data.frame(species = sce$species,
                                        gene = gene,
                                        expr = SingleCellExperiment::logcounts(sce)[gene,],
                                        cell_type = sce[[cell_type_column]])

                           } %>%
    dplyr::mutate(gene = factor(.data[["gene"]], levels = genes))

  expr %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["cell_type"]], y = .data[["expr"]], fill = .data[["species"]])) +
    ggplot2::geom_violin(draw_quantiles = 0.5) +
    ggplot2::scale_fill_manual(values = species_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
    ggplot2::xlab("cell type") +
    ggplot2::ylab("expression (logcounts)") +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::facet_wrap(~.data[["gene"]], scales = "free_y", ncol = 1, strip.position = "left") +
    ggplot2::scale_y_continuous(position = "right")

}


#' Plot eigengenes per cell type
#'
#' Plots the eigengenes (or other types of summarized module expression profiles) per cell type and species, and thus allows the expression patterns to be visually compared across species.
#'
#' A concept adapted from WGCNA, the eigengene summarizes the expression profile of an entire module, and it is calculated as the first principal component of the module expression data (see also \code{\link{calculateEigengenes}}). Other possible ways of representing the expression profile of a module include the mean expression and the regulator expression.
#'
#' @param eigengenes Data frame of eigengenes, required columns:
#'\describe{
#' \item{cell}{Character, the cell barcode.}
#' \item{species}{Character, the name of the species.}
#' \item{\{\{cell_type_column\}\}}{Character, cell type annotation.}
#' \item{module}{Character, transcriptional regulator and in case the eigengene was calculated for the positively or negatively regulated targets only, the direction of regulation (format: nameOfRegulator(+) or nameOfRegulator(-)).}
#' \item{\{\{expr_column\}\}}{Numeric, summarized module expression profiles (typically the eigengene, the mean expression of the module, or the expression of the regulator).}
#' }
#' @param expr_column Character specifying which column of \code{eigengenes} the heatmap should be colored by. This column is expected to contain summarized module expression profiles, typically the eigengene (default: "eigengene"), the mean expression of the module, or the expression of the regulator.
#' @param cell_type_column Character, the name of the cell type annotation column in \code{eigengenes}.
#' @param species_colors Character vector, colors per species.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A \code{\link{ggplot}} object.
#' @export
#'
#' @examples plotEigengenesPerCellType(eigengenes)
plotEigengenesPerCellType <- function(eigengenes, expr_column = "eigengene", cell_type_column = "cell_type", species_colors = NULL, font_size = 14) {

  if (!is.data.frame(eigengenes))
    stop("The argument \"eigengenes\" should be a data frame.")

  if (any(!c("cell", "species", "module") %in% colnames(eigengenes)))
    stop("The argument \"eigengenes\" should contain the columns \"cell\", \"species\" and \"module\".")

  if (!inherits(expr_column, "character") || length(expr_column) != 1 || !expr_column %in% colnames(eigengenes))
    stop("The argument \"expr_column\" should be a string specifying the column of \"eigengenes\" that contains the summarized module expression profiles to be plotted.")

  if (!inherits(cell_type_column, "character") || length(cell_type_column) != 1 || !cell_type_column %in% colnames(eigengenes))
    stop("The argument \"cell_type_column\" should be a string specifying the metadata column of \"eigengenes\" that contains the cell type labels.")

  if (!is.null(species_colors) && (!inherits(species_colors, "character") || any(!areColors(species_colors))))
    stop("The argument \"species_colors\" should be a character vector of valid color representations.")

  if (!is.null(species_colors) && length(species_colors) != length(unique(eigengenes$species)))
    stop("The argument \"species_colors\" should contain as many colors as there are unique species in the \"species\" column of \"eigengenes\".")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  if (is.null(species_colors))
    species_colors <- species_color_ramp(seq(0, 1, len = length(unique(eigengenes$species))))

  eigengenes %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["cell_type"]], y = .data[[expr_column]], fill = .data[["species"]])) +
    ggplot2::geom_violin(draw_quantiles = 0.5) +
    ggplot2::scale_fill_manual(values = species_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
    ggplot2::xlab("cell type") +
    ggplot2::ylab(expr_column) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::facet_wrap(~.data[["module"]], scales = "free_y", ncol = 1, strip.position = "left") +
    ggplot2::scale_y_continuous(position = "right")

}
