#' Plot expression profiles as a heatmap
#'
#' Plots the expression profiles of one or more genes as a heatmap, with the cells ordered by pseudotime or cell type.
#'
#' The function creates a heatmap, with columns corresponding to cells, rows corresponding to genes and colors corresponding to scaled and centered expression levels.
#'
#' The expression data are taken from the logcounts assay of the input \code{sce} object. The colors to represent the expression levels can be controlled by the parameter \code{heatmap_colors}.
#'
#' The cells can be order by pseudotime, cell type or any other variable that makes sense in the given context. The name of metadata column in the input \code{sce} object containing the desired cell property can be specified by the parameter \code{order_by} (default: "pseudotime"). This cell property is visualized as a rug plot under the heatmap, the colors for the rug can be controlled by the parameter \code{annotation_colors}.
#'
#' Expression levels are clipped to the range of mean \eqn{\pm} 3\eqn{\sigma} per gene, this can be changed via the parameter \code{clip}. Clipping aids visualization by preventing the outlier data points from squishing the rest of the data into a small color range. If clipping is not desired, please set \code{clip} to Inf.
#'
#' @param genes Character vector, the names of the genes for which the expression profiles should be plotted.
#' @param sce \code{\link{SingleCellExperiment}} object containing the expression data (logcounts and metadata) for all network genes. Required metadata column:
#'\describe{
#' \item{\{\{order_by\}\}}{Numeric or character, the column by which the cells should be ordered on the plot. Typically contains the inferred pseudotime, or the cell type labels.}
#' }
#' @param order_by Character specifying the column in the metadata of \code{sce} by which the cells should be ordered on the plot. This column typically contains the inferred pseudotime (default: "pseudotime"), or the cell type labels.
#' @param heatmap_colors Character vector, the heatmap colors for the expression levels. The vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{\link{scale_color_gradientn}}.
#' @param annotation_colors Character vector, the colors for the variable specified by \code{order_by}. If the variable is discrete, the vector should contain as many colors as there are unique values of the variable, if the variable is continuous, the vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{\link{scale_color_gradientn}}
#' @param clip Numeric specifying the degree of clipping. For each gene, the expression level values that are more standard deviations away from the mean than \code{clip} are treated as NA. The default is 3 meaning that expression levels are clipped to the range of mean \eqn{\pm} 3\eqn{\sigma} per gene.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A heatmap as a \code{\link{ggplot}} object showing the expression profiles of the input genes across all cells.
#' @export
#'
#' @examples
#' plotExprHeatmap(c("POU5F1","DNMT3B","TERF1","TDGF1","L1TD1","VIM","MAP1B","MARCKS","PTN","CDH2"),
#'                 sce)
#' @family functions to plot gene expression profiles
plotExprHeatmap <- function(genes, sce, order_by = "pseudotime", heatmap_colors = NULL, annotation_colors = NULL, clip = 3, font_size = 14) {

  # check input data
  if (!inherits(genes, "character"))
    stop("The argument \"genes\" should be a character vector.")

  if (!inherits(sce, "SingleCellExperiment"))
    stop("The argument \"sce\" should be of class \"SingleCellExperiment\".")

  if (!("logcounts" %in% names(SummarizedExperiment::assays(sce))))
    stop("The argument \"sce\" should contain the assay \"logcounts\".")

  if (!inherits(order_by, "character") || length(order_by) != 1)
    stop("The argument \"order_by\" should be a string specifying the metadata column of \"sce\" by which the cells should be ordered on the plot.")

  if (is.null(sce[[order_by]]))
    stop("The argument \"sce\" should contain the metadata column specified by \"order_by\".")

  genes_not_found <- genes[!genes %in% rownames(sce)]

  if (length(genes_not_found) > 0)
    stop(paste0("The following input genes cannot be found in \"sce\": ", paste(genes_not_found, collapse = ", "), "."))

  if (!is.null(heatmap_colors) && (!inherits(heatmap_colors, "character") || any(!areColors(heatmap_colors))))
    stop("The argument \"heatmap_colors\" should be a character vector of valid color representations.")

  if (!is.null(annotation_colors) && (!inherits(annotation_colors, "character") || any(!areColors(annotation_colors))))
    stop("The argument \"annotation_colors\" should be a character vector of valid color representations.")

  if (!is.null(annotation_colors) && !inherits(sce[[order_by]], "numeric") && length(annotation_colors) != length(unique(sce[[order_by]])))
    stop("The argument \"annotation_colors\" should contain as many colors as there are unique values in the metadata column of \"sce\" specified by \"order_by\".")

  if (!inherits(clip, "numeric") || length(clip) != 1 || clip <= 0)
    stop("The argument \"clip\" should be a positive numeric value.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  # avoid NSE notes in R CMD check
  gene_name = NULL

  # get the number of genes
  n_genes <- length(genes)

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

  # clipping threshold
  max_abs_expr <- min(max(abs(expr$expr)), clip)

  # if no heatmap colors are provided, take the default
  if (is.null(heatmap_colors))
    heatmap_colors <- eigengene_colors

  # heatmap
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

  # if the cells are ordered by pseudotime (or any other continuous variable), add continuous color scale
  if (is.numeric(expr$order_by)) {

    # if no annotation colors are provided, take the default
    if (is.null(annotation_colors))
      annotation_colors <- pseudotime_colors

    return(
      p + ggplot2::scale_color_gradientn(colors = annotation_colors, name = order_by)
    )

  # if the cells are ordered by cell type (or any other discrete variable), add discrete color scale
  } else {

    # if no annotation colors are provided, take the default
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
#' The function takes a data frame containing any of these summarized module expression profiles as input (normally the output of \code{\link{calculateEigengenes}}) and plots it as a heatmap, with columns corresponding to cells, rows corresponding to modules and colors corresponding to scaled and centered expression levels.
#'
#' The column of the data frame containing the chosen type of summarized expression profile can be specified by the parameter \code{expr_column} (default: "eigengene"). The colors to represent the expression levels can be controlled by the parameter \code{heatmap_colors}.
#'
#' The cells can be order by pseudotime, cell type or any other variable that makes sense in the given context, the name of the column containing the desired cell property can be specified by the parameter \code{order_by} (default: "pseudotime").  This cell property is visualized as a rug plot under the heatmap, the colors for the rug can be be controlled by the parameter \code{annotation_colors}.
#'
#' If the eigengene (or mean expression) has been calculated for positively and negatively regulated targets separately, then these will appear as separate rows of the heatmap.
#'
#' Expression levels are clipped to the range of mean \eqn{\pm} 3\eqn{\sigma} per gene, this can be changed via the parameter \code{clip}. Clipping aids visualization by preventing the outlier data points from squishing the rest of the data into a small color range. If clipping is not desired, please set \code{clip} to Inf.
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
#' @param heatmap_colors Character vector, the heatmap colors for the variable specified by \code{expr_column}. The vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{\link{scale_color_gradientn}}.
#' @param annotation_colors Character vector, the colors for the variable specified by \code{order_by}. If the variable is discrete, the vector should contain as many colors as there are unique values of the variable, if the variable is continuous, the vector can contain any number of colors that will be passed on to and converted into a continuous scale by \code{\link{scale_color_gradientn}}.
#' @param clip Numeric specifying the degree of clipping. For each gene, the expression level values that are more standard deviations away from the mean than \code{clip} are treated as NA. The default is 3 meaning that expression levels are clipped to the range of mean \eqn{\pm} 3\eqn{\sigma} per gene.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A heatmap as a \code{\link{ggplot}} object showing the summarized expression profiles of the modules across all cells.
#' @export
#'
#' @examples plotEigengeneHeatmap(eigengenes)
#' @family functions to plot eigengene profiles
#' @references
#' Zhang, B., & Horvath, S. (2005). A general framework for weighted gene co-expression network analysis. Statistical Applications in Genetics and Molecular Biology, 4, 17-60. https://doi.org/10.2202/1544-6115.1128
plotEigengeneHeatmap <- function(eigengenes, expr_column = "eigengene", order_by = "pseudotime", heatmap_colors = NULL, annotation_colors = NULL, clip = 3, font_size = 14) {

  # check input data
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

  if (!inherits(clip, "numeric") || length(clip) != 1 || clip <= 0)
    stop("The argument \"clip\" should be a positive numeric value.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  # get the number of modules
  n_modules <- length(unique(eigengenes$module))

  # scale
  eigengenes <- eigengenes %>%
    dplyr::group_by(.data[["module"]]) %>%
    dplyr::mutate({{expr_column}} := scale(.data[[expr_column]])) %>%
    dplyr::ungroup()

  # clipping threshold
  max_abs_expr <- min(max(abs(eigengenes[[expr_column]])), clip)

  # if no heatmap colors are provided, take the default
  if (is.null(heatmap_colors))
    heatmap_colors <- eigengene_colors

  # heatmap
  p <- eigengenes %>%
    ggplot2::ggplot(ggplot2::aes(x = stats::reorder(.data[["cell"]], as.numeric(.data[[order_by]])), y = .data[["module"]], fill = .data[[expr_column]])) +
    ggplot2::geom_tile() +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::scale_fill_gradientn(colors = heatmap_colors,  limits = c(-max_abs_expr, max_abs_expr), guide = ggplot2::guide_colorbar(order = 1)) +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(2, 2, 2, 2)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.85, 0)), limits = rev) +
    ggplot2::geom_rug(ggplot2::aes(color = .data[[order_by]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(1 / (3 * n_modules), "npc")) +
    ggplot2::xlab("cells") +
    ggplot2::ylab("modules")

  # if the cells are ordered by pseudotime (or any other continuous variable), add continuous color scale
  if (is.numeric(eigengenes[[order_by]])) {

    # if no annotation colors are provided, take the default
    if (is.null(annotation_colors))
      annotation_colors <- pseudotime_colors

    return(
      p + ggplot2::scale_color_gradientn(colors = annotation_colors)
    )

  # if the cells are ordered by cell type (or any other discrete variable), add discrete color scale
  } else {

    # if no annotation colors are provided, take the default
    if (is.null(annotation_colors))
      annotation_colors <- cell_type_color_ramp(seq(0, 1, length = length(unique(eigengenes[[order_by]]))))

    p + ggplot2::scale_color_manual(values = annotation_colors) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 8.5))) +
      ggplot2::theme(legend.key.width = ggplot2::unit(0.035, "npc"))

  }

}


#' Plot expression profiles along a pseudotime trajectory
#'
#' Plots the expression profiles of one or more genes along a pseudotime trajectory per species.
#'
#' The function plots the expression profiles as smoothed curves colored by species and faceted by gene.
#'
#' The species, pseudotime and cell type information are taken from the metadata slot of the input \code{sce} object, and the expression data are taken from the logcounts assay of the input \code{sce} object.
#'
#' The smoothed expression profiles are fitted per species and gene using \code{\link{loess}} with the formula "expression ~ pseusotime". The 95\% confidence intervals of the fitted lines are calculated using \code{\link{predict}} and shown as lightly colored areas around the lines.
#'
#' If a cell type metadata column is specified by the parameter \code{cell_type_column}, the cell are shown as a rug plot along the pseudotime axis colored by cell types. The colors for the rug can be be controlled by the parameter \code{cell_type_colors}.
#'
#' @param genes Character vector, the names of the genes for which the expression profiles should be plotted.
#' @param sce \code{\link{SingleCellExperiment}} object containing the expression data (raw counts, logcounts and metadata) for all network genes. Required metadata columns:
#'\describe{
#' \item{species}{Character, the name of the species.}
#' \item{\{\{pseudotime_column\}\}}{Numeric, inferred pseudotime.}
#' \item{\{\{cell_type_column\}\}}{Character, cell type annotation (optional).}
#' }
#' @param pseudotime_column Character, the name of the pseudotime column in the metadata of \code{sce} (default: "pseudotime").
#' @param cell_type_column Character, the name of the cell type annotation column in the metadata of \code{sce} (default: "cell_type", if there is no cell type annotation available or the user wants to omit the cell type rug plot, this parameter should be set to NULL).
#' @param species_colors Character vector, colors per species.
#' @param cell_type_colors Character vector, colors per cell type.
#' @param font_size Numeric, font size (default: 14).
#' @param ncol Integer, the number of columns the genes (facets) should be organized into (default: 1).
#'
#' @return A curve plot as a \code{\link{ggplot}} object showing the expression profiles of the genes along the pseudotime trajectory per species.
#' @export
#'
#' @examples plotExprAlongPseudotime(regulators, sce)
#' @family functions to plot gene expression profiles
plotExprAlongPseudotime <- function(genes, sce, pseudotime_column = "pseudotime", cell_type_column = "cell_type", species_colors = NULL, cell_type_colors = NULL, font_size = 14, ncol = 1) {

  # check input data
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

  if (length(ncol) != 1 || (!inherits(ncol, "integer") && !(inherits(ncol, "numeric") & ncol == round(ncol))) || ncol < 1)
    stop("The argument \"ncol\" should be a positive integer.")

  # avoid NSE notes in R CMD check
  gene_name = species_name = NULL

  # if no species colors are provided, take the default
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

  # initialize plot
  p <- expr %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["pseudotime"]], y = .data[["expr"]])) +
    ggplot2::xlab(pseudotime_column) +
    ggplot2::ylab("expression (logcounts)") +
    ggplot2::theme_bw(base_size = font_size)

  # if cell type information is present, add rug to the plot
  if (!is.null(cell_type_column)) {

    # if no cell type colors are provided, take the default
    if (is.null(cell_type_colors))
      cell_type_colors <- cell_type_color_ramp(seq(0, 1, len = length(unique(sce[[cell_type_column]]))))

    # calculate rug positions
    rug_positions <- loess_fit %>%
      dplyr::group_by(.data[["gene"]]) %>%
      dplyr::summarize(pos = min(.data[["lwr"]]) - 0.2*(max(.data[["upr"]]) - min(.data[["lwr"]]))) %>%
      tibble::deframe()

    # rug plot
    p <- p +
      ggplot2::geom_rug(data = expr %>% dplyr::mutate(expr = rug_positions[.data[["gene"]]]),
                        ggplot2::aes(color = .data[["cell_type"]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(0.1, "npc"), linewidth = 0.15) +
      ggplot2::scale_color_manual(values = cell_type_colors, name = cell_type_column) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
      ggnewscale::new_scale_color()

  }

  # smooth curves
  p +
    ggplot2::geom_line(data = loess_fit, ggplot2::aes(color = .data[["species"]], group = .data[["species"]]), linewidth = 1) +
    ggplot2::geom_ribbon(data = loess_fit, ggplot2::aes(ymin = .data[["lwr"]], ymax = .data[["upr"]], fill = .data[["species"]]), alpha = 0.1) +
    ggplot2::scale_color_manual(values = species_colors) +
    ggplot2::scale_fill_manual(values = species_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
    ggplot2::facet_wrap(~.data[["gene"]], scales = "free_y", ncol = ncol, strip.position = "left") +
    ggplot2::scale_y_continuous(position = "right")

}


#' Plot module eigengenes along a pseudotime trajectory
#'
#' Plots the eigengenes (or other types of summarized module expression profiles) along a pseudotime trajectory per species.
#'
#' A concept adapted from WGCNA, the eigengene summarizes the expression profile of an entire module, and it is calculated as the first principal component of the module expression data (see also \code{\link{calculateEigengenes}}). Other possible ways of representing the expression profile of a module include the mean expression and the regulator expression.
#'
#' The function takes a data frame containing any of these summarized module expression profiles as input (normally the output of \code{\link{calculateEigengenes}}). The column containing the chosen type of summarized expression profile can be specified by the parameter \code{expr_column} (default: "eigengene").
#'
#' The values of this chosen metric are plotted as smoothed curves colored by species and faceted by module. The smoothed expression profiles are fitted per species and module using \code{\link{loess}} with the formula "expression ~ pseusotime". The 95\% confidence intervals of the fitted lines are calculated using \code{\link{predict}} and shown as lightly colored areas around the lines.
#'
#' If a cell type column is specified by the parameter \code{cell_type_column}, the cell are shown as a rug plot along the pseudotime axis colored by cell types. The colors for the rug can be be controlled by the parameter \code{cell_type_colors}.
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
#' @param cell_type_column Character, the name of the cell type annotation column in \code{eigengenes} (default: "cell_type", if there is no cell type annotation available or the user wants to omit the cell type rug plot, this parameter should be set to NULL).
#' @param species_colors Character vector, colors per species.
#' @param cell_type_colors Character vector, colors per cell type.
#' @param font_size Numeric, font size (default: 14).
#' @param ncol Integer, the number of columns the modules (facets) should be organized into (default: 1).
#'
#' @return A curve plot as a \code{\link{ggplot}} object showing the summarized expression profiles of the modules along the pseudotime trajectory per species.
#' @export
#'
#' @examples
#' plotEigengenesAlongPseudotime(eigengenes_per_species)
#' plotEigengenesAlongPseudotime(eigengenes_per_species, expr_column = "mean_expr")
#' plotEigengenesAlongPseudotime(eigengenes_per_species, cell_type_column = NULL)
#' @family functions to plot eigengene profiles
#' @references
#' Zhang, B., & Horvath, S. (2005). A general framework for weighted gene co-expression network analysis. Statistical Applications in Genetics and Molecular Biology, 4, 17-60. https://doi.org/10.2202/1544-6115.1128
plotEigengenesAlongPseudotime <- function(eigengenes, expr_column = "eigengene",  pseudotime_column = "pseudotime", cell_type_column = "cell_type", species_colors = NULL, cell_type_colors = NULL, font_size = 14, ncol = 1) {

  # check input data
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

  if (length(ncol) != 1 || (!inherits(ncol, "integer") && !(inherits(ncol, "numeric") & ncol == round(ncol))) || ncol < 1)
    stop("The argument \"ncol\" should be a positive integer.")

  # avoid NSE notes in R CMD check
  module_name = species_name = NULL

  # get all modules and species
  eigengenes <- eigengenes %>%
    dplyr::rename(expr = expr_column,
                  pseudotime = pseudotime_column)
  module_names <- unique(eigengenes$module)
  species_names <- unique(eigengenes$species)

  # if no species colors are provided, take the default
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

  # initialize plot
  p <- eigengenes %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["pseudotime"]], y = .data[["expr"]])) +
    ggplot2::xlab(pseudotime_column) +
    ggplot2::ylab(expr_column) +
    ggplot2::theme_bw(base_size = font_size)

  # if cell type information is present, add rug to the plot
  if (!is.null(cell_type_column)) {

    # if no cell type colors are provided, take the default
    if (is.null(cell_type_colors))
      cell_type_colors <- cell_type_color_ramp(seq(0, 1, len = length(unique(eigengenes[[cell_type_column]]))))

    # calculate rug positions
    rug_positions <- loess_fit %>%
      dplyr::group_by(.data[["module"]]) %>%
      dplyr::summarize(pos = min(.data[["lwr"]]) - 0.2*(max(.data[["upr"]]) - min(.data[["lwr"]]))) %>%
      tibble::deframe()

    # rug plot
    p <- p +
      ggplot2::geom_rug(data = eigengenes %>% dplyr::mutate(expr = rug_positions[.data[["module"]]]),
                        ggplot2::aes(color = .data[[cell_type_column]]),  sides = "b", show.legend = TRUE, length = ggplot2::unit(0.1, "npc"), linewidth = 0.15) +
      ggplot2::scale_color_manual(values = cell_type_colors,  name = cell_type_column) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
      ggnewscale::new_scale_color()

  }

  # smooth curves
  p +
    ggplot2::geom_line(data = loess_fit, ggplot2::aes(color = .data[["species"]], group = .data[["species"]]), linewidth = 1) +
    ggplot2::geom_ribbon(data = loess_fit, ggplot2::aes(ymin = .data[["lwr"]], ymax = .data[["upr"]], fill = .data[["species"]]), alpha = 0.1) +
    ggplot2::scale_color_manual(values = species_colors) +
    ggplot2::scale_fill_manual(values = species_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
    ggplot2::facet_wrap(~.data[["module"]], scales = "free_y", ncol = ncol, strip.position = "left") +
    ggplot2::scale_y_continuous(position = "right")

}


#' Plot expression levels per cell type
#'
#' Plots the expression levels of one or more genes per cell type and species, and thus allows the expression patterns to be visually compared across species.
#'
#' The function produces a violin plot of expression levels per cell type and species, faceted by gene. The colors for the species can be controlled by the parameter \code{species_colors}.
#'
#' The species and cell type information are taken from the metadata slot of the input \code{sce} object, and the expression data are taken from the logcounts assay of the input \code{sce} object.
#'
#' @param genes Character vector, the names of the genes for which the expression profiles should be plotted.
#' @param sce \code{\link{SingleCellExperiment}} object containing the expression data (raw counts, logcounts and metadata) for all network genes. Required metadata columns:
#'\describe{
#' \item{species}{Character, the name of the species.}
#' \item{\{\{cell_type_column\}\}}{Character, cell type annotation.}
#' }
#' @param cell_type_column Character, the name of the cell type annotation column in the metadata of \code{sce}.
#' @param species_colors Character vector, colors per species.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A violin plot as a \code{\link{ggplot}} object showing the expression levels of the genes per cell type and species.
#' @export
#'
#' @examples plotExprPerCellType(regulators, sce)
#' @family functions to plot gene expression profiles
plotExprPerCellType <- function(genes, sce, cell_type_column = "cell_type", species_colors = NULL, font_size = 14) {

  # check input data
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

  # avoid NSE notes in R CMD check
  gene = NULL

  # if no species colors are provided, take the default
  if (is.null(species_colors))
    species_colors <- species_color_ramp(seq(0, 1, length = dplyr::n_distinct(sce$species)))

  # extract metadata and logcounts of the genes of interest from the SCE object
  expr <- foreach::foreach(gene = genes,
                           .combine = dplyr::bind_rows,
                           .multicombine = TRUE) %do% {

                             data.frame(species = sce$species,
                                        gene = gene,
                                        expr = SingleCellExperiment::logcounts(sce)[gene,],
                                        cell_type = sce[[cell_type_column]])

                           } %>%
    dplyr::mutate(gene = factor(.data[["gene"]], levels = genes))

  # violin plot
  expr %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["cell_type"]], y = .data[["expr"]], fill = .data[["species"]])) +
    ggplot2::geom_violin(draw_quantiles = 0.5) +
    ggplot2::scale_fill_manual(values = species_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
    ggplot2::xlab(cell_type_column) +
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
#' The function takes a data frame containing any of these summarized module expression profiles as input (normally the output of \code{\link{calculateEigengenes}}). The column containing the chosen type of summarized expression profile can be specified by the parameter \code{expr_column} (default: "eigengene").
#'
#' The values of the chosen metric are plotted as violin plots per cell type and species, faceted by module. The colors for the species can be controlled by the parameter \code{species_colors}.
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
#' @return A violin plot as a \code{\link{ggplot}} object showing the summarized expression levels of the modules per cell type and species.
#' @export
#'
#' @examples plotEigengenesPerCellType(eigengenes)
#' @family functions to plot eigengene profiles
#' @references
#' Zhang, B., & Horvath, S. (2005). A general framework for weighted gene co-expression network analysis. Statistical Applications in Genetics and Molecular Biology, 4, 17-60. https://doi.org/10.2202/1544-6115.1128
plotEigengenesPerCellType <- function(eigengenes, expr_column = "eigengene", cell_type_column = "cell_type", species_colors = NULL, font_size = 14) {

  # check input data
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

  # if no species colors are provided, take the default
  if (is.null(species_colors))
    species_colors <- species_color_ramp(seq(0, 1, len = length(unique(eigengenes$species))))

  # violin plot
  eigengenes %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["cell_type"]], y = .data[[expr_column]], fill = .data[["species"]])) +
    ggplot2::geom_violin(draw_quantiles = 0.5) +
    ggplot2::scale_fill_manual(values = species_colors) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 0.7))) +
    ggplot2::xlab(cell_type_column) +
    ggplot2::ylab(expr_column) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::facet_wrap(~.data[["module"]], scales = "free_y", ncol = 1, strip.position = "left") +
    ggplot2::scale_y_continuous(position = "right")

}
