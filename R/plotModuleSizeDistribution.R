#' Plot the distribution of module sizes
#'
#' Plots the distribution of module sizes (i.e. the number of target genes assigned to each regulator) as a histogram.
#'
#' @param modules Data frame of the pruned modules, required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Module size, the numer of target genes assigned to a regulator.}
#' }
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A histogram as a \code{\link{ggplot}} object showing the distribution of module sizes.
#' @export
#'
#' @examples plotModuleSizeDistribution(pruned_modules)
plotModuleSizeDistribution <- function(modules, font_size = 14) {

  # check input data
  if (!is.data.frame(modules))
    stop("The argument \"modules\" should be a data frame.")

  if (any(!c("regulator", "module_size") %in% colnames(modules)))
    stop("The argument \"modules\" should contain the columns \"regulator\" and \"module_size\".")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  # create histogram
  modules %>%
    dplyr::distinct(.data[["regulator"]], .data[["module_size"]]) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["module_size"]])) +
    ggplot2::geom_histogram(color = "black", linewidth = 0.2, fill = "grey70", bins = max(modules$module_size) - min(modules$module_size) + 1) +
    ggplot2::theme_bw(base_size = font_size) +
    ggplot2::xlab("module size [number of genes]")

}
