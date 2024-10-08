% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotExpr.R
\name{plotEigengenesPerCellType}
\alias{plotEigengenesPerCellType}
\title{Plot eigengenes per cell type}
\usage{
plotEigengenesPerCellType(
  eigengenes,
  expr_column = "eigengene",
  cell_type_column = "cell_type",
  species_colors = NULL,
  font_size = 14
)
}
\arguments{
\item{eigengenes}{Data frame of eigengenes, required columns:
\describe{
\item{cell}{Character, the cell barcode.}
\item{species}{Character, the name of the species.}
\item{\{\{cell_type_column\}\}}{Character, cell type annotation.}
\item{module}{Character, transcriptional regulator and in case the eigengene was calculated for the positively or negatively regulated targets only, the direction of regulation (format: nameOfRegulator(+) or nameOfRegulator(-)).}
\item{\{\{expr_column\}\}}{Numeric, summarized module expression profiles (typically the eigengene, the mean expression of the module, or the expression of the regulator).}
}}

\item{expr_column}{Character specifying which column of \code{eigengenes} the heatmap should be colored by. This column is expected to contain summarized module expression profiles, typically the eigengene (default: "eigengene"), the mean expression of the module, or the expression of the regulator.}

\item{cell_type_column}{Character, the name of the cell type annotation column in \code{eigengenes}.}

\item{species_colors}{Character vector, colors per species.}

\item{font_size}{Numeric, font size (default: 14).}
}
\value{
A violin plot as a \code{\link{ggplot}} object showing the summarized expression levels of the modules per cell type and species.
}
\description{
Plots the eigengenes (or other types of summarized module expression profiles) per cell type and species, and thus allows the expression patterns to be visually compared across species.
}
\details{
A concept adapted from WGCNA, the eigengene summarizes the expression profile of an entire module, and it is calculated as the first principal component of the module expression data (see also \code{\link{calculateEigengenes}}). Other possible ways of representing the expression profile of a module include the mean expression and the regulator expression.

The function takes a data frame containing any of these summarized module expression profiles as input (normally the output of \code{\link{calculateEigengenes}}). The column containing the chosen type of summarized expression profile can be specified by the parameter \code{expr_column} (default: "eigengene").

The values of the chosen metric are plotted as violin plots per cell type and species, faceted by module. The colors for the species can be controlled by the parameter \code{species_colors}.
}
\examples{
plotEigengenesPerCellType(eigengenes)
}
\references{
Zhang, B., & Horvath, S. (2005). A general framework for weighted gene co-expression network analysis. Statistical Applications in Genetics and Molecular Biology, 4, 17-60. https://doi.org/10.2202/1544-6115.1128
}
\seealso{
Other functions to plot eigengene profiles: 
\code{\link{plotEigengeneHeatmap}()},
\code{\link{plotEigengenesAlongPseudotime}()}
}
\concept{functions to plot eigengene profiles}
