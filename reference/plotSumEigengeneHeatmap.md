# Plot summarized eigengenes as a heatmap

Plots the mean eigengene (or the mean of other types of summarized
module expression profile) per cell type and species as a heatmap.

## Usage

``` r
plotSumEigengeneHeatmap(
  eigengenes,
  expr_column = "eigengene",
  cell_type_column = "cell_type",
  heatmap_colors = NULL,
  species_colors = NULL,
  cell_type_colors = NULL,
  clip = 3,
  font_size = 14,
  label_size = font_size/3
)
```

## Arguments

- eigengenes:

  Data frame of eigengenes, required columns:

  cell

  :   Character, the cell barcode.

  species

  :   Character, the name of the species.

  {{cell_type_column}}

  :   Character, cell type annotation.

  module

  :   Character, transcriptional regulator and in case the eigengene was
      calculated for the positively or negatively regulated targets
      only, the direction of regulation (format: nameOfRegulator(+) or
      nameOfRegulator(-)).

  {{expr_column}}

  :   Numeric, summarized module expression profiles (typically the
      eigengene, the mean expression of the module, or the expression of
      the regulator).

- expr_column:

  Character specifying the column of `eigengenes` by which the heatmap
  should be colored. This column is expected to contain summarized
  module expression profiles, typically the eigengene (default:
  "eigengene"), the mean expression of the module, or the expression of
  the regulator.

- cell_type_column:

  Character, the name of the cell type annotation column in `eigengenes`
  (default: "cell_type").

- heatmap_colors:

  Character vector, the heatmap colors for the summarized expression
  levels. The vector can contain any number of colors that will be
  passed on to and converted into a continuous scale by
  `scale_color_gradientn`.

- species_colors:

  Character vector, colors per species.

- cell_type_colors:

  Character vector, colors per cell type.

- clip:

  Numeric specifying the degree of clipping. For each gene, the
  expression level values that are more standard deviations away from
  the mean than `clip` are treated as NA. The default is 3 meaning that
  expression levels are clipped to the range of mean \\\pm\\ 3\\\sigma\\
  per gene.

- font_size:

  Numeric, font size (default: 14).

- label_size:

  Numeric, font size of the cell type labels (default: `font_size` / 5).

## Value

A heatmap as a `ggplot` object showing the summarized expression levels
of the modules per cell type and species.

## Details

A concept adapted from WGCNA, the eigengene summarizes the expression
profile of an entire module, and it is calculated as the first principal
component of the module expression data (see also
[`calculateEigengenes`](https://hellmann-lab.github.io/CroCoNet/reference/calculateEigengenes.md)).
Other possible ways of representing the expression profile of a module
include the mean expression and the regulator expression.

The function takes a data frame containing any of these summarized
module expression profiles as input (normally the output of
[`calculateEigengenes`](https://hellmann-lab.github.io/CroCoNet/reference/calculateEigengenes.md))
and plots it as a heatmap, with columns corresponding to cell types,
rows corresponding to modules and species, and colors corresponding to
the mean module expression across all cells of the given cell type and
species.

The column of the data frame containing the chosen type of summarized
expression profile can be specified by the parameter `expr_column`
(default: "eigengene"). The colors to represent the expression levels
can be controlled by the parameter `heatmap_colors`.

If the eigengene (or mean expression) has been calculated for positively
and negatively regulated targets separately, then these will appear as
separate rows of the heatmap.

The cell types and species are visualized as annotation graphics above
the heatmap using `cell_type_colors` and on the right hand side of the
heatmap using `species_colors`, respectively.

Expression levels are clipped to the range of mean \\\pm\\ 3\\\sigma\\
per module, this can be changed via the parameter `clip`. Clipping aids
visualization by preventing the outlier data points from squishing the
rest of the data into a small color range. If clipping is not desired,
please set `clip` to Inf.

## References

Zhang, B., & Horvath, S. (2005). A general framework for weighted gene
co-expression network analysis. Statistical Applications in Genetics and
Molecular Biology, 4, 17-60. https://doi.org/10.2202/1544-6115.1128

## See also

Other functions to plot eigengene profiles:
[`plotEigengeneHeatmap()`](https://hellmann-lab.github.io/CroCoNet/reference/plotEigengeneHeatmap.md),
[`plotEigengenesAlongPseudotime()`](https://hellmann-lab.github.io/CroCoNet/reference/plotEigengenesAlongPseudotime.md),
[`plotEigengenesViolin()`](https://hellmann-lab.github.io/CroCoNet/reference/plotEigengenesViolin.md),
[`plotSumEigengenesLine()`](https://hellmann-lab.github.io/CroCoNet/reference/plotSumEigengenesLine.md)

## Examples

``` r
plotSumEigengeneHeatmap(eigengenes)
```
