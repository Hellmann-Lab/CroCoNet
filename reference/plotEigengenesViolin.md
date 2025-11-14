# Plot eigengene distributions per cell type and species

Plots the distribution of eigengenes (or other types of summarized
module expression profiles) per cell type and species, and thus allows
the expression patterns to be visually compared across species.

## Usage

``` r
plotEigengenesViolin(
  eigengenes,
  expr_column = "eigengene",
  cell_type_column = "cell_type",
  species_colors = NULL,
  font_size = 14
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

  Character specifying which column of `eigengenes` the heatmap should
  be colored by. This column is expected to contain summarized module
  expression profiles, typically the eigengene (default: "eigengene"),
  the mean expression of the module, or the expression of the regulator.

- cell_type_column:

  Character, the name of the cell type annotation column in
  `eigengenes`.

- species_colors:

  Character vector, colors per species.

- font_size:

  Numeric, font size (default: 14).

## Value

A violin plot as a `ggplot` object showing the expression distributions
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
[`calculateEigengenes`](https://hellmann-lab.github.io/CroCoNet/reference/calculateEigengenes.md)).
The column containing the chosen type of summarized expression profile
can be specified by the parameter `expr_column` (default: "eigengene").

The values of the chosen metric are plotted as violin plots per cell
type and species, faceted by module. The colors for the species can be
controlled by the parameter `species_colors`.

## References

Zhang, B., & Horvath, S. (2005). A general framework for weighted gene
co-expression network analysis. Statistical Applications in Genetics and
Molecular Biology, 4, 17-60. https://doi.org/10.2202/1544-6115.1128

## See also

Other functions to plot eigengene profiles:
[`plotEigengeneHeatmap()`](https://hellmann-lab.github.io/CroCoNet/reference/plotEigengeneHeatmap.md),
[`plotEigengenesAlongPseudotime()`](https://hellmann-lab.github.io/CroCoNet/reference/plotEigengenesAlongPseudotime.md),
[`plotSumEigengeneHeatmap()`](https://hellmann-lab.github.io/CroCoNet/reference/plotSumEigengeneHeatmap.md),
[`plotSumEigengenesLine()`](https://hellmann-lab.github.io/CroCoNet/reference/plotSumEigengenesLine.md)

## Examples

``` r
plotEigengenesViolin(eigengenes)
```
