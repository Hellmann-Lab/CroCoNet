# Calculate module eigengenes

Calculates the module eigengene for each input module.

## Usage

``` r
calculateEigengenes(
  pruned_modules,
  sce,
  direction_of_regulation = "+_only",
  module_names = as.character(unique(pruned_modules$regulator)),
  per_species = FALSE,
  pseudotime_column = "pseudotime",
  cell_type_column = "cell_type",
  n_cores = 1L
)
```

## Arguments

- pruned_modules:

  Data frame of the pruned modules, required columns:

  regulator

  :   Character, transcriptional regulator.

  target

  :   Character, target gene of the transcriptional regulator (member of
      the regulator's pruned module).

  direction

  :   Character specifying the direction of interaction between the
      regulator and the target, either "+" or "-" (only required if
      `direction_of_regulation` is set to "+\_only" or
      "+-\_separately").

- sce:

  `SingleCellExperiment` object containing the expression data
  (logcounts and metadata) for all network genes. Required metadata
  columns:

  species

  :   Character, the name of the species.

  {{pseudotime_column}}

  :   Numeric, inferred pseudotime (optional).

  {{cell_type_column}}

  :   Character, cell type annotation (optional).

- direction_of_regulation:

  Character specifying how positively and negatively regulated targets
  of the same transcriptional regulator should be treated, one of
  "+\_only", "+-\_separately", "all_together". If "+\_only", the
  eigengene is calculated only for the positively regulated targets, the
  negatively regulated targets are removed. If "+-\_separately", the
  eigengene is calculated separately for the positively and negatively
  regulated targets. If "all_together", the eigengene is calculated for
  all targets, irrespective of the direction of regulation.

- module_names:

  Character vector, the names of the modules for which the eigengenes
  should be calculated (default: all unique names in the column
  `regulator` of `pruned_modules`).

- per_species:

  Logical, if FALSE (default), the eigengenes are calculated across all
  cells, if TRUE, the eigengenes are calculated per species.

- pseudotime_column:

  Character, the name of the pseudotime column in the metadata of `sce`
  (default: "pseudotime", if there is no pseudotime column, it should be
  set to NULL).

- cell_type_column:

  Character, the name of the cell type annotation column in the metadata
  of `sce` (default: "cell_type", if there is no cell type column, it
  should be set to NULL).

- n_cores:

  Integer, the number of cores (default: 1).

## Value

A data frame of eigengenes with the following columns:

- cell:

  Character, the cell barcode.

- species:

  Character, the name of the species.

- {{pseudotime_column}}:

  Numeric, inferred pseudotime (only present if `pseudotime_column` is
  not NULL).

- {{cell_type_column}}:

  Character, cell type annotation (only present if `cell_type_column` is
  not NULL).

- module:

  Character, transcriptional regulator and in case the eigengene was
  calculated for the positively or negatively regulated targets only,
  the direction of interaction (format: nameOfRegulator(+) or
  nameOfRegulator(-)).

- eigengene:

  Numeric, the eigengene (i.e. the first principal component of the
  scaled and centered logcounts) of the module. In case `per_species` is
  TRUE, it is calculated for each species separately.

- mean_expr:

  Numeric, the mean of the scaled and centered logcounts across all
  genes in the module.

- regulator_expr:

  Numeric, the scaled and centered logcounts of the regulator.

## Details

A concept adapted from WGCNA, the eigengene summarizes the expression
profile of an entire module, and it is calculated as the first principal
component of the module expression data. Effectively, it is a weighted
mean of the individual genes' expression profiles.

As the first step, the logcounts are subsetted to keep only module
member genes and the resulting count matrix is scaled and centered per
gene. Next, singular value decomposition is performed on the scaled and
centered count matrix using [`svd`](https://rdrr.io/r/base/svd.html),
with `nu` = 1 and `nv` = 1 (only 1 left and 1 right singular vector
computed). The right singular vector is taken as the eigengene. Finally,
this eigengene is aligned along the average expression of the module: if
the correlation of the two vectors is negative, the eigengene is
negated, if the correlation is positive, the eigengene is kept as it is.

If a module contains both activated and repressed targets of the
transcriptional regulator, calculating the eigengene (or any other
summary expression profiles) across both directions of regulation does
not make biological sense and leads to the dilution of signal. It makes
more sense to calculate the eigengene either for the activated targets
only (`direction_of_regulation` = "+\_only") or for the activated and
repressed targets separately (`direction_of_regulation` =
"+-\_separately"). In these cases, the `module` column of the output
will specify not only the name of the transcriptional regulator but also
the direction of regulation (format: nameOfRegulator(+) or
nameOfRegulator(-)). If an eigengene across all targets is desired
irrespective of the direction of regulation, `direction_of_regulation`
should be set to "all_together".

If the aim is to compare the eigengenes across species, it is
recommended to calculate the eigengenes per species by setting
`per_species` to TRUE. In this case, the scaling, centering and SVD is
performed for each species separately.

If the user plans to plot the eigengene along pseudotime and/or cell
types, the corresponding columns of the `sce` object can be specified by
the arguments of `pseudotime_column` and `cell_type_column`, and then
the pseudotime and cell type information of each cell will be added to
the output.

## References

Zhang, B., & Horvath, S. (2005). A general framework for weighted gene
co-expression network analysis. Statistical Applications in Genetics and
Molecular Biology, 4, 17-60. https://doi.org/10.2202/1544-6115.1128

## Examples

``` r
eigengenes <- calculateEigengenes(pruned_modules, sce)
eigengenes_per_dir <- calculateEigengenes(pruned_modules, sce, "+-_separately")
eigengenes_per_species <- calculateEigengenes(pruned_modules, sce, per_species = TRUE)
```
