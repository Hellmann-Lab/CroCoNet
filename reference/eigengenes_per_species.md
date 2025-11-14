# Eigengenes per species

Eigengenes of the pruned modules calculated per species using the
activated targets in each module. An eigengene summarizes the expression
profile of the module as a whole, mathematically it is the first
principal component of the module expression data (i.e. the scaled and
centered logcounts subsetted for the activated targets and the regulator
of the given module). In this case, the principal component was
calculated for each species separately.

## Usage

``` r
eigengenes_per_species
```

## Format

A data frame with 6300 rows and 8 columns:

- cell:

  Character, the cell barcode.

- species:

  Character, the name of the species.

- pseudotime:

  Numeric, inferred pseudotime.

- cell_type:

  Character, cell type annotation.

- module:

  Character, transcriptional regulator and direction of regulation (in
  this case always nameOfRegulator(+)).

- eigengene:

  Numeric, the eigengene (i.e. the first principal component of the
  scaled and centered logcounts) of the module.

- mean_expr:

  Numeric, the mean of the scaled and centered logcounts across all
  genes in the module.

- regulator_expr:

  Numeric, the scaled and centered logcounts of the regulator.
