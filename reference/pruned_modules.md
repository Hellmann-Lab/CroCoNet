# Pruned modules

Pruned modules created by the dynamic pruning of the initial modules via
the "UIK_adj_kIM" method. The initial module members were filtered in
successive steps based on their adjacency to the regulator and their
intramodular connectivitiy alternately. In each step, the cumulative sum
curve based on one of these two characteristics was calculated per
module, then the targets below the knee point of the curve were kept.
This process was continued until the module sizes became as small as
possible without the median falling below a pre-defined minimum of 20
genes.

## Usage

``` r
pruned_modules
```

## Format

A data frame with 225 rows and 9 columns:

- regulator:

  Character, transcriptional regulator.

- target:

  Character, target gene of the transcriptional regulator (member of the
  regulator's pruned module).

- weight:

  Numeric, consensus edge weight/adjacency, the weighted average of
  replicate-wise adjacencies.

- rho:

  Approximate Spearman's correlation coefficient of the 2 genes'
  expression profiles that form the edge.

- p.adj:

  BH-corrected approximate p-value of rho.

- direction:

  Character specifying the direction of regulation between the regulator
  and the target, either "+" or "-".

- module_size:

  Module size, the numer of target genes assigned to a regulator.
