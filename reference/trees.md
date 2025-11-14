# Trees of the pruned modules

Neighbor-joining trees representing the similarities of connectivity
patterns across the 9 primate replicates for the original (i.e. not
jackknifed) pruned modules. For each of the modules, the trees were
inferred based on the preservation statistic cor.kIM (correlation of
intramodular connectivities): first, the preservation scores were
calculated between all possible replicate pairs, then they were
converted into a distance matrix of replicates, and finally trees were
reconstructed based on this distance matrix using the neighbor-joining
algorithm. The result is a single tree per module where the tips
represent the replicates and the branch lengths represent the
dissimilarity of connectivity patterns between these replicates.

## Usage

``` r
trees
```

## Format

A named list with 12 elements containing the neighbor-joining trees as
\[phylo\] objects.
