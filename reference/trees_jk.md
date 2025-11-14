# Trees of the original and jackknifed pruned modules

Neighbor-joining trees representing the similarities of connectivity
patterns across the 9 primate replicates for the original and all
jackknifed versions of the pruned modules. The jackknifed versions of
the modules were created by removing each target gene assigned to a
module (the regulators were never excluded). For each of these module
versions, the trees were inferred based on the preservation statistic
cor.kIM (correlation of intramodular connectivities): first, the
preservation scores were calculated between all possible replicate
pairs, then they were converted into a distance matrix of replicates,
and finally trees were reconstructed based on this distance matrix using
the neighbor-joining algorithm. The result is a single tree for each
original or jackknife module where the tips represent the replicates and
the branch lengths represent the dissimilarity of connectivity patterns
between these replicates.

## Usage

``` r
trees_jk
```

## Format

A named list with 232 elements containing the neighbor-joining trees as
\[phylo\] objects.
