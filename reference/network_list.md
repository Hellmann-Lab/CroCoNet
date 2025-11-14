# List of networks

List of networks per replicate with edge weights re-scaled between 0 and
1, and gene pairs (edges) with overlapping annotations removed. The
networks were inferred using GRNBoost2 based on a subset of the primate
neural differentiation scRNA-seq dataset. All 300 genes in the subsetted
data were used as potential regulators. To circumvent the stochastic
nature of the algorithm, GRNBoost2 was run 10 times on the same count
matrices, then the results were averaged across runs, and rarely
occurring edges were removed altogether. In addition, edges inferred
between the same gene pair but in opposite directions were also
averaged. Edge weights were scaled by the maximum edge weight across all
replicates. Gene pairs that have overlapping annotations in any of the
species' genomes were removed from all networks.

## Usage

``` r
network_list
```

## Format

A named list of 9 \[igraph\] objects. Each network contains 300 nodes,
and has 1 node attribute and 3 edge attributes:

Node attributes:

- name:

  Name of the node (gene).

Edge attributes:

- weight:

  Edge weight, the importance score calculate by GRNBoost2 rescaled
  between 0 and 1.

- genomic_dist:

  Numeric, the genomic distance of the 2 genes that form the edge (Inf
  if the 2 genes are annotated on different chromosomes/contigs).
