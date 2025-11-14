# Consensus network

Consensus network of the 9 primate replicates in the example dataset.
For each edge, the consensus adjacency was calculated as the weighted
average of replicate-wise adjacencies using weights that correct for 1)
the phylogenetic distances between species and 2) the different numbers
of replicates per species. If an edge was not detected in certain
replicate, the adjacency of that replicate was regarded as 0 for the
calculation of the consensus. The directionality of each edge was
determined based on a modified Spearman's correlation between the
corresponding 2 genes' expression profiles (positive expression
correlation - activating interaction, negative expression correlation -
repressing interaction). The correlations were calculated per replicate,
then the mean correlation was taken across all replicates.

## Usage

``` r
consensus_network
```

## Format

An \[igraph\] object with 300 nodes, and has 1 node attribute and 3 edge
attributes:

Node attributes:

- name:

  Name of the node (gene).

Edge attributes:

- weight:

  Consensus edge weight/adjacency, the weighted average of
  replicate-wise adjacencies.

- rho:

  Approximate Spearman's correlation coefficient of the 2 genes'
  expression profiles that form the edge.

- p.adj:

  BH-corrected approximate p-value of rho.

- direction:

  Direction of the interaction between the 2 genes that form the edge
  ("+" or "-").
