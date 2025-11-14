# List of raw networks

List of networks per replicate with raw edge weights. The networks were
inferred using GRNBoost2 based on a subset of the primate neural
differentiation scRNA-seq dataset. All 300 genes in the subsetted data
were used as potential regulators. To circumvent the stochastic nature
of the algorithm, GRNBoost2 was run 10 times on the same count matrices,
then the results were averaged across runs, and rarely occurring edges
were removed altogether. In addition, edges inferred between the same
gene pair but in opposite directions were also averaged.

## Usage

``` r
network_list_raw
```

## Format

A named list of 9 \[igraph\] objects. Each network contains 300 nodes,
and has 1 node and 2 edge attributes:

Node attributes:

- name:

  Name of the node (gene).

Edge attributes:

- weight:

  Edge weight, the importance score calculate by GRNBoost2.
