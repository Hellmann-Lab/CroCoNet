# Create consensus network

Integrates networks across different replicates and different species
into a single consensus network in a phylogeny-aware manner.

## Usage

``` r
createConsensus(network_list, replicate2species, tree = NULL)
```

## Arguments

- network_list:

  A named list of `igraph` objects containing the networks of all
  replicates.

- replicate2species:

  A data frame specifying which species each replicate belongs to,
  required columns:

  replicate

  :   Character, name of the replicate.

  species

  :   Character, name of the species.

- tree:

  Object of class `phylo`, the phylogenetic tree of the species.

## Value

Consensus network in an `igraph` format with the following edge
attributes:

- weight:

  Numeric, consensus edge weight/adjacency, the weighted average of
  replicate-wise edge weights.

- n_supporting_replicates:

  Integer, the number of replicates where the edge was detected (only
  added if more than 10% of the edges are not present in all
  replicates).

- supporting_replicates:

  Character, the list of replicates where the edge was detected (only
  added if more than 10% of the edges are not present in all
  replicates).

## Details

The input networks should all contain the same nodes (genes). The output
consensus network contains the same nodes and all edges that were
detected in at least 1 of the replicates.

For each edge, the consensus edge weight (adjacency) is calculated as
the weighted mean of replicate-wise edge weights. The weighted mean
corrects for 1) the phylogenetic distances between species (if the
phylogenetic tree is provided) and 2) the different numbers of
replicates per species. As a result, the approach downweighs the edge
weights of the replicates that 1) belong to closely related species or
2) belong to species with many replicates, so that an imbalanced
sampling across the phylogenetic tree does not bias the consensus
network.

If an edge is not present in one of the replicates, the edge weight in
that replicate is regarded as 0 for the calculation of the weighted
mean. If there are more than 10% such missing/zero-weight edges, the
number of replicates and the names of replicate where an edge was
detected are saved as 2 new edge attributes ("n_supporting_replicates"
and "supporting_replicates", respectively) in the output consensus
network.

In the next steps of the pipeline, this consensus network can be used to
assign modules jointly for all species while avoiding species bias (see
[`assignInitialModules`](https://hellmann-lab.github.io/CroCoNet/reference/assignInitialModules.md)
and
[`pruneModules`](https://hellmann-lab.github.io/CroCoNet/reference/pruneModules.md)).

## Examples

``` r
consensus_network <- createConsensus(network_list, replicate2species, tree)
```
