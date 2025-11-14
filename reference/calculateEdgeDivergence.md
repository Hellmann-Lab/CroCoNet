# Calculate the divergence of intramodular edges

Selects the top `N` strongest edges of the chosen module(s) based on the
consensus network and calculates a divergence score for each edge.

## Usage

``` r
calculateEdgeDivergence(
  module_names,
  pruned_modules,
  consensus_network,
  network_list,
  replicate2species,
  N = Inf,
  n_cores = 1L
)
```

## Arguments

- module_names:

  Character or character vector, the name(s) of the module(s) of
  interest.

- pruned_modules:

  Data frame of the pruned modules, required columns:

  regulator

  :   Character, transcriptional regulator.

  target

  :   Character, target gene of the transcriptional regulator (member of
      the regulator's pruned module).

- consensus_network:

  `igraph` object, the consensus network across all species and
  replicates.

- network_list:

  A list of `igraph` objects containing the networks per replicate.

- replicate2species:

  A data frame that specifies which species each replicate belongs to,
  required columns:

  replicate

  :   Character, name of the replicate.

  species

  :   Character, name of the species.

- N:

  Integer, the number of strongest edges to subset. If set to Inf
  (default), all edges in the module are considered.

- n_cores:

  Integer, the number of cores (default: 1).

## Value

A data frame of the selected edges with 5 columns:

- regulator:

  Character, transcriptional regulator.

- from, to:

  Character, the 2 member genes of the regulator's module that form the
  edge.

- consensus_weight:

  Numeric, consensus edge weight/adjacency (the weighted average of
  replicate-wise adjacencies).

- f_statistic:

  Numeric, measue of edge divergence. It is calculated as the
  F-statistic from the ANOVA of edge weights with species as groups.

- p-value:

  Numeric, the p-value of the F-statistic.

## Details

In the CroCoNet approach, networks are reconstructed per replicate and
combined into a single phylogeny-aware consensus network which is the
basis of the module assignment.

This function selects the intramodular edges of the input module(s) in
the consensus network. If `N` is set to `Inf`, all intramodular edges
are considered for the divergence calculation, if `N` is smaller than
the module size, the edges are ordered by their consensus edge weight
and only the top `N` edges are kept per module.

For each edge that was kept, an edge divergence score is calculated
based on its edge weights in the networks of individual replicates. The
edge weights are compared across species using an ANOVA. and the
F-statistic (i.e. the variation across the species means / variation
within the species) and the p-value of this F-statistic are output as
the measures of edge divergence.

## Examples

``` r
POU5F1_mod_edge_divergence <- calculateEdgeDivergence("POU5F1",
                                                                 pruned_modules,
                                                                 consensus_network,
                                                                 network_list,
                                                                 replicate2species)
```
