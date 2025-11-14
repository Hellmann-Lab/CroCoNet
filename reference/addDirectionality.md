# Add directionality of interactions

Determines the directionality of each network edge (positively
correlated/coexpressed or negatively correlated/anti-coexpressed) based
on a modified Spearman's correlation between the expression profiles of
the 2 genes that form the edge.

## Usage

``` r
addDirectionality(network, sce, assay = "logcounts", n_cores = 1L)
```

## Arguments

- network:

  An `igraph` object containing the consensus network or the network of
  a replicate.

- sce:

  A `SingleCellExperiment` object containing the expression data either
  for all replicates (in case `network` is the consensus network) or for
  the replicate of interest (in case `network` is a replicate-wise
  network). If `sce` contains the expression data of all replicates, it
  is also expected to have a metadata column "replicate" specifying
  which replicate each cell belongs to; this will be used to define the
  blocking levels.

- assay:

  Character, the name of the assay in `sce` that should be used for the
  calculation of gene-gene correlations (default: "logcounts").

- n_cores:

  Integer, the number of cores (default: 1).

## Value

An `igraph` object, the input `network` extended by the information
about the direction of interactions. In addition to the original edge
attributes, it contains 3 new attributes:

- rho:

  Numeric, the approximate Spearman's correlation coefficient between
  the expression profiles of the 2 genes that form the edge.

- p.adj:

  Numeric, BH-corrected approximate p-value of rho.

- direction:

  Character, the direction of the interaction between the 2 genes that
  form the edge ("+" = positively correlated/coexpressed or "-" =
  negatively correlated/anti-coexpressed).

## Details

If the networks were inferred using a method that does not distinguish
coexpressed and anti-coexpressed gene pairs, it might be useful to add
this information for lines of analysis where it makes sense to separate
the activated and repressed target genes of a regulator (e.g. for the
calculation of eigengenes, see
[`calculateEigengenes`](https://hellmann-lab.github.io/CroCoNet/reference/calculateEigengenes.md)).
If the network inference method output edges with both positive and
negative edge weights in the first place (e.g. correlation-based
methods), the edge attribute "direction" is already created during the
step
[`normalizeEdgeWeights`](https://hellmann-lab.github.io/CroCoNet/reference/normalizeEdgeWeights.md)
and does not have to be calculated again.

Here the directionality of a geneA-geneB edge refers to the
characteristic whether geneA and geneB are coexpressed or
anti-coexpressed and NOT whether geneA regulates geneB or geneB
regulates geneA. The network remains undirected in a graph theoretical
sense.

The calculation of directionality relies on the approximate version of
the Spearman's rho, significance testing and blocking implemented by
`correlatePairs`. The results are summarized as 3 new edge attributes in
the `igraph` object: rho (approximate Spearman's correlation
coefficient), p.adj (BH-corrected approximate p-value) and direction
("+" or "-").

## Examples

``` r
consensus_network <- network_list %>%
 createConsensus(replicate2species, tree) %>%
 addDirectionality(sce)
#> Loading required namespace: SingleCellExperiment
```
