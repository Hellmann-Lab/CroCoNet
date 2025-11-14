# Reconstruct trees across all modules

Reconstructs neighbor-joining trees based on pairwise distance measures
between replicates.

## Usage

``` r
reconstructTrees(dist_list, n_cores = 1L)
```

## Arguments

- dist_list:

  A named list containing the distance measures as data frames for all
  modules or jackknifed module versions, the output of
  [`convertPresToDist`](https://hellmann-lab.github.io/CroCoNet/reference/convertPresToDist.md).
  Required columns for the data frames:

  regulator

  :   Character, transcriptional regulator.

  module_size

  :   Integer, the numer of target genes assigned to a regulator.

  type

  :   Character, module type ("orig" = original or "jk" = jackknifed,
      only needed if the distances were calculated with jackknifing).

  id

  :   Character, the unique ID of the module version (format:
      nameOfRegulator_jk_nameOfGeneRemoved in case of module type "jk"
      and nameOfRegulator_orig in case of module type "orig", only
      needed if the distances were calculated with jackknifing).

  gene_removed

  :   Character, the name of the gene removed by jackknifing (NA in case
      of module type "orig", only needed if the distances were
      calculated with jackknifing).

  replicate1, replicate2

  :   Character the names of the replicates compared.

  species1, species2

  :   Character, the names of the species `replicate1` and `replicate2`
      belong to, respectively.

  dist

  :   Numeric, the distance measure to be used for the tree
      reconstruction.

- n_cores:

  Number of cores.

## Value

A named list containing the neighbor-joining trees as `phylo` objects
for all modules or jackknifed module versions in the input `dist_list`.
All trees contain a component `species` that specifies which species
each tip belongs to and a component `info` that stores metadata of the
module/jackknifed module version in a data frame format.

## Details

As part of the CroCoNet approach, pairwise module preservation scores
are calculated between replicates, both within and across species, to
quantify how similar module connectivity patterns are between the
networks of two replicates (see
[`calculatePresStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculatePresStats.md)).
These preservation scores are then converted into distance measures (see
[`convertPresToDist`](https://hellmann-lab.github.io/CroCoNet/reference/convertPresToDist.md)).

This function first sorts the distance measures of each
module/jackknifed module version into a distance matrix of all
replicates, then based on this distance matrix reconstructs a tree using
the neighbor-joining algorithm.

The procedure results in a single tree per module/jackknifed module
version (depending on whether the preservation statistics were
calculated with or without jackknifing, see see the parameter
`jackknife` in the function
[`calculatePresStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculatePresStats.md)).
The tips of the trees represent the replicates and the branch lengths
represent the dissimilarity of module topology between the networks of 2
replicates. The trees are output as a list of `phylo` objects.

In the next steps of the pipeline, statistics based on these trees can
be used to identify conserved and diverged modules and pinpoint target
genes within these modules that contribute the most to
conservation/divergence (see
[`calculateTreeStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculateTreeStats.md),
[`fitTreeStatsLm`](https://hellmann-lab.github.io/CroCoNet/reference/fitTreeStatsLm.md),
[`findConservedDivergedModules`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedModules.md)
and
[`findConservedDivergedTargets`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedTargets.md)).

## Examples

``` r
trees_jk <- reconstructTrees(dist_jk)
```
