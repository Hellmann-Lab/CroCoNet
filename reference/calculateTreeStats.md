# Calculate tree statistics across all modules

Calculates various tree statistics (total tree length, within-species
diversity, diversity, monophyleticity and subtree length of each
species) for all modules/jackknifed module versions based on a list of
tree representations.

## Usage

``` r
calculateTreeStats(tree_list, n_cores = 1L)
```

## Arguments

- tree_list:

  A named list of `phylo` objects containing the tree representations of
  all modules/jackknifed module versions, with the tips of the trees
  corresponding to replicates. All trees are expected to have a
  component `species` that specifies which species each tip belongs to.
  The trees can also contain a component `info` that stores metadata of
  the module/jackknifed module version in a data frame format.

- n_cores:

  Integer, the number of cores (default: 1).

## Value

A data frame of tree statistics with the following columns:

- regulator:

  Character, transcriptional regulator.

- module_size:

  Module size, the number of target genes assigned to a regulator (only
  present if the column is also present in the component `info` of the
  trees).

- type:

  Character, module type (orig = original or jk = jackknifed, only
  present if the input trees were reconstructed with jackknifing).

- id:

  Character, the unique ID of the module version (format:
  nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and
  nameOfRegulator_orig in case of module type 'orig', only present if
  the input trees were reconstructed with jackknifing).

- gene_removed:

  Character, the name of the gene removed by jackknifing (NA in case of
  module type 'orig', only present if the input trees were reconstructed
  with jackknifing).

- within_species_diversity:

  Numeric, the sum of the species-wise diversities.

- total_tree_length:

  Numeric, the total length of all branches in the tree.

- {{species}}\_diversity:

  Numeric, as many columns as there are species, each of them containing
  the the total length of the branches connecting the replicates of the
  given species.

- {{species}}\_monophyl:

  Logical, as many columns as there are species, each of them indicating
  whether the tree is monophyletic for the replicates of the given
  species.

- {{species}}\_subtree_length:

  Numeric, as many columns as there are species, each of them containing
  the sum of the branch lengths in the subtree that is defined by the
  replicates of the species and includes the internal branch connecting
  these replicates to the rest of the tree. NA if the tree is not
  monophyletic for the replicates of the given species.

## Details

As part of the CroCoNet approach, pairwise module preservation scores
are calculated between replicates, both within and across species (see
[`calculatePresStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculatePresStats.md))
and neighbor-joining trees are reconstructed based on these preservation
scores per module (see
[`convertPresToDist`](https://hellmann-lab.github.io/CroCoNet/reference/convertPresToDist.md)
and
[`reconstructTrees`](https://hellmann-lab.github.io/CroCoNet/reference/reconstructTrees.md)).
The tips of the resulting tree represent the replicates and the branch
lengths represent the dissimilarity of module connectivity patterns
between the networks of 2 replicates.

Various useful statistics can be defined based on these module trees:

- Total tree length: The sum of all branch lengths in the tree; measures
  module variability both within and across species.

- Diversity of a species: The total length of the branches connecting
  the replicates of the given species to each other; measures module
  variability within this particular species.

- Within-species diversity: The sum of the diversity values across all
  species; measures module variability within species in general.

- Monophyleticity of a species: Indicates whether the tree is
  monophyletic for the replicates of the given species. Only if a module
  tree is monophyletic for a species of interest can the module be
  tested for divergence between this species and all others.

- Subtree length of a species: The sum of the branch lengths in the
  subtree that is defined by the replicates of the species and includes
  the internal branch connecting these replicates to the rest of the
  tree. Undefined if the tree is not monophyletic for the species of
  interest.

In the later steps of the pipeline, these tree-based statistics can be
used to 1) identify modules that are conserved, diverged overall or
diverged between a species and all others (see
[`findConservedDivergedModules`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedModules.md)),
and 2) pinpoint individual target genes within these modules that
contribute the most to the conservation/divergence (see
[`findConservedDivergedTargets`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedTargets.md)).

## Examples

``` r
tree_stats_jk <- calculateTreeStats(trees_jk)
```
