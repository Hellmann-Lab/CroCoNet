# Tree-based statistics of the original and jackknifed pruned modules

Tree-based statistics calculated for the original and all jackknifed
versions of the pruned modules. For each module version, a tree was
reconstructed based on the module preservation scores within and across
species. The tips of resulting tree represent the replicates and the
branch lengths represent the dissimilarity of connectivity patterns
between the replicates. Branches between replicates of different species
carry information about cross-species differences, while the branches
between replicates of the same species carry information about the
within-species diversity.

## Usage

``` r
tree_stats_jk
```

## Format

A data frame with 232 rows and 16 columns:

- regulator:

  Character, transcriptional regulator.

- module_size:

  Module size, the numer of target genes assigned to a regulator.

- type:

  Character, module type (orig = original or jk = jackknifed).

- id:

  Character, the unique ID of the module version (format:
  nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and
  nameOfRegulator_orig in case of module type 'orig').

- gene_removed:

  The name of the gene removed by jackknifing (NA in case of module type
  'orig').

- within_species_diversity:

  Numeric, the sum of the human, gorilla and cynomolgus diversity.

- total_tree_length:

  Numeric, the sum of all branch lengths in the tree.

- human_diversity:

  Numeric, the sum of the branch lengths in the subtree that contains
  only the human tips.

- gorilla_diversity:

  Numeric, the sum of the branch lengths in the subtree that contains
  only the gorilla tips.

- cynomolgus_diversity:

  Numeric, the sum of the branch lengths in the subtree that contains
  only the cynomolgus tips.

- human_monophyl:

  Logical indicating whether the tree is monophyletic for the human
  replicates.

- gorilla_monophyl:

  Logical indicating whether the tree is monophyletic for the gorilla
  replicates.

- cynomolgus_monophyl:

  Logical indicating whether the tree is monophyletic for the cynomolgus
  replicates.

- human_subtree_length:

  Numeric, the sum of the branch lengths in the subtree that is defined
  by the human replicates and includes the internal branch connecting
  the human replicates to the rest of the tree. NA if the tree is not
  monophyletic for the human replicates.

- gorilla_subtree_length:

  Numeric, the sum of the branch lengths in the subtree that is defined
  by the gorilla replicates and includes the internal branch connecting
  the gorilla replicates to the rest of the tree. NA if the tree is not
  monophyletic for the gorilla replicates.

- cynomolgus_subtree_length:

  Numeric, the sum of the branch lengths in the subtree that is defined
  by the cynomolgus replicates and includes the internal branch
  connecting the cynomolgus replicates to the rest of the tree. NA if
  the tree is not monophyletic for the cynomolgus replicates.
