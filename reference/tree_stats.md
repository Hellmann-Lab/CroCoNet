# Tree-based statistics of the pruned modules

Tree-based statistics per module. The trees were reconstructed based on
module preservation scores (cor.kIM) within and across species. The tips
of resulting tree represent the replicates and the branch lengths
represent the dissimilarity of connectivity patterns between the
replicates. The total tree length is the sum of all branch lengths in
the tree and it measures module variability both within and across
species, while the within-species_diversity is the sum of the human,
gorilla and cynomolgus within-species branch lengths and it measures
module variability within species only. These statistics were calculated
for all possible jackknifed versions of the modules, each of which was
created by removing a target gene assigned to the given module (the
regulators were never excluded). Then the statistics were summarized per
module by taking the median and its 95% confidence interval across all
jackknifed module versions.

## Usage

``` r
tree_stats
```

## Format

A data frame with 7 rows and 10 columns:

- regulator:

  Character, transcriptional regulator.

- module_size:

  Module size, the numer of target genes assigned to a regulator.

- total_tree_length:

  The median of total tree lengths across all jackknifed versions of the
  module.

- var_total_tree_length:

  The variance of total tree lengths across all jackknifed versions of
  the module.

- lwr_total_tree_length:

  The lower bound of the 95% confidence interval of the total tree
  length calculated by jackknifing.

- upr_total_tree_length:

  The upper bound of the 95% confidence interval of the total tree
  length calculated by jackknifing.

- within_species_diversity:

  The median of within-species_diversity values across all jackknifed
  versions of the module.

- var_within_species_diversity:

  The variance of within-species_diversity values across all jackknifed
  versions of the module.

- lwr_within_species_diversity:

  The lower bound of the 95% confidence interval of the
  within-species_diversity calculated by jackknifing.

- upr_within_species_diversity:

  The upper bound of the 95% confidence interval of the
  within-species_diversity calculated by jackknifing.
