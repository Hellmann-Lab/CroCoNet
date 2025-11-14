# Cross-species conservation measures of the target genes in the POU5F1 module

Cross-species conservation measures of the 21 target genes in the POU5F1
module in the subsetted early neuronal differentiation dataset.

## Usage

``` r
POU5F1_target_conservation
```

## Format

A data frame with 23 rows and 12 columns:

- focus:

  Character, the focus of interest in terms of cross-species
  conservation ("overall").

- regulator:

  Character, transcriptional regulator.

- module_size:

  Integer, the numer of target genes assigned to a regulator.

- type:

  Character, module type (orig = original, jk = jackknifed or summary =
  the summary across all jackknife module versions).

- id:

  Character, the unique ID of the module version (format:
  nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk',
  nameOfRegulator_orig in case of module type 'orig' and
  nameOfRegulator_summary in case of module type 'summary').

- gene_removed:

  Character, the name of the gene removed by jackknifing (NA in case of
  module types 'orig' and 'summary').

- total_tree_length:

  Numeric, total tree length per jackknife module version.

- within_species_diversity:

  Numeric, within-species diversity per jackknife module version.

- fit:

  Numeric, the fitted total tree length at the within-species diversity
  value of a jackknife module version.

- lwr_fit:

  Numeric, the lower bound of the prediction interval of the fit.

- upr_fit:

  Numeric, the upper bound of the prediction interval of the fit.

- residual:

  Numeric, the residual of the jackknife module version in the linear
  model. It is calculated as the difference between the observed and
  expected (fitted) total tree lengths.

- contribution:

  Numeric, the contribution score of the removed target gene.

## Details

To determine which target genes are the most responsible for the
divergence of the POU5F1 module, the same statistics were be used as for
the identification of conserved/diverged modules but in combination with
jackknifing. This means that each target gene in the module was removed
and all statistics, including the tree-based statistics that carry
information about cross-species conservation, were re-calculated (please
see
[`calculatePresStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculatePresStats.md),
[`reconstructTrees`](https://hellmann-lab.github.io/CroCoNet/reference/reconstructTrees.md)
and
[`calculateTreeStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculateTreeStats.md)).
Our working hypothesis is that if removing a target gene from the
diverged POU5F1 module makes the module more conserved, then this gene
contributed to the divergence of the original module in the first place.
This can be quantified using the residuals compared to the regression
line between the total tree lengths and within-species diversities of
the tree representations across all modules (please see
[`fitTreeStatsLm`](https://hellmann-lab.github.io/CroCoNet/reference/fitTreeStatsLm.md)
and
[`findConservedDivergedTargets`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedTargets.md).
The jackknifed versions of the POU5F1 module that have the lowest
residuals are the most conserved, therefore the corresponding target
genes are considered the most diverged.
