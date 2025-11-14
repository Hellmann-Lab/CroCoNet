# Cross-species conservation measures per module

Cross-species conservation measures of the 12 modules in the subsetted
early neuronal differentiation dataset. Modules that were found to be
conserved or diverged overall (across all species) are labelled as
"conserved" or "diverged" in the column `conservation`.

## Usage

``` r
module_conservation_overall
```

## Format

A data frame with 12 rows and 16 columns:

- focus:

  Character, the focus of interest in terms of cross-species
  conservation ("overall").

- regulator:

  Character, transcriptional regulator.

- module_size:

  Integer, the numer of target genes assigned to a regulator.

- total_tree_length:

  Numeric, total tree length per module (the median across all jackknife
  versions of the module).

- lwr_total_tree_length:

  Numeric, the lower bound of the confidence interval of the total tree
  length calculated based on the jackknifed versions of the module.

- upr_total_tree_length:

  Numeric, the upper bound of the confidence interval of the total tree
  length calculated based on the jackknifed versions of the module.

- within_species_diversity:

  Numeric, within-species diveristy per module (the median across all
  jackknife versions of the module).

- lwr_within_species_diversity:

  Numeric, the lower bound of the confidence interval of the
  within-species diversity calculated based on the jackknifed versions
  of the module.

- upr_within_species_diversity:

  Numeric, the upper bound of the confidence interval of the
  within-species diversity calculated based on the jackknifed versions
  of the module.

- fit:

  Numeric, the fitted total tree length at the within-species diversity
  value of the module.

- lwr_fit:

  Numeric, the lower bound of the prediction interval of the fit.

- upr_fit:

  Numeric, the upper bound of the prediction interval of the fit.

- residual:

  Numeric, the residual of the module in the linear model. It is
  calculated as the difference between the observed and expected
  (fitted) total tree lengths.

- weight:

  Numeric, the weight of the module in the linear regression, inversely
  proportional to the variance of total tree lengths.

- t_score:

  Numeric, the t-score of the module. It is calculated as the residual
  normalized by the standard error of the total tree length prediction
  at the given within-species diversity value.

- conservation:

  Character, "not_significant" if the module falls inside the prediction
  interval of the fit, "diverged" if a module has a higher total tree
  length than the upper boundary of the prediction interval, and
  "conserved" if a module has a lower total tree length than the lower
  boundary of the prediction interval.

## Details

To determine whether a module as whole is conserved or diverged overall,
module trees were reconstructed from pairwise preservation scores
between replicates and based on these trees 2 useful statistics were
calculated for each module: the total tree length and the within-species
diversity (for details please see
[`calculatePresStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculatePresStats.md),
[`reconstructTrees`](https://hellmann-lab.github.io/CroCoNet/reference/reconstructTrees.md),
[`calculateTreeStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculateTreeStats.md)).
After fitting a weighted linear model between the total tree length and
within-species diversity values of all modules, a module was considered
diverged if it fell above the prediction interval of the regression
line, while a module was considered conserved if it fell below the
prediction interval of the regression line (for details please see
[`fitTreeStatsLm`](https://hellmann-lab.github.io/CroCoNet/reference/fitTreeStatsLm.md)
and
[`findConservedDivergedModules`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedModules.md).
The degree of conservation/divergence can be further compared between
the modules categorized as conserved/diverged using 2 measures, the
residual and the t-score.
