# Find conserved and diverged modules

Finds conserved and diverged network modules based on a (weighted)
linear model between the total tree length and within-species diversity
(in case the focus of interest is conservation and overall divergence)
or between the subtree length and diversity of a species (in case the
focus of interest is species-specific divergence) of the module trees.
Modules are considered conserved or diverged if they fall outside the
95% prediction interval of the regression line.

## Usage

``` r
findConservedDivergedModules(tree_stats, lm_tree_stats, conf_level = 0.95)
```

## Arguments

- tree_stats:

  Data frame of the tree-based statistics for the pruned modules.

  Columns required in case the focus of interest is conservation and
  overall divergence:

  regulator

  :   Character, transcriptional regulator.

  module_size

  :   Integer, the numer of target genes assigned to a regulator.

  total_tree_length

  :   Numeric, total tree length per module (typically the median across
      all jackknife versions of the module).

  var_total_tree_length

  :   Numeric, the variance of total tree lengths across all jackknifed
      versions of the module (optional, only needed for weighted
      regression).

  lwr_total_tree_length

  :   Numeric, the lower bound of the confidence interval of the total
      tree length calculated based on the jackknifed versions of the
      module (optional, only needed for plotting error bars later on).

  upr_total_tree_length

  :   Numeric, the upper bound of the confidence interval of the total
      tree length calculated based on the jackknifed versions of the
      module (optional, only needed for plotting error bars later on).

  within_species_diversity

  :   Numeric, within-species diveristy per module (typically the median
      across all jackknife versions of the module).

  lwr_within_species_diversity

  :   Numeric, the lower bound of the confidence interval of the
      within-species diversity calculated based on the jackknifed
      versions of the module (optional, only needed for plotting error
      bars later on).

  upr_within_species_diversity

  :   Numeric, the upper bound of the confidence interval of the
      within-species diversity calculated based on the jackknifed
      versions of the module (optional, only needed for plotting error
      bars later on).

  Columns required in case the focus of interest is species-specific
  divergence:

  regulator

  :   Character, transcriptional regulator.

  module_size

  :   Integer, the numer of target genes assigned to a regulator.

  {{species}}\_subtree_length

  :   Numeric, the subtree length of a species per module (typically the
      median across all jackknife versions of the module).

  var\_{{species}}\_subtree_length

  :   Numeric, the variance of the species subtree length across all
      jackknifed versions of the module (optional, only needed for
      weighted regression).

  lwr\_{{species}}\_subtree_length

  :   Numeric, the lower bound of the confidence interval of species
      subtree length calculated based on the jackknifed versions of the
      module (optional, only needed for plotting error bars later on).

  upr\_{{species}}\_subtree_length

  :   Numeric, the upper bound of the confidence interval of species
      subtree length calculated based on the jackknifed versions of the
      module (optional, only needed for plotting error bars later on).

  {{species}}\_diversity

  :   Numeric, the diveristy of the species of interest per module
      (typically the median across all jackknife versions of the module)

  lwr\_{{species}}\_diversity

  :   Numeric, the lower bound of the confidence interval of the
      species-specific diversity calculated based on the jackknifed
      versions of the module (optional, only needed for plotting error
      bars later on).

  upr\_{{species}}\_diversity

  :   Numeric, the upper bound of the confidence interval of the
      species-specific diversity calculated based on the jackknifed
      versions of the module (optional, only needed for plotting error
      bars later on).

- lm_tree_stats:

  An object of class [`lm`](https://rdrr.io/r/stats/lm.html), the
  (weighted) linear model between the total tree length and
  within-species diversity (in case the focus of interest is
  conservation and overall divergence) or between the subtree length and
  diversity of a species (in case the focus of interest is divergence
  between this particular species and all others).

- conf_level:

  Confidence level of the prediction interval (default: 0.95).

## Value

Data frame of cross-species conservation measures per module. In
addition to the relevant columns of the input "tree_stats", it contains
6 or 7 new columns:

- focus:

  Character, the focus of interest in terms of cross-species
  conservation, either "overall" if the focus of interest is
  conservation and overall divergence, or the name of a species if the
  focus of interest is the divergence between this particular species
  and all others. Derived from the input linear model (`lm_tree_stats`).

- fit:

  Numeric, the fitted total tree length/subtree length of a species.

- lwr_fit:

  Numeric, the lower bound of the prediction interval of the fit.

- upr_fit:

  Numeric, the upper bound of the prediction interval of the fit.

- residual:

  Numeric, the residual of the module in the linear model. It is
  calculated as the difference between the observed and expected
  (fitted) total tree length/subtree length of a species.

- weight:

  Numeric, the weight of the module in the linear regression, inversely
  proportional to the variance of the total tree length/the subtree
  length of a species (only if `weighted_lm` is set to TRUE).

- t_score:

  Numeric, the t-score of the module. It is calculated as the residual
  normalized by the standard error of the total tree length/species
  subtree length prediction.

- conservation:

  Character, "not_significant" if the module falls inside the prediction
  interval of the fit, "diverged" if a module has a higher total tree
  length/subtree length of a species than the upper boundary of the
  prediction interval, and "conserved" if a module has a lower total
  tree length than the lower boundary of the prediction interval.

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

Various useful statistics can be defined based on these trees (see also
[`calculateTreeStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculateTreeStats.md)).
The total tree length is the sum of all branch lengths in the tree and
it measures module variability both within and across species. The
diversity of a species is the sum of the branches connecting the
replicates of this species and it measures module variability within
this particular species. The within-species diversity is the sum of the
diversity values across all species and it measures module variability
within species in general. The subtree length of a species is defined as
the sum of the branch lengths in the subtree that is defined by the
replicates of the species and includes the internal branch connecting
these replicates to the rest of the tree, and it measures module
variability within this species as well as between this species and all
others.

Overall divergence is best represented by the total tree length, while
species-specific divergence is best represented by the subtree length of
a species, however the absolute values of these statistics in themselves
are difficult to interpret. The diversity (either overall or for a
particular species) can be used as an internal reference to determine
whether the total tree length or the subtree length of a species is
larger/smaller than expected. The general relationship between total
tree length and overall diversity or between the subtree length and the
diversity of the species is captured via linear regression (see also
[`fitTreeStatsLm`](https://hellmann-lab.github.io/CroCoNet/reference/fitTreeStatsLm.md)).
In term of cross-species conservation, the most interesting modules are
then the outliers that do not follow the linear trend.

By analyzing the linear regression between total tree length and
within-species diversity, we can single out modules that are conserved
or diverged overall. The modules that fall far below the trend line are
the most conserved ones: they have a lower total tree length than
expected based on their within-species diversity, meaning that the
species are not well-separated within the tree but rather mixed among
each other. In contrast, the modules that are located far above the
trend line are the most diverged ones: they have a higher total tree
length than expected based on their within-species diversity, meaning
that they thave long cross-species branches between one or more pairs of
species.

By analyzing the linear regression between the subtree length and the
diversity of a species, we can single out modules that are diverged
between this particular species and all others: these are the modules
that fall far above the trend line, i.e. that have a longer branch
leading to this species than expected based on the diversity of the
species. The modules that fall below the trend line are in this case not
meaningful, since we can only call a module conserved, if it is
conserved across all species.

For both types of analysis, we quantitatively define which modules are
outliers, and thus conserved or diverged, by calculating the prediction
interval of the linear fit. A module is considered diverged overall/in a
species-specific manner if it has a higher total tree length/subtree
length of a species than the upper boundary of the prediction interval,
while a module is considered conserved if it has a lower total tree
length than the lower boundary of the prediction interval.

The degree of conservation/divergence can be further compared between
the modules categorized as conserved/diverged. We recommend 2
measures: 1) the residual, which is the difference between the observed
and expected total tree length/subtree length of a species, and 2) the
t-score, which is the residual normalized by the standard error of the
total tree length/species subtree length prediction.

The linear regression can be weighted by the error of the data points.
To gain this information, we can use jackknifing: each member gene of a
module is removed and all statistics are recalculated (see the parameter
`jackknife` in the function
[`calculatePresStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculatePresStats.md)).
The weight of a module in the regression is then defined to be inversely
proportional to the variance of the dependent variable (total tree
length/subtree length of a species) across all jackknife versions.

## Examples

``` r
module_conservation_overall <- findConservedDivergedModules(tree_stats, lm_overall)
```
