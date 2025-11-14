# Fit (weighted) linear model between tree-based statistics

Fits a (weighted) linear model between the total tree length and
within-species diversity (in case the focus of interest is conservation
and overall divergence) or between the subtree length and the diversity
of a species (in case the focus of interest is the species-specific
divergence).

## Usage

``` r
fitTreeStatsLm(tree_stats, focus = "overall", weighted_lm = TRUE)
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

  within_species_diversity

  :   Numeric, within-species diveristy per module (typically the median
      across all jackknife versions of the module).

  Columns required in case the focus of interest is species-specific
  divergence:

  regulator

  :   Character, transcriptional regulator.

  module_size

  :   Integer, the numer of target genes assigned to a regulator.

  {{species}}\_subtree_length

  :   Numeric, the sum of the branch lengths in the subtree that is
      defined by the replicates of the species and includes the internal
      branch connecting these replicates to the rest of the tree
      (typically the median across all jackknife versions of the
      module).

  var\_{{species}}\_subtree_length

  :   Numeric, the variance of the subtree length of a species across
      all jackknifed versions of the module (optional, only needed for
      weighted regression).

  {{species}}\_diversity

  :   Numeric, within-species diveristy per module (typically the median
      across all jackknife versions of the module).

- focus:

  Character, the focus of interest in terms of cross-species
  conservation, either "overall" if the focus of interest is
  conservation and overall divergence, or the name of a species if the
  focus of interest is the divergence between that particular species
  and all others.

- weighted_lm:

  Logical indicating whether the linear regression should be weighted or
  not (default: TRUE). If TRUE, `tree_stats` is expected to contain the
  column `var_total_tree_length` and the weights will be inversely
  proportional to these variances. If no jackknifing was performed and
  thus these variances were not calculated, please set this parameter to
  FALSE.

## Value

An object of class [`lm`](https://rdrr.io/r/stats/lm.html).

## Details

The linear models output by this function can be used to identify
conserved and diverged modules, and to identify target genes within
these modules that contribute the most to the conservation/divergence.
For details, please see
[`findConservedDivergedModules`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedModules.md)
and
[`findConservedDivergedTargets`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedTargets.md).

The focus of interest can be specified using the parameter `focus`. If
`focus` is set to "overall" (default), the linear model will be fit
between the total tree length and within-species diversity, and
subsequent analysis using
[`findConservedDivergedModules`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedModules.md)
and
[`findConservedDivergedTargets`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedTargets.md)
can identify modules and target genes that are conserved or diverged
across all species. If `focus` is set to the name of a species in the
dataset, the linear model will be fit between the subtree length and the
diversity of that species, and subsequent analysis using
[`findConservedDivergedModules`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedModules.md)
and
[`findConservedDivergedTargets`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedTargets.md)
can identify modules and target genes that are diverged between the
species and all others. Please note that if the aim is to find conserved
modules, `focus` should always be set to "overall".

The function fits the linear model corresponding to the focus of
interest by calling [`lm`](https://rdrr.io/r/stats/lm.html). If a
weighted model is desired (default), the weights are defined to be
inversely proportional to the variance of the dependent variable (total
tree length or the subtree length of a species). If no jackknifing was
performed and thus the variance is unknown, please set `weighted_lm` to
FALSE.

## Examples

``` r
lm_overall <- fitTreeStatsLm(tree_stats, focus = "overall")
```
