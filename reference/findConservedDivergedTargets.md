# Find conserved and diverged target genes within a module

Calculates a contribution score for each target gene in a
conserved/diverged module and finds genes that are particularly
responsible for the conservation/divergence.

## Usage

``` r
findConservedDivergedTargets(
  module_name,
  tree_stats_jk,
  lm_tree_stats,
  conf_level = 0.95
)
```

## Arguments

- module_name:

  Character, the name of the module of interest.

- tree_stats_jk:

  Data frame of tree statistics per jackknife module version across all
  modules.

  Columns required in case the focus of interest is conservation and
  overall divergence:

  regulator

  :   Character, transcriptional regulator.

  module_size

  :   Integer, the numer of target genes assigned to a regulator.

  type

  :   Character, module type (orig = original or jk = jackknifed).

  id

  :   Character, the unique ID of the module version (format:
      nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk'
      and nameOfRegulator_orig in case of module type 'orig').

  gene_removed

  :   Character, the name of the gene removed by jackknifing (NA in case
      of module type 'orig').

  total_tree_length

  :   Numeric, total tree length per jackknife module version.

  within_species_diversity

  :   Numeric, within-species diversity per jackknife module version.

  Columns required in case the focus of interest is species-specific
  divergence:

  regulator

  :   Character, transcriptional regulator.

  module_size

  :   Integer, the numer of target genes assigned to a regulator.

  type

  :   Character, module type (orig = original or jk = jackknifed).

  id

  :   Character, the unique ID of the module version (format:
      nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk'
      and nameOfRegulator_orig in case of module type 'orig').

  gene_removed

  :   Character, the name of the gene removed by jackknifing (NA in case
      of module type 'orig').

  {{species}}\_subtree_length

  :   Numeric, the subtree length of a species per jackknife module
      version.

  {{species}}\_diversity

  :   Numeric, diversity within the species of interest per jackknife
      module version.

- lm_tree_stats:

  Object of class [`lm`](https://rdrr.io/r/stats/lm.html), the
  (weighted) linear model fit between the total tree length and
  within-species diversity (in case the focus of interest is
  conservation and overall divergence) or between the subtree length and
  diversity of a species (in case the focus of interest is the
  divergence between that particular species and all others).

- conf_level:

  Confidence level of the prediction interval (default: 0.95).

## Value

Data frame of cross-species conservation measures per target gene for
the module of interest. In addition to the relevant columns of
`tree_stats_jk`, it contains 6 new columns:

- focus:

  Character, the focus of interest in terms of cross-species
  conservation, either "overall" if the focus of interest is
  conservation and overall divergence, or the name of a species if the
  focus of interest is the divergence between that particular species
  and all others.

- fit:

  Numeric, the fitted total tree length/subtree length of a species at
  the diversity value of a jackknife module version.

- lwr_fit:

  Numeric, the lower bound of the prediction interval of the fit.

- upr_fit:

  Numeric, the upper bound of the prediction interval of the fit.

- residual:

  Numeric, the residual of the jackknife module version in the linear
  model. It is calculated as the difference between the observed and
  expected (fitted) total tree length/subtree length of a species.

- contribution:

  Numeric, the contribution score of the removed target gene.

It also contains a new summary row with the aggregate values across all
jackknife values of the module.

## Details

To determine whether a module as whole is conserved, diverged overall or
diverged on a specific lineage, the CroConet approach relies on module
trees reconstrcuted from pairwise preservation scores between replicates
and statistics calculated based on these trees (total tree length,
within-species diversity, subtree length and diversity of a particular
species). To identify individual genes that contribute the most to
conservation/divergence, the same statistics can be used in combination
with jackknifing.

During jackknifing, each member gene of a module is removed and all
statistics are recalculated, including the tree-based statistics that
inform us about cross-species conservation (see the parameter
`jackknife` in the function
[`calculatePresStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculatePresStats.md)).
Our working hypothesis is that if removing a target from a diverged
module makes the module more conserved, then that target was responsible
for divergence in the original module, and vice versa, if removing a
target from a conserved module makes the module more diverged, then that
target was responsible for conservation in the original module.

To quantify these effects, the function takes use of the linear models
that were fitted across all modules between the total tree length and
within-species diversity (in case the focus of interest is conservation
and overall divergence) or between the subtree length and diversity of a
species (in case the focus of interest is species-specific divergence),
similarly to how the conserved and diverged modules were identified in
the first place. However, in this case it is not the aggregate statistic
for a module as a whole that is compared to the regression line, but the
statistic of each jackknife module version separately. Using the
residual of a jackknifed version (\\\hat{\epsilon}\_i\\) and of the
original module (\\\hat{\epsilon}\_\text{orig}\\), we defined the target
gene contribution score (\\c_i\\) as follows: \$\$c_i =
\frac{\hat{\epsilon}\_\text{orig} -
\hat{\epsilon}\_i}{\hat{\epsilon}\_\text{orig}}\$\$

This score is applicable both to conserved and to diverged modules: the
highest score corresponds to the most diverged target gene in case of a
diverged module and to the most conserved target gene in case of a
conserved module. A score of zero means that removing a target gene does
not affect the signal of conservation or divergence, i.e. the extent of
cross-species differences associated with this gene reflect the module
average. A negative score indicates that removing the target gene
strengthens the signal of conservation or divergence, i.e. the target
acts against the overall module trend. A positive score means that
removing a target gene diminishes the signal of conservation or
divergence, i.e. the target contributes to the observed signal when it
is present. A score of 1 (which corresponds to \\\hat{\epsilon}\_i =
0\\) means that the removal of the target gene completely eliminates the
signal, therefore it is single-handedly responsible for the observed
conservation or divergence. Finally, a score greater than 1 suggests
that removing the gene reverses the direction of the signal, meaning
that it is a strongly diverged gene within an otherwise conserved
module, or vice versa.

It is important to note, that in most cases the divergence/conservation
of a module cannot be contributed to a single target gene, but it is
rather the combined signal of all targets together. Therefore
contribution scores nearing or exceeding 1 are rare. This also means
that the "most conserved" target gene of diverged module is still most
likely diverged, just to a lesser extent than the others. Thus, it
mainly makes sense to investigate the most diverged targets of the
diverged modules and the most conserved targets of the conserved
modules, not the other way around.

Based on the linear model provided, the function calculates the fitted
values of the total tree length/species subtree length, the prediction
interval of the fit, the residual of each jackknife module version and
the original module, as well as the target gene contribution scores. The
results can be plotted using the function
[`plotConservedDivergedTargets`](https://hellmann-lab.github.io/CroCoNet/reference/plotConservedDivergedTargets.md).

## Examples

``` r
POU5F1_target_conservation <- findConservedDivergedTargets("POU5F1", tree_stats_jk, lm_overall)
```
