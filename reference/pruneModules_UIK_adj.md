# Prune modules based on the regulator-target adjacencies using dynamic filtering

Prunes the initial modules by applying a dynamic stepwise pruning based
on the regulator-target adjacencies.

## Usage

``` r
pruneModules_UIK_adj(
  initial_modules,
  min_module_size = 20L,
  max_frac_modules_lost = 0,
  exponent = 1L
)
```

## Arguments

- initial_modules:

  Data frame of initial modules, required columns:

  regulator

  :   Character, transcriptional regulator.

  target

  :   Character, member gene of the regulator's initial module.

  weight

  :   Numeric, consensus edge weight/adjacency, the weighted mean of
      replicate-wise edge weights.

- min_module_size:

  Integer, the lower threshold of module size. Modules with a smaller
  size than this threshold are removed after each pruning step (default:
  20).

- max_frac_modules_lost:

  Numeric, the threshold for the fraction of removed modules (default:
  0).

- exponent:

  Integer, the exponent the regulator-target adjacency and intramodular
  connectivity is raised to the power of during the cumulative sum curve
  calculation (default: 1, i.e. the regulator-target adjacencies and
  intramodular connectivities stay unchanged).

## Value

Data frame of the pruned modules with the following columns:

- regulator:

  Character, transcriptional regulator.

- module_size:

  Integer, the numer of genes assigned to a regulator.

- target:

  Character, target gene of the transcriptional regulator (member of the
  regulator's pruned module).

- weight:

  Numeric, consensus edge weight/adjacency, the weighted mean of
  replicate-wise edge weights.

Additional columns present in `initial_modules` will also be preserved
in `pruned_modules`.

## Details

For each module, the initial module members are filtered in successive
steps based on their regulator-target edge weight/adjacency, which
quantifies how strongly a target is connected to the regulator. In each
step, the cumulative sum curve of the regulator-target adjacencies is
calculated per module, the knee point of the curve is identified using
the Unit Invariant Knee (UIK) method (see `uik`), then only the targets
that rank higher than the knee point are kept.

The modules containing less target genes than `min_module_size` are
removed after each pruning step. The steps continue until the fraction
of removed modules becomes higher than `max_frac_modules_lost`. By
default, this cutoff is set to zero, meaning that all modules need to
stay above the specified `min_module_size`. It is recommended to set
`min_module_size` to at least 20, because the correlation-based
preservation statistics in the next steps might be coupled with high
uncertainty for modules smaller than this (see
[`calculatePresStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculatePresStats.md)).

If `exponent` is set higher than 1, the adjacencies are raised to the
equivalent power when calculating the cumulative sum curves and their
knee points.

While setting the parameter `min_module_size` prevents the modules from
becoming too small, the exact number of target genes per regulator does
not have to be pre-defined, in line with the notion that different
regulators can have an effect on different numbers of genes. There is
also no hard cutoff applied to the regulator-target adjacencies, but by
using knee point detection the target genes are filtered in a
data-driven way.

The modules are allowed to overlap, and in addition to having its own
module, a regulator can be assigned to another regulator's module as
well, in line with the notion that genes can be multifunctional and gene
regulation can be combinatorial.

## References

Christopoulos, D. (2016). Introducing Unit Invariant Knee (UIK) As an
Objective Choice for Elbow Point in Multivariate Data Analysis
Techniques. SSRN Electronic Journal.
https://doi.org/10.2139/SSRN.3043076

## See also

Other methods to prune modules:
[`pruneModules_UIK_adj_kIM()`](https://hellmann-lab.github.io/CroCoNet/reference/pruneModules_UIK_adj_kIM.md),
[`pruneModules_topN()`](https://hellmann-lab.github.io/CroCoNet/reference/pruneModules_topN.md)

## Examples

``` r
pruned_modules_UIK_adj <- pruneModules_UIK_adj(initial_modules)
#> Step 1: filtering targets based on their adjacencies to the regulator
#> Median module size after filtering: 92.5
#> Step 2: filtering targets based on their adjacencies to the regulator
#> Median module size after filtering: 35
```
