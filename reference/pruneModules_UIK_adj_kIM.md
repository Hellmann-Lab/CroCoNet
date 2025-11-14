# Prune modules based on the regulator-target adjacencies and intramodular connectivities using dynamic filtering

Prunes the initial modules by applying a dynamic stepwise pruning based
on the regulator-target adjacencies and intramodular connectivities.

## Usage

``` r
pruneModules_UIK_adj_kIM(
  initial_modules,
  consensus_network,
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

- consensus_network:

  `igraph` object, the consensus network across all species and
  replicates.

- min_module_size:

  Integer, the lower treshold of module size. Modules with a smaller
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
steps based on 2 metrics alternately: 1) the regulator-target edge
weight/adjacency, which quantifies how strongly a target is connected to
the regulator, and 2) the intramodular connectivity, which quantifies
how strongly a target is connected to all other genes in the module. In
each step, the cumulative sum curve based on one of these two
characteristics is calculated per module, the knee point of the curve is
identified using the Unit Invariant Knee (UIK) method (see `uik`), then
only the targets that rank higher than the knee point are kept. The
first step is based on the regulator-target adjacencies, the second step
is based on the intramodular connectivities, the third step is again
based on the regulator-target adjacencies and so on. The intramodular
connectivity is recalculated in each relevant filtering step based on
the then-current module assignment.

The modules containing less target genes than `min_module_size` are
removed after each pruning step. The steps continue until the fraction
of removed modules becomes higher than `max_frac_modules_lost`. By
default, this cutoff is set to zero, meaning that all modules need to
stay above the specified `min_module_size`. It is recommended to set
`min_module_size` to at least 20, because the correlation-based
preservation statistics in the next steps might be coupled with high
uncertainty for modules smaller than this (see
[`calculatePresStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculatePresStats.md)).

If `exponent` is set higher than 1, the regulator-target adjacencies and
intramodular connectivities are raised to the equivalent power when
calculating the cumulative sum curves and their knee points.

Pruning based on intramodular connectivities in addition to the
regulator-target adjacencies ensures that the chosen targets co-vary not
just with the main regulator but also with the rest of the module. These
intramodular connections between targets can carry important information
about combinatorial regulation, feedback loops and co-functionality.

While setting the parameter `min_module_size` prevents the modules from
becoming too small, the exact number of target genes per regulator does
not have to be pre-defined, in line with the notion that different
regulators can have an effect on different numbers of genes. There are
also no hard cutoffs applied to the regulator-target adjacencies or
intramodular connectivities, but by using knee point detection the
target genes are filtered in a data-driven way.

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
[`pruneModules_UIK_adj()`](https://hellmann-lab.github.io/CroCoNet/reference/pruneModules_UIK_adj.md),
[`pruneModules_topN()`](https://hellmann-lab.github.io/CroCoNet/reference/pruneModules_topN.md)

## Examples

``` r
pruned_modules_UIK_adj_kIM <- pruneModules_UIK_adj_kIM(initial_modules, consensus_network)
#> Step 1: filtering targets based on their adjacencies to the regulator
#> Median module size after filtering: 92.5
#> Step 2: filtering targets based on their intramodular connectivities
#> Median module size after filtering: 33
```
