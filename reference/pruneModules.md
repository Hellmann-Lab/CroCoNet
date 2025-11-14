# Prune modules

Prunes the initial modules by keeping only the best targets of each
transcriptional regulator.

## Usage

``` r
pruneModules(
  initial_modules,
  method = c("UIK_adj", "UIK_adj_kIM", "topN"),
  consensus_network = NULL,
  min_module_size = 20L,
  max_frac_modules_lost = 0,
  exponent = 1L,
  N = 50L
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

- method:

  Character, the pruning method, one of "UIK_adj", "UIK_adj_kIM",
  "topN".

- consensus_network:

  `igraph` object, the consensus network across all species and
  replicates.

- min_module_size:

  Integer, the lower threshold of module size in case of the methods
  "UIK_adj" and "UIK_adj_kIM" (default: 20).

- max_frac_modules_lost:

  Numeric, the threshold for the fraction of removed modules in case of
  the methods "UIK_adj" and "UIK_adj_kIM" (default: 0).

- exponent:

  Integer, the exponent the regulator-target adjacency and intramodular
  connectivity is raised to the power of during the cumulative sum curve
  calculation in case of the methods "UIK_adj" and "UIK_adj_kIM"
  (default: 1, i.e. the regulator-target adjacencies and intramodular
  connectivities stay unchanged).

- N:

  An integer or a named integer vector specifying the desired pruned
  module size(s) in case of the method "topN" (default: 50).

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

3 methods are implemented to choose the best targets:

- topN: Takes a fixed number of targets per regulator with the highest
  regulator-target adjacencies (for details see
  [`pruneModules_topN`](https://hellmann-lab.github.io/CroCoNet/reference/pruneModules_topN.md)).

- UIK_adj: Applies a dynamic stepwise pruning based on the
  regulator-target adjacencies (for details see
  [`pruneModules_UIK_adj`](https://hellmann-lab.github.io/CroCoNet/reference/pruneModules_UIK_adj.md)).

- UIK_adj_kIM: Applies a dynamic stepwise pruning based on the
  regulator-target adjacencies and intramodular connectivities (for
  details see
  [`pruneModules_UIK_adj_kIM`](https://hellmann-lab.github.io/CroCoNet/reference/pruneModules_UIK_adj_kIM.md)).

## Examples

``` r
pruned_modules <- pruneModules(initial_modules, "topN", N = 30)
pruned_modules <- pruneModules(initial_modules, "UIK_adj")
#> Step 1: filtering targets based on their adjacencies to the regulator
#> Median module size after filtering: 92.5
#> Step 2: filtering targets based on their adjacencies to the regulator
#> Median module size after filtering: 35
pruned_modules <- pruneModules(initial_modules, "UIK_adj_kIM", consensus_network)
#> Step 1: filtering targets based on their adjacencies to the regulator
#> Median module size after filtering: 92.5
#> Step 2: filtering targets based on their intramodular connectivities
#> Median module size after filtering: 33
```
