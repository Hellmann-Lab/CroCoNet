# Prune modules based on the regulator-target adjacencies by keeping a fixed number of top targets

Prunes the initial modules by keeping a fixed number of targets per
transcriptional regulator with the highest regulator-target adjacencies.

## Usage

``` r
pruneModules_topN(initial_modules, N = 50L)
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

- N:

  Either an integer specifying a single desired pruned module size for
  all modules or a named integer vector specifying the desired pruned
  module size for each regulator (default: 50).

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

Each pruned module output by the function contains the regulator and its
`N` best target genes. When choosing the best targets, the genes are
ranked based on how strongly they are connected to the regulator
(regulator-target edge weight/adjacency).

Based on prior biological knowledge, `N` can be set to a different value
for different regulators, however, in most cases it will just be the
same desired module size for all modules. Fixing the sizes of all
modules to the same number is a simple but widespread approach.

The modules are allowed to overlap, and in addition to having its own
module, a regulator can be assigned to another regulator's module as
well, in line with the notion that genes can be multifunctional and gene
regulation can be combinatorial.

## See also

Other methods to prune modules:
[`pruneModules_UIK_adj()`](https://hellmann-lab.github.io/CroCoNet/reference/pruneModules_UIK_adj.md),
[`pruneModules_UIK_adj_kIM()`](https://hellmann-lab.github.io/CroCoNet/reference/pruneModules_UIK_adj_kIM.md)

## Examples

``` r
pruned_modules_top50 <- pruneModules_topN(initial_modules, 50)
```
