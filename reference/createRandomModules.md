# Create random modules

Creates a matching random module for each actual module.

## Usage

``` r
createRandomModules(pruned_modules, network_genes, seed = 42)
```

## Arguments

- pruned_modules:

  Data frame of pruned modules, required columns:

  regulator

  :   Character, transcriptional regulator.

  target

  :   Character, target gene of the transcriptional regulator (member of
      the regulator's pruned module).

- network_genes:

  Character vector of all genes in the network.

- seed:

  Integer, the seed to use for the random sampling (default: 42).

## Value

Data frame of the random modules with the following columns:

- regulator:

  Character, transcriptional regulator.

- module_size:

  Integer, the number of genes assigned to a regulator (only present if
  the column is also present in the input `pruned_modules`).

- target:

  Character, member gene of the regulator's random module.

## Details

The function outputs a random module for each module in
`pruned_modules`. The random modules have the same regulators and
contain the same number of target genes as the original modules, but
these target genes are randomly drawn from `network_genes`.

In the next steps of the pipeline, the actual modules are compared to
these random modules in terms of the statistics calculated to check
whether the 2 groups of modules behave in general differently (see
[`plotPresStatDistributions`](https://hellmann-lab.github.io/CroCoNet/reference/plotPresStatDistributions.md),
[`plotPresStats`](https://hellmann-lab.github.io/CroCoNet/reference/plotPresStats.md),
[`plotTreeStatDistributions`](https://hellmann-lab.github.io/CroCoNet/reference/plotTreeStatDistributions.md)
and
[`plotTreeStats`](https://hellmann-lab.github.io/CroCoNet/reference/plotTreeStats.md))
and to remove those individual actual modules that show too similar
characteristics to the random modules (see
[`filterModuleTrees`](https://hellmann-lab.github.io/CroCoNet/reference/filterModuleTrees.md)).

## Examples

``` r
random_modules <- createRandomModules(pruned_modules, genes)
```
