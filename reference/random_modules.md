# Random modules

Random modules that match the actual (pruned) modules in size. These
random modules have the same regulators and contain the same number of
target genes as the actual modules, but the target genes were randomly
drawn from all genes in the network.

## Usage

``` r
random_modules
```

## Format

A data frame with 225 rows and 3 columns:

- regulator:

  Character, transcriptional regulator.

- target:

  Member of the regulator's random module.

- module_size:

  Module size, the numer of target genes assigned to a regulator.
