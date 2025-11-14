# Plot the distribution of module sizes

Plots the distribution of module sizes (i.e. the number of target genes
assigned to each regulator) as a histogram.

## Usage

``` r
plotModuleSizeDistribution(modules, font_size = 14)
```

## Arguments

- modules:

  Data frame of the pruned modules, required columns:

  regulator

  :   Character, transcriptional regulator.

  module_size

  :   Module size, the numer of target genes assigned to a regulator.

- font_size:

  Numeric, font size (default: 14).

## Value

A histogram as a `ggplot` object showing the distribution of module
sizes.

## Examples

``` r
plotModuleSizeDistribution(pruned_modules)
```
