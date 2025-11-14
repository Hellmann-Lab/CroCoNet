# Linear model for the characterization of overall module conservation

A weighted linear model between the total tree length and within-species
diversity of the 12 modules in the subsetted early neuronal
differentiation dataset. It captures the general relationship between
the two tree-based statistics in this dataset and thus also describes
the expected total tree length for any given within-species diversity.
By identifying modules that differ the most from this expectation, i.e.
identifying the outlier data points, it can be used to pinpoint
conserved and overall diverged modules. In combination with jackknifing,
it can also be used to identify target genes within these modules that
contribute the most to the conservation/divergence. The linear model was
fit by calling [`lm`](https://rdrr.io/r/stats/lm.html) with weights
inversely proportional to the variance of the total tree lengths.

## Usage

``` r
lm_overall
```

## Format

An object of class [`lm`](https://rdrr.io/r/stats/lm.html) fitted on the
object `tree_stats` with the formula "total_tree_length ~
within_species_diversity".
