# Select power for scale-free transfrmation

Select power for scale-free transfrmation

## Usage

``` r
selectPower(scale_free_fit, R2_cutoff = 0.85, rescue = TRUE)
```

## Arguments

- scale_free_fit:

  A data frame containing the \\R^2\\ values of the scale-free model fit
  for each network and soft-thresholding power.

- R2_cutoff:

  Numeric, the \\R^2\\ cutoff each network needs to pass when selecting
  the power (default: 0.85).

- rescue:

  Logical indicating whether to still select a power if none of the
  powers achieve \\R^2\\ values greater than the specified `R2_cut` for
  all networks. If TRUE (the default), the power with the highest
  minimum \\R^2\\ across the networks is output in case the cutoff is
  not met, if FALSE, the function gives an error in case the cutoff is
  not met.

## Value

Numeric, the best power to use for the scale-free transformation.
