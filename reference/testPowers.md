# Test soft-thresholding power values

Applies a range of soft-thresholding powers to the networks provided and
determines the goodness of fit for each soft power and network.

## Usage

``` r
testPowers(network_list, powers = 1:20, n_cores = 1L)
```

## Arguments

- network_list:

  A named list that contains the networks of each replicate and the
  consensus network in an `igraph` format.

- powers:

  A numeric or integer vector specifying the soft-thresholding powers to
  test.

- n_cores:

  Integer, the number of cores (default: 1).

## Value

A data frame containing the \\R^2\\ values of the scale-free model fit
for each network and soft-thresholding power.

## Examples

``` r
scale_free_fit <- testPowers(network_list, c(1, seq(2, 10, by = 2)))
```
