# Apply soft-thresholding

Apply soft-thresholding

## Usage

``` r
applyPower(network_list, power, n_cores = 1L)
```

## Arguments

- network_list:

  A named list of `igraph` objects containing the networks of all
  replicates.

- power:

  Numeric, the power to raise the edge weights to.

- n_cores:

  Integer, the number of cores (default: 1).

## Value

A named list of `igraph` objects containing the networks of all
replicates after the scale-free transformation on the edge weights
