# Convert a list of igraphs to a list of data tables

Takes care of the conversion between the igraph and data table formats
of the networks.

## Usage

``` r
convertToDT(network_list)
```

## Arguments

- network_list:

  A list of igraph objects

## Value

A list of data tables with the columns 'from', 'to', 'weight' and any
additional edge attributes the input igraphs contain.
