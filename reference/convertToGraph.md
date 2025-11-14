# Convert a list of data frames to a list of igraphs

Takes care of the conversion between data table and igraph formats of
the networks.

## Usage

``` r
convertToGraph(dt_list, network_genes, n_cores = 1L)
```

## Arguments

- dt_list:

  A list of data tables that contain the columns 'from', 'to' and
  'weight'. Each row should correspond to an edge in the network with
  'from' and 'to' as the end nodes and 'weight' as the edge weight.
  Additional columns will be converted to edge attributes in the
  igraphs.

- network_genes:

  Character vector of all genes in the network.

- n_cores:

  Number of cores.

## Value

A list of igraphs.
