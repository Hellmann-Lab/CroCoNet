# Normalize edge weights between 0 and 1

Transforms the edge weights of all networks to be in the range \[0, 1\].

## Usage

``` r
normalizeEdgeWeights(
  network_list,
  signed = FALSE,
  min_weight = NULL,
  max_weight = NULL,
  n_cores = 1L
)
```

## Arguments

- network_list:

  A named list of `igraph` objects containing the networks of all
  replicates.

- signed:

  Logical indicating whether a signed network is desired (default:
  FALSE, see `Details`).

- min_weight, max_weight:

  Numeric, the theoretical minimum and maximum values of the edge
  weights. If set to NULL (default), the normalization is performed
  using the empirical minimum and maximum. For correlation-based edge
  weights, please set `min_weight` and `max_weight` to -1 and 1,
  respectively.

- n_cores:

  Integer, the number of cores (default: 1).

## Value

A named list of `igraph` objects containing the networks of all
replicates, with the edge weights normalized between 0 and 1. If
originally there were both positive and negative edge weights present in
the data, a new edge attribute is added to all `igraph` objects:

- direction:

  Character, the direction of the interaction between the 2 genes that
  form the edge ("+" or "-" depending on the sign of the edge weight
  before normalization).

## Details

Normalizing the edge weights between 0 and 1 makes them interpretable as
adjacencies and ensures that network concepts such as connectivity are
applicable.

There are 2 approaches for the normalization:

- If `signed` is set to FALSE (unsigned network, the default), gene
  pairs with high negative edge weights are considered as connected as
  gene pairs with high positive edge weights. Therefore the negative
  edge weights are first replaced by their absolute values, then all
  edge weights are scaled by the maximum weight across all networks:
  \$\$w\_{new} = \frac{\|w\|}{max(\|w\|)}\$\$ After the transformation,
  the edge weights/adjacencies around 0 correspond to the former low
  positive and low negative values, while the edge weights/adjacencies
  around 1 correspond to the former high positive and high negative
  values.

- If `signed` is set to TRUE (signed network), gene pairs with high
  negative edge weights are considered unconnected. Therefore all edge
  weights are transformed between 0 and 1 using a min-max normalization:
  \$\$w\_{new} = \frac{w - min(w)}{max(w) - min(w)}\$\$ After the
  transformation, the edge weights/adjacencies around 0 correspond to
  the former high negative values and the edge weights/adjacencies
  around 1 correspond to the former high positive values.

If the theoretical minimum and maximum edge weights are known, these can
be provided using the parameters `min_weight` and `max_weight`. For
example, if the networks were inferred by calculating correlations
between gene expression profiles, `min_weight` should be set to -1 and
`max_weight` should be set to 1. If `min_weight` and `max_weight` are
left at NULL, the minimum and maximum edge weights are calculated
empirically using the data.

After normalization by either of the approaches above, it is not
possible to tell anymore which edges had a positive and which edges had
a negative weight originally. Since this can be a useful piece of
information (it can specify the mode of regulation: activation or
repression), the information about the sign is stored as a new edge
attribute "direction" in the output `igraph` objects ("+" if the
original edge weight was positive and "-" is the original edge weight
was negative). If all original edge weights were positive, the edge
attribute "direction" is not added.

## See also

["Signed or unsigned: which network type is preferable?" by Peter
Langfelder](https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/)

## Examples

``` r
network_list_norm <- normalizeEdgeWeights(network_list_raw)
```
