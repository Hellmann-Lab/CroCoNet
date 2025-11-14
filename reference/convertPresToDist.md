# Convert preservation statistics to distances across all modules

Converts preservation statistics between replicates to distance measures
ranging from 0 to 1.

## Usage

``` r
convertPresToDist(pres_stats, stat, min_stat = -1, max_stat = 1, n_cores = 1L)
```

## Arguments

- pres_stats:

  Data frame of the preservation statistics, required columns:

  regulator

  :   Character, transcriptional regulator.

  module_size

  :   Integer, the numer of target genes assigned to a regulator.

  type

  :   Character, module type ("orig" = original or "jk" = jackknifed,
      only needed if the preservation statistics were calculated with
      jackknifing).

  id

  :   Character, the unique ID of the module version (format:
      nameOfRegulator_jk_nameOfGeneRemoved in case of module type "jk"
      and nameOfRegulator_orig in case of module type "orig", only
      needed if the preservation statistics were calculated with
      jackknifing).

  gene_removed

  :   Character, the name of the gene removed by jackknifing (NA in case
      of module type "orig", only needed if the preservation statistics
      were calculated with jackknifing).

  replicate1, replicate2

  :   Character the names of the replicates compared.

  species1, species2

  :   Character, the names of the species `replicate1` and `replicate2`
      belong to, respectively.

  {{stat}}

  :   Numeric, the preservation statistic specified by the parameter
      `stat`.

- stat:

  Character, the name of the column containing the preservation
  statistic that is to be converted into a distance measure.

- min_stat, max_stat:

  Numeric, the theoretical minimum and maximum value of `stat` (default:
  -1 and 1, respectively). For the preservation statistics implemented
  by CroConet and all other correlation-based statistics, please leave
  `min_stat` and `max_stat` at -1 and 1, respectively. For custom
  preservation statistics, a different value might have to be used. If
  set to NULL or ±Inf, the conversion to a distance measure is performed
  using the empirical minimum and maximum.

- n_cores:

  Integer, the number of cores (default: 1).

## Value

A named list containing the distance measures as data frames for all
modules or jackknifed module versions in the input `pres_stats`. The
data frames contain 1 new column in addition to the relevant columns of
`pres_stats`:

- dist:

  Numeric, the distance measure ranging from 0 to 1 calculated based on
  `stat`.

While `pres_stats` contains only non-redundant replicate pairs, the data
frames in the output contain all possible replicate pairs (i.e.
replicateA-replicateB and replicateB-replicateA are 2 separate entries
with the same distance).

## Details

As part of the CroCoNet approach, pairwise module preservation scores
are calculated between replicates, both within and across species (see
[`calculatePresStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculatePresStats.md))
to gain information about the cross-species differences but also about
the within-species diversity of the modules. These correlation-based
preservation statistics quantify how well the module topology is
preserved between the networks of two replicates.

This function converts a chosen preservation statistic (*p*) specified
by the argument `stat` into a distance measure (*d*) using the following
formula: \$\$d = \frac{max(p) - p}{max(p) - min(p)}\$\$ If the
theoretical minimum and maximum of the preservation statistic are known,
these can be provided using the parameters `min_stat` and `max_stat`. As
the preservation statistics implemented by CroCoNet are
correlation-based, they all range between -1 and 1, and thus `min_stat`
and `max_stat` should be set to -1 and 1, respectively (default). If a
custom preservation statistic is used, `min_stat` and `max_stat` might
have to be set to different values. If they are set to NULL or ±Inf, the
minimum and maximum of the preservation statistic are calculated
empirically using the data.

The function also splits up the distance measures into modules/jackknife
module versions (depending on whether the preservation statistics were
calculated with or without jackknifing, see see the parameter
`jackknife` in the function
[`calculatePresStats`](https://hellmann-lab.github.io/CroCoNet/reference/calculatePresStats.md))
and outputs a list of data frames per module/jackknife module version.
Modules/jackknifed module versions where the preservation statistic for
any of the replicate pairs is NA are removed.

In the next step of the pipeline, these distance measures are used to
reconstruct a neighbor-joining tree per module/jackknife module version
that represent the dissimilarity of module topology across all
replicates (see
[`reconstructTrees`](https://hellmann-lab.github.io/CroCoNet/reference/reconstructTrees.md)).

## Examples

``` r
dist_jk <- convertPresToDist(pres_stats_jk, "cor_kIM")
```
