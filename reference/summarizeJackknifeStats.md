# Summarize jackknife statistics

Calculates the estimate and confidence interval of a module-level
statistic based on all values obtained by jackknifing the module.

## Usage

``` r
summarizeJackknifeStats(
  stats_df,
  stats = setdiff(colnames(stats_df), c("regulator", "type", "id", "gene_removed",
    "module_size", "replicate", "species", "replicate1", "replicate2", "species1",
    "species2")),
  summary_method = ifelse(endsWith(stats, "monophyl"), "mean", "median"),
  conf_level = 0.95
)
```

## Arguments

- stats_df:

  Data frame of the statistics per jackknife module version:

  regulator

  :   Character, transcriptional regulator.

  module_size

  :   Integer, the numer of target genes assigned to a regulato.r

  type

  :   Character, module type (orig = original or jk = jackknifed).

  id

  :   Character, the unique ID of the module version (format:
      nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk'
      and nameOfRegulator_orig in case of module type 'orig').

  gene_removed

  :   Character, the name of the gene removed by jackknifing (NA in case
      of module type 'orig').

  replicate

  :   Character, the name of the replicate (optional, only needed if the
      statistics are calculated per replicate).

  species

  :   Character, the name of the species (optional, only needed if the
      statistics are calculated per replicate or per species).

  replicate1, replicate2

  :   Character, the names of the replicates compared (optional, only
      needed if the statistics are calculated per replicate pair).

  species1, species2

  :   Character, the names of the species compared (optional, only
      needed if the statistics are calculated per replicate pair or
      species pair).

  {{nameOfStat}}

  :   Numeric, integer or logical, one or more columns containing the
      values of the statistic(s) specified in 'stats'.

- stats:

  Character, the name of the statistic(s) that need to be summarized.

- summary_method:

  Character or character vector, the summary method ("mean" or "median")
  to be used for the statistics in 'stats'. If only one value is
  provided, the same summary method will be used for all statistics, if
  a vector is provided, each element of the vector will be used for the
  corresponding element in 'stats'. By default, "mean" for statistics
  measuring monophyleticity and "median" for all other statistics.

- conf_level:

  Numeric, confidence level of the interval (default: 0.95).

## Value

A data frame of the statistics per module:

- regulator:

  Character, transcriptional regulator.

- module_size:

  Module size, the numer of target genes assigned to a regulator.

- replicate:

  The name of the replicate (only if the column is present in the input
  'stats_df').

- species:

  The name of the species (only if the column is present in the input
  'stats_df').

- replicate1, replicate2:

  The names of the replicates compared (only if the column is present in
  the input 'stats_df').

- species1, species2:

  The names of the species compared (only if the column is present in
  the input 'stats_df').

- {{nameOfStat}}:

  Numeric, one or more columns containing the estimates (mean or median)
  of the statistic(s) specified in 'stats'.

- var\_{{nameOfStat}}:

  Numeric, one or more columns containing the variances of the
  statistic(s) specified in 'stats'.

- lwr\_{{nameOfStat}}:

  Numeric, one or more columns containing the lower bounds of the
  confidence interval for the statistic(s) specified in 'stats'.

- upr\_{{nameOfStat}}:

  Numeric, one or more columns containing the upper bounds of the
  confidence interval for the statistic(s) specified in 'stats'.

## Details

To gain information about the confidence of various statistics,
jackknifing can be used: each member gene of a module is removed and the
statistics are re-calculated. This way, the median or mean and its
confidence interval across the jackknifed versions can be used to
estimate a statistic of interest for the module as a whole. For
continuous statistics the median, for boolean statistics the mean is
recommended as the summary method.

## Examples

``` r
pres_stats <- summarizeJackknifeStats(pres_stats_jk)
tree_stats <- summarizeJackknifeStats(tree_stats_jk,
                                      c("total_tree_length", "within_species_diversity"))
```
