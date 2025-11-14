# Calculate preservation statistics

Calculates one or both of two preservation statistics - correlation of
adjacencies (cor_adj) and correlation of intramodular connectivities
(cor_kIM) - for all modules and replicate pairs with or without
jackknifing.

## Usage

``` r
calculatePresStats(
  pruned_modules,
  network_list,
  stats = c("cor_adj", "cor_kIM"),
  replicate2species = NULL,
  jackknife = TRUE,
  signed = FALSE,
  n_cores = 1L,
  corr_method = "spearman"
)
```

## Arguments

- pruned_modules:

  Data frame of pruned modules, required columns:

  regulator

  :   Character, transcriptional regulator.

  target

  :   Character, target gene of the transcriptional regulator (member of
      the regulator's pruned module).

- network_list:

  A named list of `igraph` objects containing the networks of all
  replicates.

- stats:

  Character or character vector specifying which preservation statistics
  to calculate (one or more of "cor_adj", "cor_kIM", default:
  c("cor_adj", "cor_kIM")).

- replicate2species:

  A data frame specifying which species each replicate belongs to,
  required columns:

  replicate

  :   Character, name of the replicate.

  species

  :   Character, name of the species.

  If NULL (default), the output will contain no species information.

- jackknife:

  Logical specifying whether jackknifing should be performed or not
  (default: TRUE).

- signed:

  Logical indicating whether the networks in `network_list` are signed
  (default: FALSE, see also
  [`normalizeEdgeWeights`](https://hellmann-lab.github.io/CroCoNet/reference/normalizeEdgeWeights.md)).
  If set to FALSE and `network_list` contains the edge attribute
  "direction", the edge weights in the network with direction "-" are
  negated for the calculation of "cor_adj".

- n_cores:

  Integer, the number of cores (default: 1).

- corr_method:

  Character, the method for the calculation of correlation, one of
  "spearman", "pearson", "kendall" (default: "spearman").

## Value

Data frame of the preservation statistics with the following columns:

- regulator:

  Character, transcriptional regulator.

- module_size:

  Integer, the number of target genes assigned to a regulator (only
  present if the column is also present in `pruned_modules`).

- type:

  Character, module type ("orig" = original or "jk" = jackknifed, only
  present if parameter `jackknife` is set to TRUE).

- id:

  Character, the unique ID of the module version (format:
  nameOfRegulator_jk_nameOfGeneRemoved in case of module type "jk" and
  nameOfRegulator_orig in case of module type "orig", only present if
  parameter `jackknife` is set to TRUE).

- gene_removed:

  Character, the name of the gene removed by jackknifing (NA in case of
  module type "orig", only present if parameter `jackknife` is set to
  TRUE).

- replicate1, replicate2:

  Character, the names of the replicates compared.

- species1, species2:

  Character, the names of the species `replicate1` and `replicate2`
  belong to, respectively (only present if `replicate2species` is not
  NULL).

- {{nameOfStat}}:

  Numeric, one or more columns containing the preservation statistics
  specified by the parameter `stats`.

## Details

The function calculates two preservation statistics adapted from WGCNA
(cor_adj and cor_kIM) for each module and each replicate pair. Both
statistics quantify how well the topology of a module is preserved
between the networks of two replicates. The statistic cor_adj is the
correlation of all edge weights within the module in the network of
replicate1 VS in the network of replicate2, while the statistic cor_kIM
is the correlation of the intramodular connectivities per module member
gene in the network of replicate1 VS in the network of replicate2.

All statistics assume a joint module assignment (typically derived from
the consensus network) but compare topological properties directly
between the replicate-wise networks. In this approach, a module is
always defined as the same set of genes, but the adjacencies/edge
weights among these genes could differ from replicate to replicate;
poorly preserved modules are expected to have many, while well-preserved
modules are expected to have few such differences.

If `jackknife` is set to TRUE (the default), the function creates all
possible jackknifed versions of each input module by removing each
target gene assigned to that module (the regulator is never excluded),
then it calculates the preservation statistics for all of these
jackknifed module versions in addition to the original module. This way,
a confidence interval can be calculated for each module and statistic
(see
[`summarizeJackknifeStats`](https://hellmann-lab.github.io/CroCoNet/reference/summarizeJackknifeStats.md)).
Later on in the pipeline, jackknifing can also provide information about
which target genes within a conserved/diverged module are particularly
responsible for the conservation/divergence (see
[`findConservedDivergedTargets`](https://hellmann-lab.github.io/CroCoNet/reference/findConservedDivergedTargets.md)).
If jackknifing is not desired, please set `jackknife` to FALSE. This
will also substantially reduce running times.

## References

Langfelder, P., Luo, R., Oldham, M. C., & Horvath, S. (2011). Is my
network module preserved and reproducible? PLoS Computational Biology,
7(1), 1001057.

## Examples

``` r
pres_stats_jk <- calculatePresStats(pruned_modules, network_list, "cor_kIM", replicate2species)
```
