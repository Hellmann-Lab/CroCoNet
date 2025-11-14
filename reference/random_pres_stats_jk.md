# Preservation statistics of the the original and jackknifed random modules

Correlation of intramodular connectivities (cor_kIM) per replicate pair
for the original and all jackknifed versions of the random modules. The
jackknifed versions of the modules were created by removing each target
gene assigned to a module (the regulators were never excluded). The
preservation statistic cor.kIM was then calculated for the original
module as well as each jackknife module version by comparing each
replicate to all others, both within and across species. cor.kIM
quantifies how well the connectivity patterns are preserved between the
networks of two replicates, mathematically it is the correlation of the
intramodular connectivities per module member gene in the network of the
1st replicate VS the intramodular connectivities per module member gene
in the network of the 2nd replicate.

## Usage

``` r
random_pres_stats_jk
```

## Format

A data frame with 8352 rows and 8 columns:

- regulator:

  Character, transcriptional regulator.

- module_size:

  Integer, the numer of target genes assigned to a regulator.

- type:

  Character, module type (orig = original or jk = jackknifed).

- id:

  Character, the unique ID of the module version (format:
  nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and
  nameOfRegulator_orig in case of module type 'orig').

- gene_removed:

  Character, the name of the gene removed by jackknifing (NA in case of
  module type 'orig').

- replicate1, replicate2:

  Character the names of the replicates compared.

- species1, species2:

  Character, the names of the species `replicate1` and `replicate2`
  belong to, respectively.

- cor_adj:

  Numeric, correlation of adjacencies.

- cor_kIM:

  Numeric, correlation of intramodular connectivities.
