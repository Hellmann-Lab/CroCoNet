# Distance measures of the the original and jackknifed random modules

Distance measures per replicate pair for the original and all jackknifed
versions of the random modules.

## Usage

``` r
random_dist_jk
```

## Format

A named list with 232 elements containing the distance measures per
(original or jackknifed) module version. Each element is a data frame
with 72 rows and 10 columns:

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

- dist:

  Numeric, distance measure ranging from 0 to 1, calculated based on the
  correlation of intramodular connectivities.

## Details

The jackknifed versions of the modules were created by removing each
target gene assigned to a module (the regulators were never excluded).
The distance measures were calculated based on the correlation of
intramodular connectivities: \$\$dist = \frac{1 - cor.kIM}{2}\$\$ (a
correlation of 1 corresponds to a distance of 0, whereas a correlation
of -1 corresponds to a distance of 1). Each element of the list
corresponds a jackknife module version and contains the distance
measures between all possible pairs of replicates for this module
version.
