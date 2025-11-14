# Distance measures of the pruned modules

Distance measures per replicate pair for the original (i.e. not
jackknifed) pruned modules.

## Usage

``` r
dist
```

## Format

A named list with 12 elements containing the distance measures per
module. Each element is a data frame with 72 rows and 7 columns:

- regulator:

  Character, transcriptional regulator.

- module_size:

  Integer, the numer of target genes assigned to a regulator.

- replicate1, replicate2:

  Character the names of the replicates compared.

- species1, species2:

  Character, the names of the species `replicate1` and `replicate2`
  belong to, respectively.

- dist:

  Numeric, distance measure ranging from 0 to 1, calculated based on the
  correlation of intramodular connectivities.

## Details

The distance measures were calculated based on the correlation of
intramodular connectivities: \$\$dist = \frac{1 - cor.kIM}{2}\$\$ (a
correlation of 1 corresponds to a distance of 0, whereas a correlation
of -1 corresponds to a distance of 1). Each element of the list
corresponds to a module and contains the distance measures between all
possible pairs of replicates for this module.
