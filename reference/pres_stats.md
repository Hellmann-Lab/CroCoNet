# Preservation statistics of the pruned modules

Correlation of intramodular connectivities (cor_kIM) per replicate pair
and module. The preservation statistic cor.kIM quantifies how well the
connectivity patterns are preserved between the networks of two
replicates, mathematically it is the correlation of the intramodular
connectivities per module member gene in the network of the 1st
replicate VS the intramodular connectivities per module member gene in
the network of the 2nd replicate. This statistic was calculated for all
possible jackknifed versions of the modules, each of which was created
by removing a target gene assigned to the given module (the regulators
were never excluded). Each jackknifed module was compared between all
posible pairs of replicates, both within and across species, resulting
in a cor.kIM value per jackknifed module version and replicate pair.
Finally, the cor.kIM values were summarized per module and replicate
pair by taking the median and its 95% confidence interval across all
jackknifed module versions.

## Usage

``` r
pres_stats
```

## Format

A data frame with 252 rows and 10 columns:

- regulator:

  Character, transcriptional regulator.

- module_size:

  Module size, the numer of target genes assigned to a regulator.

- replicate1, replicate2:

  The names of the replicates compared.

- species1, species2:

  The names of the species 'replicate1' and 'replicate2' belongs to,
  respectively.

- cor_kIM:

  The median of cor.kIM across all jackknifed versions of the module.

- var_cor_kIM:

  The variance of cor.kIM across all jackknifed versions of the module.

- lwr_cor_kIM:

  The lower bound of the 95% confidence interval of cor.kIM calculated
  by jackknifing.

- upr_cor_kIM:

  The upper bound of the 95% confidence interval of cor.kIM calculated
  by jackknifing.

- cor_adj:

  The median of cor.adj across all jackknifed versions of the module.

- var_cor_adj:

  The variance of cor.adj across all jackknifed versions of the module.

- lwr_cor_adj:

  The lower bound of the 95% confidence interval of cor.adj calculated
  by jackknifing.

- upr_cor_adj:

  The upper bound of the 95% confidence interval of cor.adj calculated
  by jackknifing.
