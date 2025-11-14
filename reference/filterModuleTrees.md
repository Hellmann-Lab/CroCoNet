# Filter modules based on their tree representations

Filters modules using two tree-based statistics (the total tree length
and within-species diversity) by comparing the values of the actual
modules to the values of the random modules.

## Usage

``` r
filterModuleTrees(tree_stats, random_tree_stats, p_cutoff = 0.95)
```

## Arguments

- tree_stats:

  Data frame of the tree-based statistics for the actual (pruned)
  modules. Required columns:

  regulator

  :   Character, transcriptional regulator.

  total_tree_length

  :   Numeric, total tree length per module (typically the median across
      all jackknife versions of a module.

  within_species_diversity

  :   Numeric, within-species diveristy per module (typically the median
      across all jackknife versions of a module.

- random_tree_stats:

  Data frame of the tree-based statistics for the random modules.
  Required columns:

  regulator

  :   Character, transcriptional regulator.

  total_tree_length

  :   Numeric, total tree length per module (typically the median across
      all jackknife versions of a module).

  within_species_diversity

  :   Numeric, within-species diveristy per module (typically the median
      across all jackknife versions of a module).

- p_cutoff:

  Numeric, the cutoff for the probability of a module's statistics being
  drawn from the distribution of the actual modules and not the
  distribution of random modules. Modules are retained only if this
  probability exceeds the cutoff.

## Value

Data frame of the tree-based statistics for the actual (pruned) modules
after filtering.

## Details

The within-species diversity measures the variability of module
connectivity patterns within species and the total tree length measures
the variability of module connectivity patterns both within and across
species. The within-species diversity is a subset of the total tree,
hence we expect a linear relationship between the two statistics. This
line explains the detection robustness: modules that have both low
within-species diversity and low total tree length are well-preserved
both within and across species, meaning that these modules could be
robustly detected in all replicates, whereas modules, for which both
metrics are high, are poorly preserved not just across but also within
species, indicating a high detection uncertainty ("wobbliness"). Random
modules are expected to fall at the wobbly end of the spectrum, while
actual modules are expected to fall at the robust end of the spectrum.

The function removes the modules that are uncertain to detect, i.e.
modules the are too similar to the random modules in terms of their
total tree length and within-species diversity. It calculates how
probable it is that the statistics of an actual module come from the
distribution of all actual modules (\\p\_{actual}\\) and how probable it
is that they come from the distribution of all random modules
(\\p\_{random}\\) using the probability density functions of the two
bivariate normal distributions. If \$\$\frac{p\_{actual}}{p\_{actual} +
p\_{random}} \> p\_{cutoff}\$\$ is fulfilled, the module is kept,
otherwise it is removed.
