# Initial modules

Initial modules created by assigning the top 250 targets to each of 7
transcriptional regulators involved in the early neuronal
differentiation of primates (`regulators`). For each regulator, the top
250 targets were selected based edge weights in the consensus network:
all targets of the regulator were ranked based on their edge weight to
the regulator (regulator-taregt adjacency) and the 250 targets with the
highest regulator-target adjacencies were kept.

## Usage

``` r
initial_modules
```

## Format

A data frame with 1750 rows and 8 columns:

- regulator:

  Character, transcriptional regulator.

- target:

  Target gene of the transcriptional regulator (member of the
  regulator's initial module).

- weight:

  Consensus edge weight/adjacency, the weighted average of
  replicate-wise adjacencies.

- rho:

  Approximate Spearman's correlation coefficient of the 2 genes'
  expression profiles that form the edge.

- p.adj:

  BH-corrected approximate p-value of rho.

- direction:

  Direction of the interaction between the 2 genes that form the edge
  ("+" or "-").
