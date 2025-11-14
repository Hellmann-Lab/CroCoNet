# Convert a tree to a data frame

Converts a tree to a data frame where each row corresponds to a branch
in the tree.

## Usage

``` r
getTreeDf(tree)
```

## Arguments

- tree:

  Object of class \[phylo\] with replicates as tips. It is expected to
  have a component 'species' that specifies which species each tip
  belongs to.

## Value

A data frame of tree branches.
