# Load networks as a list of igraph objects

Loads the TSV files containing the gene-gene edge weights and summarizes
them as a list of `igraph` objects per replicate (i.e. cell line or
biological replicate within a species).

## Usage

``` r
loadNetworks(
  path,
  replicate_names = NULL,
  rep = 1L,
  directed = TRUE,
  min_occurrence = 2L,
  n_cores = 1L
)
```

## Arguments

- path:

  Character, the path where the files containing the gene-gene edge
  weights can be found.

- replicate_names:

  Character vector, the names of the replicates that were used for
  network reconstruction. These are expected to be the base names of the
  TSV files (format: nameOfReplicate.tsv or nameOfReplicate_index.tsv).
  If set to NULL (default), the names are deduced from the file names by
  stripping anything starting with "\_" and ".tsv".

- rep:

  Integer, the number of output files per replicate, e.g. the number of
  different subsamplings or the number of independent runs in case of a
  stochastic network inference algorithm (default: 1). If `rep` \> 1,
  the TSV files are expected to be indexed from 1 to `rep` for each
  replicate (format: nameOfReplicate_index.tsv).

- directed:

  Logical indicating whether the network inference method produces a
  directed output, i.e. whether geneA-geneB and geneB-geneA can both be
  present among the edges (default: TRUE). For GRNBoost2 this has to be
  left as TRUE, while for correlation-based methods this has to most
  likely be set to FALSE.

- min_occurrence:

  Integer, the minimum number of occurrences an edges has to have across
  runs/subsamplings to be kept (default: 2). Disregarded if `rep` = 1.

- n_cores:

  Integer, the number of cores (default: 1).

## Value

A named list of `igraph` objects containing the networks per replicate.
The list names are taken from `replicate_names`. Each `igraph` object
contains the edge attribute "weight" (taken from the 3rd column of the
input TSV files). Starting from the 4th column of the input TSV files,
all columns are preserved as edge attributes with unchanged names.

## Details

The functions reads and formats the TSV files in the directory specified
by `path`. The base names of the TSV files of interest can be specified
using the parameter `replicate_names` (in this case the files are
expected to be named using the format: nameOfReplicate.tsv or
nameOfReplicate_index.tsv). If `replicate_names` is set to NULL
(default), all TSV files in `path` are read and the replicate names are
deduced from the file names by stripping anything starting with "\_" and
".tsv".

As the first 3 columns, each input TSV file is expected to contain the
first gene that forms the edge, the second gene that forms the edge and
the edge weight. There can be additional columns as well, these will be
preserved as edge attributes in the `igraph` objects.

If the edge weights were calculated on different subsamplings of cells
per replicate, or a network inference algorithm involving stochastic
steps (e.g. GRNBoost2) was run several times on each replicate, the edge
weights are averaged across these subsamplings/runs, and a single
combined `igraph` object is returned per replicate. The number of
subsamplings/runs has to be specified by the parameter `rep`. The TSV
files corresponding to the same replicate but different
subsamplings/runs are expected to have the same base name with different
indices: nameOfReplicate_index.tsv.

If the network inference method produces an output with directed edges,
i.e. geneA-geneB and geneB-geneA can both be present, the parameter
`directed` should be set to TRUE (this is the case e.g. for GRNBoost2).
In this case the edge weights inferred between the same gene pair but in
opposite directions are averaged. If the edge in one of the directions
is missing, it is regarded as a 0. For correlation-based methods (e.g.
`correlatePairs`), `directed` has to most likely be set to FALSE. In
both cases, `loadNetworks` returns an undirected network.

Rarely occurring edges can be removed altogether by specifying
`min_occurrence`. This is not relevant if the network inference was done
only once per replicate (`rep` = 1), therefore in this case the value of
`min_occurrence` is ignored. If the network inference was done several
times for each replicate (`rep` \> 1), the highest possible number of
occurrences for each edge is 2×`rep` in case of a directed network
inference method and `rep` in case of an undirected network inference
method. If an edge occurs less often than the specified value of
`min_occurrence`, the edge is removed. This can be helpful to 1) denoise
the networks and 2) decrease the computational power needed for the next
steps.

## References

Moerman, T., Santos, S. A., González-Blas, C. B., Simm, J., Moreau, Y.,
Aerts, J., & Aerts, S. (2019). GRNBoost2 and Arboreto: efficient and
scalable inference of gene regulatory networks. Bioinformatics , 35(12),
2159–2161.

Lun, A. T. L., McCarthy, D. J., & Marioni, J. C. (2016). A step-by-step
workflow for low-level analysis of single-cell RNA-seq data with
Bioconductor. F1000Research, 5, 2122.

## Examples

``` r
network_list_raw <- loadNetworks(path = system.file("extdata", package = "CroCoNet"),
                                 rep = 10)
```
