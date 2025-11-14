# SCE object of the primate neural differentiation dataset

A subset of the primate neural differentiation scRNA-seq dataset in an
SCE format. The data was collected during the early neural
differentiation of human, gorilla and cynomolgus macaque iPS cells with
3 human, 2 gorilla and 4 cynomolgus cell lines (replicates). Using a
directed differentiation protocol, cells were differentiated into neural
progenitor cells (NPCs) over the course of 9 days, and scRNA-seq data
was obtained at six time points (days 0, 1, 3, 5, 7 and 9) during this
process. The SCE object contains the raw and log-normalized counts as
well as the metadata for 300 genes and 900 cells (100 cells per
replicate).

## Usage

``` r
sce
```

## Format

An SCE object with 300 rows, 900 columns, 9 metadata columns and 2
assays.

Metadata columns:

- species:

  Name of the species.

- replicate:

  Name of the replicate/cell line.

- day:

  The day when the cell was collected.

- n_UMIs:

  Number of UMIs detected.

- n_genes:

  Number of genes detected.

- perc_mito:

  Percent of mitochondrial reads.

- sizeFactor:

  Size factor for scaling normalization, calculated first per replicate
  using \[scran::computeSumFactors\] and \[scran::quickCluster\], then
  adjusted by \[batchelor::multiBatchNorm\] to remove systematic
  differences in covergae across replicates.

- pseudotime:

  Pseudotime inferred by \[SCORPIUS::infer_trajectory\].

- cell_type:

  Cell type labels predicted by \[SingleR::classifySingleR\] using the
  embryoid body dataset from Rhodes et al. 2022 as reference.

Assays:

- counts:

  Raw counts.

- logcounts:

  Log-normalized counts created by first calculating size factors per
  replicate using \[scran::computeSumFactors\] and
  \[scran::quickCluster\], then adjusted them to remove systematic
  differences in covergae across replicates and log-normalizing by
  \[batchelor::multiBatchNorm\].
