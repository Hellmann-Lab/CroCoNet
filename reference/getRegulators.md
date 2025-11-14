# Get transcriptional regulators

Selects relevant transcriptonal regulators by intersecting highly
variable genes with known transcriptional regulators from the specified
source(s). These regulators can be used in subsequent steps of the
CroCoNet pipeline to assemble modules around them.

## Usage

``` r
getRegulators(sce, source = "jaspar_core", ...)
```

## Arguments

- sce:

  `SingleCellExperiment` object containing the expression data
  (logcounts and metadata) for all network genes. Required metadata
  columns:

  species

  :   Character, the name of the species.

  replicate

  :   Character, the name of the replicate/cell line.

- source:

  Character or character vector specifying the source of known
  transcriptional regulators. It can be one or more of "jaspar_core"
  (JASPAR 2024 vertebrate core, the default), "jaspar_unvalidated"
  (JASPAR 2024 unvalidated), and "image" (IMAGE database from Madsen et
  al. 2018), or alternatively, it can be a custom list of regulators
  provided by the user.

- ...:

  Additional arguments passed to `getTopHVGs`.

## Value

A character vector of known transcriptional regulators that are highly
variable in at least one species of the data.

## Details

This approach combines information from the dataset at hand with prior
biological knowledge to pinpoint transcriptional regulators that could
be interesting to investigate within the regulatory networks.

As for the information from the dataset, the function identifies highly
variable genes in the data using `modelGeneVar` and `getTopHVGs`.
Briefly, for each species a trend is fitted between the variance and
mean of the log-expression values across all genes, and the fitted value
of a gene is regarded as the technical component of variation, while the
residual from the trend is regarded as the biological component of
variation. All genes with a positive biological component are selected
as highly variable genes in each species. If more/less lenient criteria
are desired, further arguments can be passed to `getTopHVGs` (e.g. `n`,
`prop`, `var.threshold`, or `fdr.threshold`). As the final set of highly
variable genes, the union of highly variable genes is taken across
species.

As for the prior biological knowledge, the function extracts a list of
known transcriptional regulators that have at least 1 annotated motif in
the motif database(s) specified in `source`. This database can be JASPAR
2024 vertebrate core ("jaspar_core"), JASPAR 2024 unvalidated
("jaspar_unvalidated"), the IMAGE database from Madsen et al. 2018
("image"), or any combination of the above. Alternatively, the user can
provide a custom list of regulators selected by their preferred method.

Finally, the highly variable genes and known transcriptional regulators
are intersected to provide the final list of regulators that are
relevant for the biological process at hand. In the next steps of the
pipeline, co-expression modules are assigned to each of these regulators
(see
[`assignInitialModules`](https://hellmann-lab.github.io/CroCoNet/reference/assignInitialModules.md)
and
[`pruneModules`](https://hellmann-lab.github.io/CroCoNet/reference/pruneModules.md)),
which can then be studied in terms of their cross-species conservation.

## References

Rauluseviciute, I., Riudavets-Puig, R., Blanc-Mathieu, R.,
Castro-Mondragon, J. A., Ferenc, K., Kumar, V., Lemma, R. B., Lucas, J.,
Chèneby, J., Baranasic, D., Khan, A., Fornes, O., Gundersen, S.,
Johansen, M., Hovig, E., Lenhard, B., Sandelin, A., Wasserman, W. W.,
Parcy, F., & Mathelier, A. (2023). JASPAR 2024: 20th anniversary of the
open-access database of transcription factor binding profiles. Nucleic
Acids Research, 52(D1), D174–D182.

Madsen, J. G. S., Rauch, A., Van Hauwaert, E. L., Schmidt, S. F.,
Winnefeld, M., & Mandrup, S. (2018). Integrated analysis of motif
activity and gene expression changes of transcription factors. Genome
Research, 28(2), 243–255.

## Examples

``` r
regulators <- getRegulators(sce)
regulators2 <- getRegulators(sce, n = 2000)
regulators3 <- getRegulators(sce, c("jaspar_core", "jaspar_unvalidated", "image"))
```
