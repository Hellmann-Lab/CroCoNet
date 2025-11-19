# CroCoNet

  
`CroCoNet` (**Cro**ss-species **Co**mparison of **Net**works) is a
computational pipeline and R package to quantitatively compare
gene-regulatory networks across species. The approach hinges on
contrasting network variability within and across species in order to
distinguish divergence from detection uncertainty. By enabling robust
network analyses, CroCoNet contributes towards finding meaningful
cross-species differences in gene regulation.

In this vignette, we demonstrate the complete CroCoNet workflow using a
single-cell RNA-seq dataset of early primate neural differentiation as
an example. The data files for this analysis are available in the
project’s [Zenodo archive](https://zenodo.org/records/17610308), in the
compressed folder “neural_differentiation_dataset.zip”. All file paths
used in the examples below are relative to this master folder.

Through this tutorial, users will learn how to:

- Prepare and format input data for CroCoNet (including preprocessing
  and network inference)

- Process networks and build the consensus network

- Assign and characterize modules

- Calculate module preservation scores within and across species

- Quantify module conservation across all species and along specific
  lineages

- Visualize modules using network- and tree-based representations  
    

## Input data

### Experimental design

`CroCoNet` was designed to analyze scRNA-seq data from multiple species,
with multiple replicates per species. For a meaningful analysis, the
conditions, pseudotemporal trajectory and cell type composition should
be as comparable across species as possible. To reliably infer networks,
the dataset should in addition have sufficient expression variance
across cells (several different conditions, time-series experiment over
a dynamic process or a mixture of cell types).

To generate the example data, three iPS cell lines from three humans,
two iPS cell lines from one gorilla and four iPS cell lines from two
cynomolgus macaques were differentiated into neural progenitor cells
(NPCs) over the course of nine days. Samples from six different time
points were sequenced using scRNA-seq. The protocol and timeline were
the same for all species to make the results comparable.

At several steps of the CroCoNet workflow, the replicates need to be
matched to species. To make this possible, a data frame
`replicate2species` has to be created that specifies which species each
replicate belongs to. How this looks like for the example data is shown
below:

``` r
replicate2species <- readRDS("processed_data/replicate2species.rds")
```

| replicate | species    |
|:----------|:-----------|
| H1c1      | human      |
| H2c1      | human      |
| H3c1      | human      |
| G1c1      | gorilla    |
| G1c2      | gorilla    |
| C1c1      | cynomolgus |
| C1c2      | cynomolgus |
| C2c1      | cynomolgus |
| C2c2      | cynomolgus |

### Initial processing

The pipeline assumes that the user has already done the standard
processing steps on the transcriptomic data, including mapping, QC and
normalization. In addition, it is required to have a shared feature
space across all species. In case of closely related species, genes can
be matched by transferring the annotation of one species to the genomes
of others using the tool *Liftoff* \[1\]. It is of course also possible
to use orthology information or even protein language model-based
macrogenes \[2\].

It is also recommended to perform cell type annotation and pseudotime
inference if applicable. After these steps, a SingleCellExperiment (SCE)
object should be created that contains the raw counts, normalized counts
and metadata (including information on species, replicate, cell type
assignment and inferred pseudotime).

In case of the primate differentiation dataset, we created a shared
feature space by transferring the human genome annotation to the gorGor6
and macFas6 genomes via *Liftoff*. We performed the nornalization using
`scran` \[3\], pseudotime inference using `SCORPIUS` \[4\] and cell type
annotation using `singleR` \[5\] with the embryoid body dataset from
Rhodes et al. \[6\] as the reference. We found that the pseudotime
trajectories are reasonably well-aligned across the 3 species,
confirming comparability.

![Pseudotime trajectory of the early primate neural differentiation
dataset coloured by sampling day, pseudotime, species and cell
type.](trajectory.png)

Pseudotime trajectory of the early primate neural differentiation
dataset coloured by sampling day, pseudotime, species and cell type.

## Network inference

Before applying the pipeline, the user must first infer co-expression
networks per replicate using a method of their choice, for example
*GRNBoost2* \[7\] (metric: gradient boosting-based importance scores) or
*correlatePairs* \[3\] (metric: Spearman’s correlation). The networks
should not be restricted to connections between transcriptional
regulators and their target genes, because `CroCoNet` relies also on
connections between the target genes themselves.

If the replicates differ substantially in cell type composition, it is
recommended to downsample them so that the proportions of cell types are
consistent across replicates. This helps disentangle species-specific
effects from cell type differences. It can be beneficial to create
multiple downsampled versions of each replicate and run network
inference on each, ensuring that valuable data are not discarded.

If the chosen network inference algorithm involves stochastic steps (as
is the case for *GRNBoost2*), the results can be made more robust by
running the algorithm multiple times on the same data with different
random seeds.

The network reconstructions should be stored as TSV files where each row
corresponds to an edge. There should be at least 3 columns containing 1)
the 1st gene that forms the edge, 2) the 2nd gene that forms the edge,
and 3) the edge weight. The TSV files should follow the naming
convention *nameOfReplicate.tsv* (if there is a single network per
replicate) or *nameOfReplicate_index.tsv* (if multiple network outputs
are generated per replicate).

In case of example dataset, we inferred networks using GRNBoost2 with
all genes as potential regulators. We run the algorithm 10 times on the
count matrix of each replicate, producing a total of 9×10 TSV files
stored in the directory `network_inference/output/`.  
  

## Loading and processing the networks

In the first part of the workflow, the output files of the network
inference are loaded as `igraph` objects, the edges are filtered and the
edge weights are normalized (if desired), and finally, a phylogeny-aware
consensus network is calculated across all replicates and species. If
the network inference algorithm does not distinguish positively and
negatively correlated gene pairs, this information can also be added at
this stage.

### Loading the networks as igraphs

As the first step, the TSV files containing the network reconstructions
are loaded and summarized as a list of `igraph` \[8\] objects per
replicate with the help of the `loadNetworks` function.

If there are several versions of network reconstruction per replicate
due to different subsamplings or repeated runs (specified via the
parameter `rep`), the edge weights are averaged across these versions,
and a single combined `igraph` object is returned for each replicate. If
the network inference method produces an output with directed edges,
i.e. geneA-geneB and geneB-geneA can both be present (specified via the
parameter `directed`), the edge weights inferred between the same gene
pair but in opposite directions are also averaged. Edges that fall below
the detection threshold in a network version/direction are assigned an
edge weight of 0 for the calculation of the mean. If an edge meets the
detection criteria in too few versions/directions (specified via the
parameter `min_occurrence`), it is assigned a final weight of 0 in the
corresponding replicate, regardless of its mean weight. This filtering
approach helps to 1) denoise the networks and 2) decrease the
computational power needed for the next steps. The steps relying on a
detection threshold are not applicable to correlation-based networks.

In case of the example dataset, the 9×10 TSV files located in the
`network_inference/output/` directory were loaded as a list of 9
`igraph` objects (one network reconstruction for each of the 9
replicates). The parameter `rep` was set to 10 as GRNBoost2 had been run
10 times on each replicate, and the parameter `directed` was left at the
default value of TRUE, as GRNBoost2 produces a directed output. As a
result, the edge weights were averaged across the 10 runs and 2
directions for each edge and replicate. Since the parameter
`min_occurrece` was left at its default value of 2, edges that occurred
zero or one times across all runs and directions were assigned a weight
of 0.

``` r
network_list_raw <- loadNetworks("network_inference/output/", 
                                 replicate_names = replicate2species$replicate,
                                 rep = 10,
                                 directed = TRUE,
                                 min_occurrence = 2)
```

The resulting igraph objects store the weights for each edge:

| from |   to   | weight |
|:----:|:------:|:------:|
| A1BG | ACAT2  | 0.0041 |
| A1BG | ACBD4  | 0.0107 |
| A1BG | ACTR10 | 0.0187 |
| A1BG | ACVR2A | 0.0184 |
| A1BG | ADAM9  | 0.0545 |
| A1BG | ADARB1 | 0.0317 |
| A1BG |  ADM5  | 0.0404 |
| A1BG |  ADO   | 0.0180 |

A few example edges from the raw network of the replicate C2c2

### Removing gene pairs with overlapping annotations

Mapping and counting is problematic for overlapping genomic features:
parts of the reads from one gene can be assigned to the other gene,
leading to correlated expression profiles simply due to genomic
position. This has only a marginal effect on the results of a
differential expression analysis, but can cause false positive edges
with very high edge weights in case of a network analysis. Such
potential artefacts can be circumvented by removing all edges between
genes that have overlapping annotations in the genome(s) of at least one
species.

To do this, the genome annotations are needed in a list format:

``` r
gtf_list <- list(human = plyranges::read_gff("genomes/hg38.gtf"),
                 gorilla = plyranges::read_gff("genomes/gorGor6.gtf"),
                 cynomolgus = plyranges::read_gff("genomes/macFas6.gtf"))
```

Then the function `removeOverlappingGenePairs` determines the genomic
positions of the network genes based on the GTFs, identifies gene pairs
with overlapping annotations in each genome and removes these gene pairs
from all networks:

``` r
network_list_raw_filt <- removeOverlappingGenePairs(network_list_raw, 
                                                    gtf_list = gtf_list, 
                                                    replicate2species = replicate2species, 
                                                    gene_col = "gene_name")
```

During this step, the genomic distance is also added for each gene pair
of the network (the distance is regarded as Inf if the two genes are
located on different chromosomes):

| from |   to   | weight | genomic_dist |
|:----:|:------:|:------:|:------------:|
| A1BG | ACAT2  | 0.0041 |     Inf      |
| A1BG | ACBD4  | 0.0107 |     Inf      |
| A1BG | ACTR10 | 0.0187 |     Inf      |
| A1BG | ACVR2A | 0.0184 |     Inf      |
| A1BG | ADAM9  | 0.0545 |     Inf      |
| A1BG | ADARB1 | 0.0317 |     Inf      |
| A1BG |  ADM5  | 0.0404 |   8270693    |
| A1BG |  ADO   | 0.0180 |     Inf      |

A few example edges from the network of the replicate C2c2 after
removing gene pairs with overlapping annotation and adding genomic
distances

### Normalizing edge weights

Normalizing the edge weights between 0 and 1 is recommended, because it
makes them interpretable as adjacencies and ensures that network
concepts such as connectivity are applicable.

There are 2 approaches for the normalization implemented as part of the
`nornmalizeEdgeWeights` function:

- **Unsigned network** (default): Gene pairs with high negative edge
  weights are considered as connected as gene pairs with high positive
  edge weights. Therefore the negative edge weights are first replaced
  by their absolute values, then all edge weights are scaled by the
  maximum weight across all networks:
  $$a = \frac{|w|}{max\left( |w| \right)}$$ where $a$ is the edge weight
  (adjacency) after normalization and $w$ is the edge weight before
  normalization. After the transformation, the adjacencies around 0
  correspond to the former low positive and low negative values, while
  the adjacencies around 1 correspond to the former high positive and
  high negative values.

- **Signed network**: Gene pairs with high negative edge weights are
  considered unconnected. Therefore all edge weights are transformed
  between 0 and 1 using a min-max normalization:
  $$a = \frac{w - min(w)}{max(w) - min(w)}$$ After the transformation,
  the adjacencies around 0 correspond to the former high negative values
  and the adjacencies around 1 correspond to the former high positive
  values.

In case of the example dataset, we opted for an unsigned network:

``` r
network_list <- normalizeEdgeWeights(network_list_raw_filt,
                                     signed = FALSE)
```

| from |   to   |  weight  | genomic_dist |
|:----:|:------:|:--------:|:------------:|
| A1BG | ACAT2  | 0.000048 |     Inf      |
| A1BG | ACBD4  | 0.000127 |     Inf      |
| A1BG | ACTR10 | 0.000222 |     Inf      |
| A1BG | ACVR2A | 0.000218 |     Inf      |
| A1BG | ADAM9  | 0.000647 |     Inf      |
| A1BG | ADARB1 | 0.000376 |     Inf      |
| A1BG |  ADM5  | 0.000480 |   8270693    |
| A1BG |  ADO   | 0.000214 |     Inf      |

A few example edges from the network of the replicate C2c2 after
normalizing edge weights

Using the resulting adjacencies, we also calculated the connectivity
($k_{i}$) of each gene, which is the sum of adjacencies between the
given gene and all others:

$$k_{i} = \sum\limits_{j \neq i}a_{ij}$$ Both the adjacencies and
connectivites have a highly left-skewed distribution (the adjacencies
range between 0 and 1, with a median of 0.000213, while the
connectivities range between 0 and ~17, with a median of 0.256). This is
expected since connectivities in many biological networks, including
gene regulatory networks, have been reported to follow scale-free
distribution \[9–11\] (although the extent of this varies across systems
and metrics \[12\]).

![The distribution od adjacencies and connectivities for the networks of
the 9 replicates.](edge_weight_con_distr.png)

The distribution od adjacencies and connectivities for the networks of
the 9 replicates.

### Creating the consensus network

As the next step of the workflow, the networks across different
replicates and different species are integrated into a single consensus
network in a phylogeny-aware manner.

For each edge, the consensus adjacency is calculated as the weighted
mean of replicate-wise adjacencies. The weighted mean corrects for 1)
the phylogenetic distances between species, and 2) the different numbers
of replicates per species. As a result, the approach downweighs the edge
weights of the replicates that 1) belong to closely related species, or
2) belong to species with many replicates, so that an imbalanced
sampling across the phylogenetic tree does not bias the consensus
network.

If a phylogeny-aware consensus is desired, the phylogenetic tree has to
be provided. In case of the example dataset, this is the mammalian tree
from Bininda-Emonds et al. \[13\] subsetted for Homo sapiens (human),
Gorilla gorilla (gorilla) and Macaca Fascicularis (cynomolgus):

``` r
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fnature05634/MediaObjects/41586_2007_BFnature05634_MOESM362_ESM.txt", here(wd, "mammalian_tree.txt"))
tree <- read.tree(here(wd, "mammalian_tree.txt"))[["treemammalST_bestDates="]] %>%
  keep.tip(c("Homo_sapiens","Gorilla_gorilla", "Macaca_fascicularis"))
tree$tip.label <- c("cynomolgus", "gorilla", "human")
```

![Phylogenetic tree of the species.](tree.png)

Phylogenetic tree of the species.

Based on the phylogenetic tree and the species assignment of each
replicate, the consensus network is computed using the function
`createConsensus`:

``` r
consensus_network <- createConsensus(network_list, 
                                     replicate2species = replicate2species, 
                                     tree = tree)
```

| from |  to   |  weight   | n_supporting_replicates | supporting_replicates |
|:----:|:-----:|:---------:|:-----------------------:|:---------------------:|
| A1BG | ACAT2 | 0.0000426 |            3            |    H3c1,C2c1,C2c2     |
| A1BG | ACBD4 | 0.0000145 |            1            |         C2c2          |
| A1BG | ACBD6 | 0.0011186 |            1            |         C1c1          |
| A1BG | ACBD7 | 0.0000081 |            1            |         C1c1          |
| A1BG | ACER3 | 0.0000665 |            2            |       H1c1,H3c1       |
| A1BG | ACIN1 | 0.0001691 |            2            |       G1c2,C1c1       |
| A1BG | ACLY  | 0.0000402 |            2            |       H3c1,C2c1       |
| A1BG | ACO1  | 0.0000244 |            1            |         C2c1          |

A few example edges from the consensus network

### Adding the direction of correlation

If the network inference method outputs both positive and negative edge
weights (e.g. correlation-based methods), the directionality of each
edge is saved during the step “Normalizing edge weights”. However, if
the networks are inferred using a method that does not distinguish
between positively and negatively correlated gene pairs, the
directionality can be added in this extra step. The information on
directionality is useful because for certain analyses results are more
meaning if positively and negatively correlated genes are separated
(e.g. eigengene calculation).

The directionality of each edge is determined by calculating a modified
Spearman’s correlation between the expression profiles of the two genes
forming the edge and assigning the sign of the modified Spearman’s rho
as the direction. The calculation uses the approximate Spearman’s rho,
significance testing, and blocking implemented in `correlatePairs`
\[3\].

In case of the example dataset, the networks were inferred using
GRNBoost2 which does not provide information about the direction of
regulation, therefore it is necessary to perform this step. To get the
expression profiles of the network genes, the SCE object needs to be
loaded:

``` r
sce <- readRDS("processed_data/sce.rds")
```

Next, edge directionality is calculated with the function
`addDirectionality`:

``` r
consensus_network <- addDirectionality(consensus_network, sce)
```

| from |  to   |  weight   | n_supporting_replicates | supporting_replicates |   rho   | p.adj | direction |
|:----:|:-----:|:---------:|:-----------------------:|:---------------------:|:-------:|:-----:|:---------:|
| A1BG | ACAT2 | 0.0000426 |            3            |    H3c1,C2c1,C2c2     | 0.0293  | 0.20  |    \+     |
| A1BG | ACBD4 | 0.0000145 |            1            |         C2c2          | -0.0202 | 0.46  |    \-     |
| A1BG | ACBD6 | 0.0011186 |            1            |         C1c1          | 0.0063  | 0.83  |    \+     |
| A1BG | ACBD7 | 0.0000081 |            1            |         C1c1          | -0.0109 | 0.75  |    \-     |
| A1BG | ACER3 | 0.0000665 |            2            |       H1c1,H3c1       | 0.0170  | 0.54  |    \+     |
| A1BG | ACIN1 | 0.0001691 |            2            |       G1c2,C1c1       | 0.0234  | 0.35  |    \+     |
| A1BG | ACLY  | 0.0000402 |            2            |       H3c1,C2c1       | 0.0442  | 0.04  |    \+     |
| A1BG | ACO1  | 0.0000244 |            1            |         C2c1          | 0.0221  | 0.43  |    \+     |

A few example edges from the consensus network after adding the
direction of correlation

  
  

## Module assignment

Once the consensus network is calculated, it can be used to assign
co-expression modules jointly for all species. The center of each module
is a transcriptional regulator and the modules are assigned in 2 main
steps: 1) large initial modules are created by selecting a fixed number
of target genes per regulator, and 2) the initial modules are pruned
using a dynamic filtering approach to keep only the best targets of each
regulator.

### Selecting key transcriptional regulators

First, the transcriptional regulators have to be selected that form the
cores of the co-expression modules. To identify regulators that are
relevant for the biological process under study, it is a useful approach
to combine prior biological knowledge with information derived from the
dataset. This step is implemented by the function `getRegulators`.

As for the prior biological knowledge, we recommend to select only
transcriptional regulators with known binding motifs based on motif
databases. To get this information, the user can choose the following
databases: JASPAR 2024 vertebrate core \[14\], JASPAR 2024 unvalidated
\[14\], the IMAGE database \[15\], or any combination of the above.
Alternatively, the user can provide a custom list of regulators selected
by their preferred method.

As for the information from the dataset, we recommend to select only
transcriptional regulators that are among the highly variable genes in
the data. For defining highly variable genes, `getRegulators` relies on
`scran::getTopHVGs` \[3\]. Briefly, for each species a trend is fitted
between the variance and mean of the log-expression values across all
genes, and the fitted value of a gene is regarded as the technical
component of variation, while the residual from the trend is regarded as
the biological component of variation. All genes with a positive
biological component are selected as highly variable genes in each
species (alternatively, a more stringent variance cutoff, an adjusted
p-value cutoff, a fixed number of genes, or a fixed percentage of genes
can also be defined). As the final set of highly variable genes, the
union of highly variable genes is taken across species.

In case of the example dataset, we used transcriptional regulators
that 1) had an annotated motif in any of the 3 motif databases (JASPAR
2024 vertebrate core, JASPAR 2024 unvalidated and IMAGE) and 2) had a
positive biological component of variance in any of the species:

``` r
regulators <- getRegulators(sce, 
                            source = c("jaspar_core", "jaspar_unvalidated", "image"))
```

This gave us 836 transcriptional regulators that included well-known
pluripotency factors and early neural regulators:

``` r
length(regulators)
#> [1] 836
```

``` r
c("NANOG", "POU5F1", "PAX6", "NEUROD4") %in% regulators
#> [1] TRUE TRUE TRUE TRUE
```

### Assigning initial modules

Once the regulators that provide the starting point of the module
assignment are selected, an initial module can be defined around each of
them using the function `assignInitialModules`. Each of these modules
contain the regulator and its N best target genes (N is recommended to
be between 1000 and 5000). When choosing the best targets, the genes are
ranked based on how strongly they are connected to the regulator
(regulator-target adjacency).

In case of the example dataset, we defined 4000 as the number of genes
in the initial modules:

``` r
initial_modules <- assignInitialModules(consensus_network, 
                                        regulators = regulators, 
                                        N = 4000)
```

| regulator | target | weight  | n_supporting_replicates |          supporting_replicates           |  rho  | p.adj | direction |
|:---------:|:------:|:-------:|:-----------------------:|:----------------------------------------:|:-----:|:-----:|:---------:|
|  POU5F1   | ACAT2  | 0.00080 |            4            |           H1c1,G1c1,C2c1,C2c2            | 0.125 | 1e-12 |    \+     |
|  POU5F1   |  ACCS  | 0.00016 |            4            |           G1c1,G1c2,C1c2,C2c1            | 0.092 | 5e-06 |    \+     |
|  POU5F1   | ACER3  | 0.00053 |            3            |              H1c1,G1c1,C1c2              | 0.090 | 1e-06 |    \+     |
|  POU5F1   |  ACLY  | 0.00081 |            4            |           H1c1,G1c1,C2c1,C2c2            | 0.189 | 2e-27 |    \+     |
|  POU5F1   | ACOT13 | 0.00300 |            6            |      H1c1,H3c1,G1c1,G1c2, C1c2,C2c2      | 0.086 | 7e-08 |    \+     |
|  POU5F1   | ACOT7  | 0.00475 |            8            | H1c1,H2c1,G1c1,G1c2, C1c1,C1c2,C2c1,C2c2 | 0.337 | 7e-96 |    \+     |
|  POU5F1   | ACOT8  | 0.00083 |            3            |              G1c1,C1c1,C2c2              | 0.212 | 3e-35 |    \+     |
|  POU5F1   |  ACP1  | 0.00053 |            3            |              G1c2,C1c2,C2c1              | 0.113 | 9e-11 |    \+     |

Part of the POU5F1 initial module

The resulting modules are obviously too large at this point to be
biologically meaningful. The purpose of this step is rather to select a
uniformly sized pool of genes for all modules to perform the pruning
step on.

### Module pruning

The goal of pruning is to keep only the best targets of each
transcriptional regulator and thereby create biologically meaningful
modules that consist of co-regulated genes with shared cellular
functions.

Which targets are the best can be defined based on 2 metrics:

![Network concepts used for module pruning: regulator-target adjacency
and intramodular connectivity.](adj_regulator_kIM.png)

Network concepts used for module pruning: regulator-target adjacency and
intramodular connectivity.

- **Regulator-target adjacency** ($a_{regulator}$): The adjacency
  between a transcriptional regulator and a member gene in its module.
  It ranges between 0 and 1 with a higher value representing a stronger
  (activating or repressing) regulatory relationship.

  - *Pro*: Pruning based on the regulator-target adjacencies selects the
    genes whose expression profiles co-vary the most with the
    regulator’s expression profile.

  - *Con*: Modules created based on the regulator-target adjacencies
    alone might not fulfill the expectation of being densely connected
    blocks of genes, apart from the link to the regulator, the target
    genes often only have sparse/weak connections within the module
    (danger of false positives).

- **Intramodular connectivity** (*kIM*): The sum of adjacencies between
  a specific gene and all other members of the module. If the module
  size is *n*, it ranges between 0 and *n* - 1 with a higher value
  representing a gene that is more strongly connected to the rest of the
  module.
  $$kIM_{i}^{(q)} = \sum\limits_{\substack{j \in q \\ j \neq i}}a_{ij}$$

  - *Pro*: Pruning based on intramodular connectivities in addition to
    the regulator-target adjacencies ensures that the chosen targets
    co-vary not just with the main regulator but also with the rest of
    the module. These intramodular connections between targets can carry
    important information that are not captured by the links between
    regulator and targets alone (e.g. combinatorial regulation).

  - *Con*: Genes that have a generally high connectivity but do not have
    a particularly strong regulatory relationship with the
    transcriptional regulator might end up in the final module.

The initial targets can be filtered based one or both of these metrics
using 2 approaches:

- **topN**: Keeps a fixed number of targets (default: 50) per regulator.

  - *Pro*: Quick, simple and widespread.

  - *Con*: The approach assigns the same number of target genes to each
    regulator even though biologically speaking different regulators can
    have an effect on different numbers of genes. The number of targets
    is in most cases arbitrary and does not take into account the
    structure of the data at hand.

- **UIK**: Applies dynamic stepwise pruning using knee-point detection.
  In each step, the cumulative sum curve based on the metric of choice
  is calculated per module, the knee point of the curve is identified
  using the Unit Invariant Knee (UIK) method, then only the targets that
  rank higher than the knee point are kept. The modules containing less
  target genes than a minimum module size (default: 20) are removed
  after each pruning step. The steps continue until the percentage of
  removed modules becomes too high (default: \>2%).

  - *Pro*: While setting a minimum module size prevents the modules from
    becoming too small, the exact number of target genes per regulator
    does not have to be pre-defined, in line with the notion that
    different regulators can have an effect on different numbers of
    genes. There are also no hard cutoffs applied to the
    regulator-target adjacencies or intramodular connectivities, but by
    using knee point detection the target genes are filtered in a
    data-driven way.

  - *Con*: The approach still requires a predefined minimum module size.

Compared to clustering-based module assignment approaches, both pruning
methods have the advantage that the modules are allowed to overlap, and
in addition to having its own module, a regulator can be assigned to
another regulator’s module as well. This is in line with the notion that
genes can fulfill several different functions in the cell and gene
regulation can be combinatorial.

In case of both pruning methods, when setting the fixed module size or
the threshold for the minimum module size not just biological but also
technical aspects have to be taken into consideration: for modules
smaller than ~20 genes, the correlation-based preservation statistics in
the next steps might be coupled with high uncertainty, therefore they
should be avoided.

Using the metrics and approaches described above, we implemented 3
methods for module pruning:

- **topN**: Takes a fixed number of targets per regulator with the
  highest regulator-target adjacencies.

- **UIK_adj**: Applies a dynamic stepwise pruning based on the
  regulator-target adjacencies.

- **UIK_adj_kIM**: Applies a dynamic stepwise pruning based on the
  regulator-target adjacencies and intramodular connectivities
  alternately.

In case of the example dataset, we used the method “UIK_adj_kIM”:

``` r
pruned_modules <- pruneModules(initial_modules, 
                               method = "UIK_adj_kIM", 
                               consensus_network = consensus_network)
```

The following figure shows the steps of the pruning process for the
POU5F1 module:

![Stepwise dynamic pruning approach exemplified by the POU5F1 module. In
each step, we calculated the cumulative sum curve of the
regulator-target adjacency or intramodular connectivitiy of the target
genes alternately, then kept the targets that fell below the knee point
of the curve. We continued this process for 3 iterations - until the
module sizes became as small as possible without falling below 20
genes.](POU5F1_detailed_UIK_adj_kIM.png)

Stepwise dynamic pruning approach exemplified by the POU5F1 module. In
each step, we calculated the cumulative sum curve of the
regulator-target adjacency or intramodular connectivitiy of the target
genes alternately, then kept the targets that fell below the knee point
of the curve. We continued this process for 3 iterations - until the
module sizes became as small as possible without falling below 20 genes.

| regulator | target | weight |    kIM    | n_supporting\_ replicates |            supporting\_ replicates             |  rho  | p.adj  | direction | module_size |
|:---------:|:------:|:------:|:---------:|:-------------------------:|:----------------------------------------------:|:-----:|:------:|:---------:|:-----------:|
|  POU5F1   | AP1S2  | 0.026  | 0.4140710 |             8             |   H1c1,H2c1,H3c1, G1c2,C1c1,C1c2, C2c1,C2c2    | 0.38  | 3e-127 |    \+     |     46      |
|  POU5F1   |  APOE  | 0.075  | 1.1619107 |             9             | H1c1,H2c1,H3c1, G1c1,G1c2,C1c1, C1c2,C2c1,C2c2 | 0.58  |   0    |    \+     |     46      |
|  POU5F1   |  CD24  | 0.143  | 1.9220911 |             9             | H1c1,H2c1,H3c1, G1c1,G1c2,C1c1, C1c2,C2c1,C2c2 | 0.63  |   0    |    \+     |     46      |
|  POU5F1   |  CDH2  | 0.036  | 0.6176581 |             8             |   H1c1,H2c1,H3c1, G1c1,G1c2,C1c1, C1c2,C2c2    | -0.46 | 1e-194 |    \-     |     46      |
|  POU5F1   | CENPV  | 0.025  | 0.4156734 |             8             |   H1c1,H3c1,G1c1, G1c2,C1c1,C1c2, C2c1,C2c2    | -0.33 | 6e-93  |    \-     |     46      |
|  POU5F1   | CHCHD2 | 0.031  | 0.5833867 |             7             |      H2c1,H3c1,G1c2, C1c1,C1c2,C2c1, C2c2      | 0.39  | 8e-135 |    \+     |     46      |
|  POU5F1   | CLDN6  | 0.033  | 0.5877413 |             9             | H1c1,H2c1,H3c1, G1c1,G1c2,C1c1, C1c2,C2c1,C2c2 | 0.57  |   0    |    \+     |     46      |
|  POU5F1   | COX7C  | 0.097  | 1.9427046 |             9             | H1c1,H2c1,H3c1, G1c1,G1c2,C1c1, C1c2,C2c1,C2c2 | 0.58  |   0    |    \+     |     46      |

Part of the POU5F1 pruned module

Since we use the “UIK_adj_kIM” approach, the final module sizes are not
predefined but shaped by the data. For this particular dataset, the
number of genes in the pruned modules range from 25 to 63 with a median
of 45:

``` r
plotModuleSizeDistribution(pruned_modules)
```

![Distribution of module sizes after pruning. Median size:
45.](module_size_distribution.png)

Distribution of module sizes after pruning. Median size: 45.

### Visualizing modules

The most important connections of the modules can be plotted using the
function `plotNetworks`. For each module, the resulting plots show the
strongest intramodular connections based on the consensus network and
the genes involved in these connections. The following plot depicts the
top 300 connections of 6 example modules, colored by the direction of
correlation:

``` r
example_modules <- c("NANOG", "POU5F1", "SALL4",  "NEUROD4", "PAX6", "FEZF2")
plotNetworks(example_modules, 
             pruned_modules, 
             consensus_network, 
             color_by = "direction")
```

![The 300 strongest connections of 6 example modules. The thickness of
the edges represents the consensus edge weights and the color of the
edges represents the direction of interaction (positively/negatively
correlated).](network_examples.png)

The 300 strongest connections of 6 example modules. The thickness of the
edges represents the consensus edge weights and the color of the edges
represents the direction of interaction (positively/negatively
correlated).

### Creating random modules

Along with the pruned modules, there is an option to create random
modules as well using the function `createRandomModules`. The random
modules have the same regulators and contain the same number of target
genes as the actual pruned modules, but these target genes are randomly
drawn from all network genes.

``` r
genes <- rownames(sce)
random_modules <- createRandomModules(pruned_modules, 
                                      network_genes = genes)
```

In the next steps of the pipeline, the actual modules are compared to
these random modules in terms of various statistics to check whether the
two groups of modules behave in general differently and to remove those
individual actual modules that show too similar characteristics to the
random modules.

  
  

## Characterizing modules

Since co-expression modules are often more closely associated with a
cellular function than any particular gene alone, it is informative to
investigate which cell types the different modules are active in and
which pathways they are associated to. This can be achieved for example
by calculating module eigengenes and peforming gene set enrichment
analysis.

### Calculating module eigengenes

A concept adapted from WGCNA \[16\], a module eigengene summarizes the
expression profile of an entire module, and it is calculated as the
first principal component of the module expression data. Effectively, it
is a weighted mean of the individual genes’ expression profiles. In
`CroConet`, the calculation of eigengene profiles is implemented by
`calculateEigengenes`. This function outputs other metrics as well for
representing the expression patterns of a module, namely the mean
expression across all module members and the expression of the
regulator.

If a module contains both activated and repressed targets of the
transcriptional regulator, calculating the eigengene (or any other
summary expression profiles) across both directions of regulation does
not make biological sense and leads to the dilution of signal. It makes
more sense to calculate the eigengene either for the activated targets
only (`direction_of_regulation` = “+\_only”, the default) or for the
activated and repressed targets separately (`direction_of_regulation` =
“+-\_separately”). If an eigengene across all targets is desired
irrespective of the direction of regulation, `direction_of_regulation`
should be set to “all_together”.

In case of the example dataset, we calculated module eigengenes for the
activated targets of each module:

``` r
eigengenes <- calculateEigengenes(regulators,
                                  pruned_modules, 
                                  sce,
                                  direction_of_regulation = "+_only")
```

| cell           | species | pseudotime | cell_type      | module    | eigengene | mean_expr | regulator_expr |
|:---------------|:--------|-----------:|:---------------|:----------|----------:|----------:|---------------:|
| ATCTCAGGAAAACT | human   |      0.504 | Early_Ectoderm | POU5F1(+) |    -0.008 |    -0.317 |          0.242 |
| ATCTCAGGAAAGTT | human   |      0.604 | Early_Ectoderm | POU5F1(+) |    -0.013 |    -0.529 |         -0.760 |
| ATCTCAGGAAATTG | human   |      0.570 | Early_Ectoderm | POU5F1(+) |    -0.011 |    -0.422 |         -0.015 |
| ATCTCAGGAAGATT | human   |      0.568 | Early_Ectoderm | POU5F1(+) |    -0.009 |    -0.347 |         -0.760 |
| ATCTCAGGAATACA | human   |      0.525 | Early_Ectoderm | POU5F1(+) |    -0.008 |    -0.275 |         -0.760 |
| ATCTCAGGAATCTT | human   |      0.515 | Early_Ectoderm | POU5F1(+) |    -0.008 |    -0.286 |         -0.760 |
| ATCTCAGGACAATA | human   |      0.576 | Early_Ectoderm | POU5F1(+) |    -0.010 |    -0.405 |         -0.760 |
| ATCTCAGGACCGGC | human   |      0.525 | Early_Ectoderm | POU5F1(+) |    -0.005 |    -0.190 |         -0.760 |

Eigengene profile of the activated POU5F1 module across a couple of
example cells

The eigengene profiles can then be plotted as a heatmap using
`plotEigengeneHeatmap`. Here we focus on the 6 example modules again,
that contain 3 pluripotency regulators (NANOG, POU5F1 and SALL4) and 3
neural regulators (NEUROD4, PAX6 and FEZF2):

``` r
eigengenes_examples <- eigengenes %>%
  dplyr::filter(module %in% paste0(example_modules, "(+)")) %>%
  dplyr::mutate(module = factor(module, paste0(example_modules, "(+)")))
plotEigengeneHeatmap(eigengenes_examples)
```

![The eigengenes of the example modules.](eigengenes_examples.png)

The eigengenes of the example modules.

As expected, the activated targets of the pluripotency regulators show
the highest expression levels at early pseudotime stages and become
downregulated later on, while the activated targets of the neural
regulators have low expression levels at first, and become upregulated
towards later stages of the differentiation process.

### Gene set enrichment analysis

Gene set enrichment analysis offers useful information for the
functional annotation of the modules. Various databases and tools can be
used for such an analysis, including Gene Ontology (GO) \[17, 18\] and
the related R package `topGO` \[19\], the Reactome Pathway Database
\[20\] and the R package `ReactomePA` \[21\], or STRING \[22\] and the R
package `STRINGdb` \[22\]. It is also worth mentioning `clusterProfiler`
\[23\] which provides a universal interface for enrichment analysis
based on various sources and makes it easy to analyze and compare
several different gene sets in a single run.

In the next steps, we perform pathway enrichment analysis on the 6
example modules using the Reactome Pathway Database. First, a data frame
containing the activated targets of the example modules is prepared and
the gene names are converted to ENTREZ IDs using the conversion table
stored in “sym2enrtez.rds” (this was created based on the [GENCODE
“Entrez gene ids” metadata
file](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.metadata.EntrezGene.gz)
corresponding to the human genome version used for mapping).

``` r
sym2entrez <- readRDS("RDS/sym2entrez.rds")
example_modules_4_reactome <- pruned_modules %>%
  dplyr::filter(regulator %in% example_modules & direction == "+") %>%
  dplyr::bind_rows(data.frame(regulator = example_modules,
                              target = example_modules)) %>%
  dplyr::inner_join(sym2entrez, by = c("target" = "gene_name")) %>% 
  dplyr::mutate(regulator = factor(regulator, example_modules))
```

Next, the enrichment analysis is done using the wrapper function
`clusterProfiler::compareCluster` and the testing implemented by
`ReactomePA::enrichPathway`. The universe is defined as all network
genes.

``` r
enrich_res <- clusterProfiler::compareCluster(entrez_id ~ regulator,
                                              "enrichPathway",
                                              example_modules_4_reactome,
                                              universe = unique(sym2entrez$entrez_id))
```

Finally, the results are plotted for all modules together:

``` r
clusterProfiler::dotplot(enrich_res, showCategory = 3)
```

![Reactome enrichment analysis on 6 example
modules.](Reactome_dotplot.png)

Reactome enrichment analysis on 6 example modules.

Encouragingly, the modules of the pluripotency regulators show
enrichment for several pathways related to pluripotency (e.g. *POU5F1
(OCT4), SOX2, NANOG activate genes related to proliferation* and
*Transcriptional regulation of pluripotent stem cells*), while the
modules of the neural regulators show enrichment for several pathways
related to the nerual lineage (e.g. *Axon guidance* and *Nervous system
development*).

  
  

## Module preservation within and across species

In the next part of the CroCoNet workflow, we quantify how well-preseved
the topologies of the modules are between different pairs of replicates,
both within and across species, and summarize these preservation scores
as neighbour-joining trees per module.

### Calculating module preservation statistics

For the calculation of preservation statistics, we use the joint module
assignment derived from the consensus network but compare topological
properties directly between the replicate-wise networks. This way, a
module is always defined as the same set of genes, but the adjacencies
and connectivities of these genes can differ from replicate to
replicate; poorly preserved modules are expected to have many, while
well-preserved modules are expected to have few such topological
differences.

Compared to the standard approach of cross-species module comparison
(defining modules for each species separately and then cross-tabulating
the module member genes across species), topology-based metrics can
detect way more subtle differences. In addition, the results of
cross-tabulation are strongly affected by statistical power: if we find
a low overlap across species, it is difficult to disentangle to what
extent this is due to real biological differences and to what extent
this is just due to having a low power for module detection in the first
place. In contrast, a comparison of network topology for a joint module
assignment is much less dependent on the power.

As part of the function `calculatePresStats`, two preservation
statistics (adapted from WGCNA \[24\]) are implemented to quantify
different aspects of this topological comparison:

- **Correlation of adjacencies** (cor_adj): the correlation of all edge
  weights within the module in the network of replicate1 VS in the
  network of replicate2
  $$cor.adj^{(q)} = {cor}\left( a^{{\lbrack net1\rbrack}{(q)}},a^{{\lbrack net2\rbrack}{(q)}} \right)$$
  where $a^{{\lbrack net\rbrack}{(q)}}$ is the vector of intramodular
  adjacencies in module $q$ and network $net$.

- **Correlation of intramodular connectivities** (cor_kIM): the
  correlation of the intramodular connectivities per module member gene
  in the network of replicate1 VS in the network of replicate2
  $$cor.kIM^{(q)} = {cor}\left( kIM^{{\lbrack net1\rbrack}{(q)}},kIM^{{\lbrack net2\rbrack}{(q)}} \right)$$
  where $kIM^{{\lbrack net\rbrack}{(q)}}$ is the vector of intramodular
  connetivities per gene in module $q$ and network $net$.

Both of these statistics can be calculated with module jackknifing (by
setting `jackknife` to TRUE). This means that the function creates all
possible jackknifed versions of each input module by removing each
target gene assigned to that module (the regulator itself is never
excluded), then it calculates the preservation statistics for all of
these jackknifed module versions in addition to the original module. The
advantage of using jackknfing is that a confidence interval can be
calculated for each module and statistic instead of just getting a
single value. Later on in the pipeline, jackknifing can also provide
information about which target genes within a conserved/diverged module
are particularly responsible for the conservation/divergence. If
jackknifing is not desired, the user can just set `jackknife` to FALSE
which will also substantially reduce running times.

In case of the example dataset, we calculated both preservation
statistics and opted for using jackknifing:

``` r
pres_stats_jk <- calculatePresStats(pruned_modules, 
                                    network_list,
                                    c("cor_adj", "cor_kIM"), 
                                    replicate2species,
                                    jackknifing = TRUE)
random_pres_stats_jk <- calculatePresStats(random_modules, 
                                           network_list,
                                           c("cor_adj", "cor_kIM"), 
                                           replicate2species,
                                           jackknifing = TRUE)
```

| regulator | type |       id       | gene\_ removed | replicate1 | replicate2 | species1 | species2 | cor_adj | cor_kIM |
|:---------:|:----:|:--------------:|:--------------:|:----------:|:----------:|:--------:|:--------:|:-------:|:-------:|
|  POU5F1   | orig |  POU5F1_orig   |       NA       |    H1c1    |    H2c1    |  human   |  human   |  0.589  |  0.726  |
|  POU5F1   | orig |  POU5F1_orig   |       NA       |    H1c1    |    H3c1    |  human   |  human   |  0.543  |  0.696  |
|  POU5F1   | orig |  POU5F1_orig   |       NA       |    H1c1    |    G1c1    |  human   | gorilla  |  0.474  |  0.644  |
|  POU5F1   | orig |  POU5F1_orig   |       NA       |    H1c1    |    G1c2    |  human   | gorilla  |  0.544  |  0.676  |
|  POU5F1   |  jk  | POU5F1_jk_APOE |      APOE      |    H1c1    |    H2c1    |  human   |  human   |  0.585  |  0.724  |
|  POU5F1   |  jk  | POU5F1_jk_APOE |      APOE      |    H1c1    |    H3c1    |  human   |  human   |  0.537  |  0.688  |
|  POU5F1   |  jk  | POU5F1_jk_APOE |      APOE      |    H1c1    |    G1c1    |  human   | gorilla  |  0.489  |  0.666  |
|  POU5F1   |  jk  | POU5F1_jk_APOE |      APOE      |    H1c1    |    G1c2    |  human   | gorilla  |  0.556  |  0.697  |

Preservation statistics of the POU5F1 module and one of its jackknifed
module versions for a couple of example replicate pairs

If the preservation statistics are calculated with jackknifing, a
statistic of interest can be estimated for a module as a whole by taking
its median across all jackknifed versions along with the confidence
interval if the median. This summarization is performed by the function
`summarizeJackknifePresStats`:

``` r
pres_stats <- summarizeJackknifeStats(pres_stats_jk)
random_pres_stats <- summarizeJackknifeStats(random_pres_stats_jk)
```

| regulator | replicate1 | replicate2 | species1 | species2 | cor_adj | var_cor_adj |
|:---------:|:----------:|:----------:|:--------:|:--------:|:-------:|:-----------:|
|  POU5F1   |    H1c1    |    H2c1    |  human   |  human   |  0.590  |  0.000113   |
|  POU5F1   |    H1c1    |    H3c1    |  human   |  human   |  0.544  |  0.000154   |
|  POU5F1   |    H1c1    |    G1c1    |  human   | gorilla  |  0.473  |  0.000251   |
|  POU5F1   |    H1c1    |    G1c2    |  human   | gorilla  |  0.545  |  0.000236   |

Preservation statistics of the POU5F1 module for a couple of example
replicate pairs, after summarizing the values of all jackknifed module
versions

| lwr_cor_adj | upr_cor_adj | cor_kIM | var_cor_kIM | lwr_cor_kIM | upr_cor_kIM |
|:-----------:|:-----------:|:-------:|:-----------:|:-----------:|:-----------:|
|    0.585    |    0.593    |  0.727  |  0.000224   |    0.724    |    0.731    |
|    0.539    |    0.546    |  0.695  |  0.000284   |    0.688    |    0.700    |
|    0.470    |    0.478    |  0.641  |  0.000327   |    0.638    |    0.646    |
|    0.540    |    0.547    |  0.673  |  0.000410   |    0.668    |    0.677    |

There are two expectations that a statistic is expected to fulfill if it
is a good measure of module preservation:

- The scores are expected to be higher for the actual modules than for
  the random modules.

- The scores are expected to decrease with increasing phylogenetic
  distance between the replicates compared.

As a visual check, the function `plotPresStatDistributions` plots the
distributions of the preservation scores for the actual and random
modules and for within-species and cross-species replicate pairs,
whereas the function `plotPresStats` plots the cross-species scores
against the within-species scores per species pair:

``` r
plotPresStatDistributions(pres_stats, 
                          random_pres_stats, 
                          c("cor_adj", "cor_kIM"))
```

![Distribution of preservation statistics within and across species, for
the actual and random modules.](pres_stat_distributions.png)

Distribution of preservation statistics within and across species, for
the actual and random modules.

``` r
plotPresStats(pres_stats, 
              random_pres_stats, 
              c("cor_adj", "cor_kIM"))
```

![Cross-species VS within-species preservation scores per species pair
for the actual and random modules.](pres_stat_across_vs_within.png)

Cross-species VS within-species preservation scores per species pair for
the actual and random modules.

The plots show that when applied to the example dataset, both *cor.kIM*
and *cor.adj* fulfill the expectations. To quantitatively compare which
of the two statistics performs better, the function `comparePresStats`
plots 1) the difference in preservation between each actual and the
corresponding random module, and 2) the inverse Pearson’s correlation
(*-r*) between the preservation scores and phylogenetic distance:

``` r
comparePresStats(pres_stats, 
                 random_pres_stats, 
                 tree)
```

![Difference of the \<em\>cor.kIM\</em\> / \<em\>cor.adj\</em\> scores
between each actual and the corresponding random module (left) and
inverted Pearson's correlation (-\<em\>r\</em\> ) between
\<em\>cor.kIM\</em\> / \<em\>cor.adj\</em\> and phylogenetic distance
(right).](pres_stat_comparison.png)

Difference of the *cor.kIM* / *cor.adj* scores between each actual and
the corresponding random module (left) and inverted Pearson’s
correlation (-*r* ) between *cor.kIM* / *cor.adj* and phylogenetic
distance (right).

Since the difference between the actual and random modules as well as
the strength of the phylogenetic signal are higher for *cor.kIM*, we use
this statistic for the rest of the analysis.

### Converting preservation scores to distance measures

As a preparation for tree reconstruction, the preservation statistics
need to be converted to distance measures. This is done by the function
`convertPresToDist` using the formula: $$d = \frac{1 - p}{2},$$ where
$d$ is the distance measure and $p$ is the preservation statistic. The
resulting distances range between 0 and 1: a distance of 0 corresponds
to a preservation score of 1, while a distance of 1 corresponds to a
preservation score of -1.

``` r
dist_jk <- convertPresToDist(pres_stats_jk, "cor_kIM")
random_dist_jk <- convertPresToDist(random_pres_stats_jk, "cor_kIM")
```

The resulting object is a list, where each list element stores the
distance measures of a module (an original or jackknifed module version,
if jackknifing is used) in a data frame format:

| regulator | type |       id       | gene_removed | replicate1 | replicate2 | species1 | species2 | dist  |
|:---------:|:----:|:--------------:|:------------:|:----------:|:----------:|:--------:|:--------:|:-----:|
|  POU5F1   |  jk  | POU5F1_jk_APOE |     APOE     |    H1c1    |    H2c1    |  human   |  human   | 0.138 |
|  POU5F1   |  jk  | POU5F1_jk_APOE |     APOE     |    H1c1    |    H3c1    |  human   |  human   | 0.156 |
|  POU5F1   |  jk  | POU5F1_jk_APOE |     APOE     |    H1c1    |    G1c1    |  human   | gorilla  | 0.167 |
|  POU5F1   |  jk  | POU5F1_jk_APOE |     APOE     |    H1c1    |    G1c2    |  human   | gorilla  | 0.151 |

Distance measures of the jackknifed POU5F1 module version
“POU5F1_jk_APOE” for a couple of example replicate pairs

The distance measures of a module can be organized into and plotted as a
distance matrix. For visualization, it is usually the simplest and most
informative to focus on the original (not jackknifed) modules only:

``` r
dist <- dist_jk[paste0(regulators, "_orig")]
names(dist) <- regulators
```

The function `plotDistMats` plots the distance matrix of a single module
(required input: data frame) or the distance matrices of several modules
together (required input: list of data frames), for example it can be
used to visualize the six example modules side by side:

``` r
plotDistMats(dist[example_modules])
```

![The distance matrices of 6 example modules.](dist_mat_examples.png)

The distance matrices of 6 example modules.

It is apparent based on these six examples that the distance patterns
vary quite a bit from module to module. Some modules (e.g. FEZF2) have
in general low distance values which points to a high degree of
preservation, while other modules (e.g. NEUROD4) tend to have higher
distance values which points to a worse preservation. While for most
modules the within-species distances are overall lower than the
cross-species distances (this was a selection criteria for the
preservation statistic after all), the magnitude of this difference is
not constant: in case of the SALL4 module the within-species and
cross-species scores span by and large a similar range, while in case of
the POU5F1 module the ape-cynomolgus cross-species scores are
drastically higher than the rest of the scores.

### Tree reconstruction

To be able to quantitatively compare these different patterns of
preservation, the distance measures are summarized as trees in the next
step of the workflow. Each of the resulting trees stands for a single
module (an original or jackknifed module version, if jackknifing is
used): the tips represent the replicates and the branch lengths
represent the dissimilarity of module topology between the networks of 2
replicates.

The trees are reconstructed using the neighbor-joining algorithm
implemented by `ape:nj`. The function `reconstructTrees` is a wrapper
around this, tailored for the specific purpose:

``` r
trees_jk <- reconstructTrees(dist_jk)
random_trees_jk <- reconstructTrees(random_dist_jk)
```

The output is a list of `phylo` objects per (jackknifed or original)
module.

For visualizaton, it is again recommended to retrieve the trees of the
original modules only:

``` r
trees <- trees_jk[paste0(regulators, "_orig")]
names(trees) <- regulators
```

Based on this list, the function `plotTrees` can visualize one or more
module trees analogously to `plotDistMats`. Here we focus on the trees
representations of the six example modules again:

``` r
spec_colors <- c(human = "#4DAF4A", gorilla = "#377EB8", cynomolgus = "#9a1ebd")
plotTrees(trees[example_modules], species_colors = spec_colors)
```

![The tree reconstructions of 6 example modules.](tree_examples.png)

The tree reconstructions of 6 example modules.

In line with the differences among the distance matrices, the trees also
differ from module to module. FEZF2, the module with generally low
distance measures (i.e. high preservation scores), tends to have short
branch lengths in its tree, while NEUROD4, the module with generally
high distance measures (i.e. low preservation scores), has a large tree
with long branches. In case of the SALL4 module that had similar
within-species and cross-species distance measures, the replicates of
different species are intermixed within the tree, while in case of the
POU5F1 module that has much higher ape-cynomolgus distances than ape-ape
or cynomolgus-cynomolgus distances, the ape and cynomolgus replicates
are clearly separated within the tree and connected by a long internal
branch.

### Calculating tree-based statistics

To quantify these tell-tale characteristics of the module trees, various
tree-based statistics can be calculated:

- **Total tree length**: The sum of all branch lengths in the tree;
  measures module variability both within and across species.

- **Diversity of a species**: The total length of the branches
  connecting the replicates of the given species to each other; measures
  module variability within this particular species.

- **Within-species diversity**: The sum of the diversity values across
  all species; measures module variability within species in general.

- **Monophyleticity of a species**: Indicates whether the tree is
  monophyletic for the replicates of the given species. Only if a module
  tree is monophyletic for a species of interest can the module be
  tested for divergence between this species and all others.

- **Subtree length of a species:** The sum of the branch lengths in the
  subtree that is defined by the replicates of the species and includes
  the internal branch connecting these replicates to the rest of the
  tree (i.e. the lineage of interest). Undefined if the tree is not
  monophyletic for the species of interest.

![Tree-based statistics to characterize the cross-species conservation
of network modules.](tree_based_stats.png)

Tree-based statistics to characterize the cross-species conservation of
network modules.

These statistics are all implemented in the function
`calculateTreeStats`:

``` r
tree_stats_jk <- calculateTreeStats(trees_jk)
random_tree_stats_jk <- calculateTreeStats(random_trees_jk)
```

| regulator | type |       id        | gene\_ removed | total\_ tree_length | within_species\_ diversity | human\_ diversity | gorilla\_ diversity | cynomolgus\_ diversity |
|:---------:|:----:|:---------------:|:--------------:|:-------------------:|:--------------------------:|:-----------------:|:-------------------:|:----------------------:|
|  POU5F1   | orig |   POU5F1_orig   |       NA       |        1.038        |           0.728            |       0.221       |        0.034        |         0.473          |
|  POU5F1   |  jk  | POU5F1_jk_AP1S2 |     AP1S2      |        1.042        |           0.724            |       0.220       |        0.034        |         0.470          |
|  POU5F1   |  jk  | POU5F1_jk_APOE  |      APOE      |        1.036        |           0.734            |       0.221       |        0.036        |         0.476          |
|  POU5F1   |  jk  | POU5F1_jk_CD24  |      CD24      |        1.038        |           0.733            |       0.225       |        0.039        |         0.469          |

Tree-based statistics of the POU5F1 module and some of its jackknifed
module versions

| human\_ monophyl | gorilla\_ monophyl | cynomolgus\_ monophyl | human\_ subtree_length | gorilla\_ subtree_length | cynomolgus\_ subtree_length |
|:----------------:|:------------------:|:---------------------:|:----------------------:|:------------------------:|:---------------------------:|
|       TRUE       |        TRUE        |         TRUE          |         0.258          |          0.120           |            0.661            |
|       TRUE       |        TRUE        |         TRUE          |         0.257          |          0.123           |            0.662            |
|       TRUE       |        TRUE        |         TRUE          |         0.256          |          0.113           |            0.666            |
|       TRUE       |        TRUE        |         TRUE          |         0.261          |          0.128           |            0.649            |

Similarly to the preservation statistics, the tree-based statistics can
also be summarized over all jackknifed version of a module to obtain an
overall estimate along with a confidence interval:

``` r
tree_stats <- summarizeJackknifeStats(tree_stats_jk, 
                                      c("total_tree_length", 
                                        "within_species_diversity", 
                                        "human_branch_length", 
                                        "human_diversity"))
random_tree_stats <- summarizeJackknifeStats(random_tree_stats_jk, 
                                             c("total_tree_length", 
                                               "within_species_diversity", 
                                               "human_branch_length", 
                                               "human_diversity"))
```

  
  

## Quantifying cross-species conservation of network modules

After having calculated first pairwise preservation scores between
replicates and then higher-level tree-based statistics to characterize
how different the topologies of the modules are within and across
species, we are ready to quantify the evolutionary conservation of
modules. But why do we need the within-species differences in the first
place? Why is the extent of cross-species differences not directly the
answer to evolutionary conservation?

There are several sources of module topological differences when we
compare two species: there is of course true evolutionary divergence
between the species which is what we are primarily interested in, but
this is convoluted with genetic diversity across individuals of the same
species, environmental factors and technical noise. In the end, if we
find that the topology of a module is different between two species, it
is difficult to say to what extent it is actually due to divergence.
This is why the within-species differences also become valuable sources
of information: they provide an estimate of the variability originating
from confounding sources and thus can be used to discern the signal of
real divergence.

### Overall conservation and divergence

To apply this concept for the quantification of overall conservation
across all species, we used two tree-based statistics: within-species
diversity and total tree length. The two metrics are based on the same
trees and calculated using in part the same branches, so they are
clearly not independent. Their general relationship can be captures by
linear regression, that in biological terms describes which part of the
total module variability is due to the confounding factors (genetic
inter-individual differences, environmental and technical noise).

#### Filtering the module trees

Where a module is located along the regression line informs us about
detection robustness. Modules that have both low within-species
diversity and low total tree length are well-preserved both within and
across species, meaning that these modules could be robustly detected in
all replicates, whereas modules, for which both metrics are high, are
poorly preserved not just across but also within species, indicating a
high detection uncertainty (“wobbliness”).

The general trend between the within-species diversities and total tree
lengths can be visualized using the function `plotTreeStats`. In the
plot below, both the actual and the random modules are included:

``` r
plotTreeStats(tree_stats, 
              random_tree_stats, 
              c("within_species_diversity", "total_tree_length"))
```

![The relationship of the total tree lengths and within-species
diversities across all actual and random
modules.](within_species_diversity_total_tree_length.png)

The relationship of the total tree lengths and within-species
diversities across all actual and random modules.

As we have seen before, the within-species diversities and total tree
lengths are in general higher for the random modules than for the actual
modules - thus the random modules cluster separately on the plot at the
“wobbly” end of the spectrum. However, there are still individual actual
modules that have particularly high total tree lengths and
within-species diversities, comparable with those of the random modules.
Since the tree properties of these modules are likely dominated by
technical and environmental noise, it is recommended to remove them for
further analysis.

The function `filterModuleTrees` calculates how probable it is that the
statistics of an actual comes from the distribution of all actual
modules ($p_{actual}$) and how probable it is that the statistics comes
from the distribution of all random modules ($p_{random}$) using the
probability density functions of the two bivariate normal distributions.
If

$$\frac{p_{actual}}{p_{actual} + p_{random}} > 0.95$$

is fulfilled for both tree statistics, the module is kept, otherwise it
is removed.

``` r
tree_stats_filt <- filterModuleTrees(tree_stats, random_tree_stats)
```

In case of the example dataset, this filtering removes 10 modules
(marked in red):

``` r
plotTreeStats(tree_stats, 
              random_tree_stats, 
              c("within_species_diversity", "total_tree_length")) +
  geom_point(data = dplyr::anti_join(tree_stats, tree_stats_filt), color = "red3", size = 1.2)
```

![The relationship of the total tree lengths and within-species
diversities across all actual and random modules, with the removed
actual modules
marked.](within_species_diversity_total_tree_length_filt.png)

The relationship of the total tree lengths and within-species
diversities across all actual and random modules, with the removed
actual modules marked.

#### Finding conserved and diverged modules

After removing these “wobbly” modules, cross-species conservation can be
quantified using the residuals of the linear regression model. These
residuals from the regression line reflect the variance of the total
tree length that cannot be explained by within-species diversity,
i.e. the contribution of true evolutionary divergence to the total
module variability. Positive residuals correspond to a
greater-than-expected contribution and negative residuals to a
lower-than-expected contributions, respectively. In order to identify
modules with significant deviations from the overall trend, we calculate
the 95% prediction interval of the linear fit and considered modules
above the upper bound of the prediction interval diverged and the
modules below the lower bound of the prediction interval conserved.

![The concept of using tree-based statistics to puzzle apart conserved,
diverged, robust and wobbly modules.](within_total_concept.png)

The concept of using tree-based statistics to puzzle apart conserved,
diverged, robust and wobbly modules.

The function `fitTreeStatsLm` fits a regression model between the two
statistics and the function `findConservedDivergedModules` identifies
conserved and diverged modules by comparing each module’s residual to
the 95% prediction interval of the fit.

The linear regression can be weighted by the error of the data points
derived from jackknifing. If this option is chosen (`weighted_lm` is set
to TRUE), the weight of a module in the regression is defined as
inversely proportional to the total tree length variance across all of
its jackknifed module versions.

In case of the example dataset, we used a weighted linear regression:

``` r
lm_overall <- fitTreeStatsLm(tree_stats_filt, 
                             focus = "overall",
                             weighted_lm = TRUE)
module_conservation_overall <- findConservedDivergedModules(tree_stats_filt, 
                                                            lm_overall)
```

|  focus  | regulator | module\_ size | total\_ tree_length | lwr_total\_ tree_length | upr_total\_ tree_length | within_species\_ diversity | lwr_within\_ species_diversity |
|:-------:|:---------:|:-------------:|:-------------------:|:-----------------------:|:-----------------------:|:--------------------------:|:------------------------------:|
| overall |   NANOG   |      42       |        0.599        |          0.596          |          0.604          |           0.457            |             0.456              |
| overall |  POU5F1   |      46       |        1.036        |          1.033          |          1.042          |           0.730            |             0.726              |
| overall |   SALL4   |      40       |        1.130        |          1.120          |          1.143          |           1.191            |             1.173              |
| overall |  NEUROD4  |      52       |        1.208        |          1.200          |          1.214          |           1.254            |             1.240              |
| overall |   PAX6    |      48       |        0.798        |          0.794          |          0.804          |           0.643            |             0.639              |
| overall |   FEZF2   |      40       |        0.560        |          0.554          |          0.567          |           0.434            |             0.426              |

Cross-species conservation measures for 6 example modules

| upr_within\_ species_diversity |  fit  | lwr_fit | upr_fit | residual | weight | t_score |  conservation   |
|:------------------------------:|:-----:|:-------:|:-------:|:--------:|:------:|:-------:|:---------------:|
|             0.461              | 0.654 |  0.571  |  0.737  |  -0.055  | 0.003  | -1.296  | not_significant |
|             0.733              | 0.881 |  0.815  |  0.946  |  0.156   | 0.004  |  4.650  |    diverged     |
|             1.204              | 1.264 |  1.092  |  1.437  |  -0.134  | 0.001  | -1.533  | not_significant |
|             1.262              | 1.317 |  1.178  |  1.456  |  -0.109  | 0.001  | -1.545  | not_significant |
|             0.647              | 0.808 |  0.733  |  0.884  |  -0.011  | 0.003  | -0.274  | not_significant |
|             0.442              | 0.635 |  0.523  |  0.746  |  -0.075  | 0.001  | -1.324  | not_significant |

Using this approach, we identified 20 conserved and 24 diverged modules:

``` r
table(module_conservation_overall$conservation)
#> 
#>       conserved        diverged not_significant 
#>              21              23             788
```

The function `plotConservedDivergedModules` helps to visualize the
regression line, its 95% prediction interval (grey area), and the
conserved and diverged modules (colored green and red, respectively). It
also labels the top N most conserved and most diverged modules. In the
plot below, the top 5 most conserved and most diverged modules based on
residuals are highlighted:

``` r
plotConservedDivergedModules(module_conservation_overall,
                             N = 5)
```

![Weighted linear regression between the total tree length and
within-species diversity, with the top 5 most conserved and most
diverged modules marked.](module_conservation_overall.png)

Weighted linear regression between the total tree length and
within-species diversity, with the top 5 most conserved and most
diverged modules marked.

To illustrate what makes these modules conserved or diverged, their
distance matrices and module trees can also be plotted:

``` r
top5_cons_div_modules <- bind_rows(module_conservation_overall %>%
                                    dplyr::filter(conservation == "conserved") %>%
                                    dplyr::slice_min(order_by = residual, n = 5),
                                   module_conservation_overall %>%
                                    dplyr::filter(conservation == "diverged") %>%
                                    dplyr::slice_max(order_by = residual, n = 5)) %>% 
  dplyr::pull(regulator)
```

``` r
plotDistMats(dist[top5_cons_div_modules], ncol = 5)
```

![The distance matrices of the 5 most conserved (top row) and 5 most
diverged (bottom row) modules.](dist_mat_cons_div.png)

The distance matrices of the 5 most conserved (top row) and 5 most
diverged (bottom row) modules.

``` r
plotTrees(trees[top5_cons_div_modules], species_colors = spec_colors, ncol = 5)
```

![The tree reconstructions of the 5 most conserved (top row) and 5 most
diverged (bottom row) modules.](trees_cons_div.png)

The tree reconstructions of the 5 most conserved (top row) and 5 most
diverged (bottom row) modules.

The contrast between the conserved and diverged modules in terms of
their distance matrices and module trees is clearly visible. The
conserved modules tend to have balanced distance matrices with similarly
low distance measures for all within-species and cross-species
comparisons, and correspondingly, tree reconstructions without any
species separation, with replicates from the same species being as close
to each other as replicates from different species.

In contrast, diverged modules tend to have much higher distance measures
across species than within species for 1 or more species pairs (in case
of the example modules it is always the human-cynomolgus and
gorilla-cynomolgus species pairs). This also results in trees where
these species are strikingly distinct, with short branches within the
species or clades and a long internal branch separating these subtrees
from each other.

#### Finding conserved and diverged edges

A divergence score can be calculated not just for a module as a whole,
but also for each edge within a module, which can pinpoint individual
edges or regions of neighboring edges that are particularly conserved or
diverged. This score is calculated based on the different weights of a
given edge in the networks of the different replicates. The edge weights
are compared across species using an ANOVA, and the log₁₀ of the
F-statistic (i.e. between-species variability / within-species
variability) is regarded as the measure of edge divergence.

These scores can be returned as a data frame using
`calculateEdgeDivergence` or directly plotted using `plotNetworks` with
`color_by` set to “edge_divergence”. Below we plot the 300 strongest
edges within each of the 5 most conserved and most diverged modules
colored by edge divergence:

``` r
plotNetworks(top5_cons_div_modules, 
             pruned_modules, 
             consensus_network, 
             network_list, 
             replicate2species, 
             color_by = "edge_divergence",
             ncol = 5)
```

![The 300 strongest connections of the 5 most conserved (top row) and 5
most diverged (bottom row) modules. The thickness of the edges
represents the consensus edge weights and the color of the edges
represents how different the mean edge weights are across the 3
species.](networks_cons_div.png)

The 300 strongest connections of the 5 most conserved (top row) and 5
most diverged (bottom row) modules. The thickness of the edges
represents the consensus edge weights and the color of the edges
represents how different the mean edge weights are across the 3 species.

Based on the network visualizations, the diverged modules (especially
HMGA1 and POU5F1) tend to have more diverged edges than the conserved
modules, which confirms the previous findings from a different aspect.
The genes that form the most diverged edges within a module can also be
pinned down, for example in case of the POU5F1 module, SCGB3A2 is the
main hub for these diverged edges.

#### Finding conserved and diverged target genes

Genes that contribute the most to the conservation/divergence of a
module can be identified not just using the edge divergence scores, but
also based on jackknifing.

During jackknifing, each member gene of a module is removed and all
statistics are recalculated, including the tree-based statistics (total
tree length and within-species diversity) that inform us about
cross-species conservation. Our working hypothesis is that if removing a
target gene from a diverged module makes the module more conserved, then
that target was responsible for divergence in the original module, and
vice versa, if removing a target gene from a conserved module makes the
module more diverged, then that target was responsible for conservation
in the original module.

To quantify these effects, the function `findConservedDivergedTargets`
makes use of the linear model that was fitted across all modules between
the total tree length and within-species diversity and that was used to
identify the conserved and diverged modules in the first place. However,
in this case we calculated a residual for each jackknifed module version
separately (${\widehat{\varepsilon}}_{i}$) and compared these residuals
to the residual of the original, complete module
(${\widehat{\varepsilon}}_{\text{orig}}$) in order to derive a target
gene contribution score ($c_{i}$):

$$c_{i} = \frac{{\widehat{\varepsilon}}_{\text{orig}} - {\widehat{\varepsilon}}_{i}}{{\widehat{\varepsilon}}_{\text{orig}}}$$

A positive score indicates that removing the target gene weakens the
signal of conservation or divergence, while a negative score means that
after jackknifing the signal becomes stronger. Thus, the most diverged
targets in a diverged module and the most conserved targets in a
conserved module both get high positive scores. In extreme cases, the
score could reach or exceed 1, which means that the removal of a single
target gene would abolish or even reverse the signal. However, in most
cases the detected signal is expected to originate from the joint
contribution of many genes (this was the case for all modules except one
in our example dataset).

We calculated the target contribution scores for the five most conserved
and five most diverged modules:

``` r
target_contributions_overall <- foreach(module = top5_cons_div_modules,
                                        .combine = bind_rows) %do% {
                                  
                                          findConservedDivergedTargets(module, tree_stats_jk, lm_overall)
                                          
                                        }
```

|  focus  | regulator | type |       id        | gene\_ removed | within_species\_ diversity |
|:-------:|:---------:|:----:|:---------------:|:--------------:|:--------------------------:|
| overall |  POU5F1   | orig |   POU5F1_orig   |       NA       |           0.728            |
| overall |  POU5F1   |  jk  | POU5F1_jk_AP1S2 |     AP1S2      |           0.724            |
| overall |  POU5F1   |  jk  | POU5F1_jk_APOE  |      APOE      |           0.734            |
| overall |  POU5F1   |  jk  | POU5F1_jk_CD24  |      CD24      |           0.733            |

Cross-species conservation measures for a few example target genes of
POU5F1

| total\_ tree_length |  fit  | lwr_fit | upr_fit | residual | contribution |
|:-------------------:|:-----:|:-------:|:-------:|:--------:|:------------:|
|        1.038        | 0.879 |  0.813  |  0.945  |  0.159   |      NA      |
|        1.042        | 0.876 |  0.810  |  0.942  |  0.166   |    -0.045    |
|        1.036        | 0.884 |  0.818  |  0.949  |  0.152   |    0.041     |
|        1.038        | 0.883 |  0.817  |  0.949  |  0.155   |    0.023     |

We plotted the scores to identify the most impactful target genes in
each module:

``` r
target_contributions_overall %>% 
  dplyr::filter(type == "jk") %>% 
  inner_join(module_conservation_overall %>% dplyr::select(regulator, conservation)) %>% 
  dplyr::mutate(regulator = factor(regulator, top5_cons_div_modules)) %>% 
  group_by(regulator) %>% 
  dplyr::mutate(to_label = ifelse(contribution %in% sort(contribution, decreasing = TRUE)[1:2], unique(conservation), "no_label")) %>% 
  ggplot(aes(x = regulator, y = contribution, color = to_label)) +
  geom_quasirandom(aes(size = to_label)) +
  theme_bw() +
  scale_size_manual(values = c("diverged" = 0.8, "conserved" = 0.8, "no_label" = 0.5), guide = "none") +
  scale_color_manual(values = c("conserved" = "#2B823A", "diverged" = "#AA4139", "no_label" = "black"), guide = "none") +
  geom_label_repel(data = . %>%
                     dplyr::filter(to_label != "no_label"),
                   ggplot2::aes(label = gene_removed),
                   fill = "white", size = 2.5, label.size = 0.08, segment.size = 0.08, box.padding = 0.05, label.padding = 0.05) +
  ylab("target gene contribution") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10.5, color = "black"))
```

![Contribution scores of the target genes in the five most conserved (in
green) and five most diverged (red) modules. The two top-contributing
targets are labeled for each module.](target_contributions_overall.png)

Contribution scores of the target genes in the five most conserved (in
green) and five most diverged (red) modules. The two top-contributing
targets are labeled for each module.

We took a closer look at the POU5F1 module and plotted the original and
jackknifed versions along the regression line with the help of the
function `plotConservedDivergedTargets`:

``` r
target_contribution_POU5F1 <- target_contributions_overall %>% 
  dplyr::filter(regulator == "POU5F1")
```

``` r
plotConservedDivergedTargets(target_contribution_POU5F1)
```

![The five most diverged target genes of POU5F1 identified using
tree-based statistics in combination with
jackknifing.](target_contributions_POU5F1.png)

The five most diverged target genes of POU5F1 identified using
tree-based statistics in combination with jackknifing.

All jackknifed versions are located way outside the prediction interval,
meaning that the whole module exhibits a strong signal of divergence
that cannot be explained by any of the target genes alone. The most
diverged target with the highest contribution score is SCGB3A2, in line
with the findings based on edge divergence.

To investigate whether the network divergence/conservation also
translates into expression divergence/conservation, the expression
profiles of the regulator and its targets can be plotted per species
using `plotExprAlongPseudotime` (for pseudotime trajectories) or
`plotExprPerCellType` (for distinct cell types). Below, we plot the
expression of POU5F1 and its 5 most diverged target genes across the
pseudotime range of the neural differentiation process:

``` r
top5_diverged_POU5F1_targets <- target_contributions_POU5F1 %>%
  slice_min(order_by = residual, n = 5) %>%
  pull(gene_removed)
```

``` r
ct_colors <- c(Pluripotent_Cells = "#86C1E6", Early_Ectoderm = "#F4AB62", Neurons = "#CA6102")
plotExprAlongPseudotime(c("POU5F1", top5_diverged_POU5F1_targets), 
                        sce, 
                        species_colors = spec_colors, 
                        cell_type_colors = ct_colors)
```

![The expression profiles of POU5F1 and its five most diverged target
genes across the pseudotime trajectory per
species.](target_expression_POU5F1.png)

The expression profiles of POU5F1 and its five most diverged target
genes across the pseudotime trajectory per species.

Several target genes of POU5F1 that are diverged in terms of network
topology (SCGB3A2, SPP1 and FTL) also show different expression patterns
across species. Such differences cannot be observed for POU5F1, the
regulator itself — a hint toward cis-regulatory changes. Indeed, after a
more thorough investigation we found that both SCGB3A2 and SPP1 are
known to be located close to LTR7 elements in the human genome \[25,
26\], an evolutionary young family of transposable elements that often
harbor binding sites for pluripotency factors. The LTR7 elements near
these two genes appeared along the ape lineage and are bound by POU5F1
\[27\], so they could have potentially contributed to the rewiring of
the POU5F1 module.

### Lineage-specific divergence

The approach described above characterizes the overall degree of module
conservation - all species within a tree are treated the same. This is
the most suitable approach for finding conserved modules and for finding
the highest signal of divergence irrespective of phylogeny. In case of
the example dataset, most of the diverged modules detected this way are
diverged between cynomolgus macaques and great apes (rather than within
the great apes) simply because longer phylogenetic distances generally
give rise to bigger differences.

In many cases, however, the primary interest lies in divergence on a
particular lineage, e.g. the human lineage. To address this, we define
the human lineage as the branch leading from the most recent common
ancestor (MRCA) of humans and gorillas to the human MRCA, using the
macaque as an outgroup. A minimal requirement for this analysis is that
the human replicates form a monophyletic group; trees that don’t fulfill
this requirement are excluded. The human subtree is then defined as the
subtree that includes all human replicates and the human lineage. To
again correct for the unwanted sources of variation, this subtree length
is contrasted against variation within the lineage, i.e. the human
diversity. This allows us to assess module divergence specifically on
the human lineage. An analogous approach can be applied to find diverged
modules on any lineage of interest where an outgroup can be defined.

#### Finding conserved and diverged modules

The modules that are diverged on the human lineage can be defined using
a (weighted) linear regression between the human subtree length and the
human diversity, implemented in `fitTreeStatsLm` with `focus` set to
“human” and `findConservedDivergedModules`. Similarly to the overall
conservation and divergence, the 95% prediction interval of the
regression line is calculated and the modules that have a higher human
subtree length than the upper bound of the 95% prediction interval are
considered diverged. It is important to note that the modules that fall
below the lower bound of the 95% prediction interval cannot be
interpreted as conserved in this case, because 1) the human-monophyletic
module trees have already been pre-selected, and 2) conservation is by
definition not a lineage-specific feature.

Below, this approach is applied to the example dataset:

``` r
lm_human <- fitTreeStatsLm(tree_stats_filt, 
                           focus = "human",
                           weighted_lm = TRUE)
module_conservation_human <- findConservedDivergedModules(tree_stats_filt, 
                                                          lm_human)
```

| focus | regulator | module\_ size | human\_ subtree_length | lwr_human_subtree\_ length | upr_human_subtree\_ length | human\_ diversity |
|:-----:|:---------:|:-------------:|:----------------------:|:--------------------------:|:--------------------------:|:-----------------:|
| human |    AR     |      39       |         0.362          |           0.359            |           0.367            |       0.335       |
| human |  ARID1B   |      43       |         0.426          |           0.423            |           0.434            |       0.338       |
| human |   ARID2   |      48       |         0.330          |           0.324            |           0.335            |       0.304       |
| human |   BNC2    |      49       |         0.319          |           0.316            |           0.324            |       0.286       |
| human |   CREB3   |      49       |         0.352          |           0.349            |           0.354            |       0.313       |
| human |   E2F3    |      46       |         0.390          |           0.385            |           0.395            |       0.355       |

Measures of divergence on the human lineage for six example modules

| lwr_human\_ diversity | upr_human\_ diversity |  fit  | lwr_fit | upr_fit | residual | weight | t_score |  conservation   |
|:---------------------:|:---------------------:|:-----:|:-------:|:-------:|:--------:|:------:|:-------:|:---------------:|
|         0.331         |         0.339         | 0.377 |  0.338  |  0.416  |  -0.015  | 0.010  | -0.786  | not_significant |
|         0.336         |         0.342         | 0.380 |  0.338  |  0.422  |  0.046   | 0.009  |  2.194  |    diverged     |
|         0.300         |         0.308         | 0.343 |  0.311  |  0.376  |  -0.013  | 0.014  | -0.800  | not_significant |
|         0.282         |         0.291         | 0.324 |  0.295  |  0.354  |  -0.005  | 0.018  | -0.346  | not_significant |
|         0.311         |         0.317         | 0.354 |  0.324  |  0.383  |  -0.002  | 0.018  | -0.148  | not_significant |
|         0.352         |         0.360         | 0.398 |  0.366  |  0.429  |  -0.008  | 0.016  | -0.517  | not_significant |

In total, there are 57 modules that fulfill the minimum criteria of
being monophyletic for the human replicates:

``` r
nrow(module_conservation_human)
#> [1] 57
```

Among these, we identified three modules that are diverged on the human
lineage:

``` r
table(module_conservation_human$conservation)
#> 
#>        diverged not_significant 
#>               2              55
```

The results can also be plotted using `plotConservedDivergedModules`,
the same function that has already been used in case of the overall
conservation and divergence:

``` r
plotConservedDivergedModules(module_conservation_human)
```

![Weighted linear regression between the human subtree length and human
diversity, with the diverged modules
marked.](module_conservation_human.png)

Weighted linear regression between the human subtree length and human
diversity, with the diverged modules marked.

The distance matrices and tree reconstructions of the ZNF552, ARID1B and
MXI1 modules can showcase in more detail what makes them diverged on the
human lineage:

``` r
human_div_modules <- module_conservation_human %>% 
  dplyr::filter(conservation == "diverged") %>% 
  dplyr::pull(regulator) %>% 
  as.character()
```

``` r
plotDistMats(dist[human_div_modules])
```

![The distance matrices of the three modules that were found to be
diverged on the human lineage.](dist_mat_human_div.png)

The distance matrices of the three modules that were found to be
diverged on the human lineage.

``` r
plotTrees(trees[human_div_modules], species_colors = spec_colors)
```

![The tree reconstructions of the three modules that were found to be
diverged on the human lineage.](trees_human_div.png)

The tree reconstructions of the three modules that were found to be
diverged on the human lineage.

The signal — unsurprisingly — is not as strong as in case of the modules
diverged between cynomolgus and great apes, but it is still apparent
that the human-gorilla and human-cynomolgus distance measures are higher
than the rest. This also translates into a way longer human subtree
length than expected based on the phylogeny and observed in most module
trees.

#### Finding conserved and diverged target genes

Similar to the approach in case of the overall conservation and
divergence, jackknifing can inform us about target genes that contribute
the most to human-specific divergence. In this case, the residuals of
the jackknifed module versions are calculated with regard to the
regression model between the human subtree length and the human
diversity. The target genes that result in the lowest residuals (i.e.,
the weakest signal of divergence) when removed are considered the most
diverged in the original module.

We performed this analysis for the ARID1B module:

``` r
target_contributions_ARID1B <- findConservedDivergedTargets("ARID1B", 
                                                           tree_stats_jk, 
                                                           lm_human)
```

| focus | regulator | type |         id         | gene\_ removed | human\_ diversity |
|:-----:|:---------:|:----:|:------------------:|:--------------:|:-----------------:|
| human |  ARID1B   | orig |    ARID1B_orig     |       NA       |       0.336       |
| human |  ARID1B   |  jk  |  ARID1B_jk_ACAT2   |     ACAT2      |       0.332       |
| human |  ARID1B   |  jk  | ARID1B_jk_C12orf45 |    C12orf45    |       0.336       |
| human |  ARID1B   |  jk  |  ARID1B_jk_CHCHD2  |     CHCHD2     |       0.305       |

Cross-species conservation measures for a few example target genes of
ARID1B

| human\_ subtree_length |  fit  | lwr_fit | upr_fit | residual | contribution |
|:----------------------:|:-----:|:-------:|:-------:|:--------:|:------------:|
|         0.428          | 0.378 |  0.335  |  0.420  |  0.050   |      NA      |
|         0.423          | 0.373 |  0.331  |  0.415  |  0.050   |    0.000     |
|         0.422          | 0.378 |  0.336  |  0.421  |  0.043   |    0.130     |
|         0.410          | 0.344 |  0.302  |  0.387  |  0.066   |    -0.312    |

``` r
plotConservedDivergedTargets(target_contributions_ARID1B)
```

![The top 5 most diverged target genes of ARID1B identified using
tree-based statistics in combination with
jackknifing.](target_contributions_ARID1B.png)

The top 5 most diverged target genes of ARID1B identified using
tree-based statistics in combination with jackknifing.

``` r
top5_diverged_ARID1B_targets <- target_contributions_ARID1B %>%
  slice_min(order_by = residual, n = 5) %>%
  pull(gene_removed)
```

![Contribution scores of the target genes in the three modules that were
found to be diverged on the human lineage. The two top-contributing
targets are labeled for each module.](target_contributions_human.png)

Contribution scores of the target genes in the three modules that were
found to be diverged on the human lineage. The two top-contributing
targets are labeled for each module.

``` r
plotExprAlongPseudotime(c("ARID1B", top5_diverged_ARID1B_targets), 
                        sce, 
                        species_colors = spec_colors, 
                        cell_type_colors = ct_colors)
```

![The expression profiles of ARID1B and its five most diverged target
genes across the pseudotime trajectory per
species.](target_expression_ARID1B.png)

The expression profiles of ARID1B and its five most diverged target
genes across the pseudotime trajectory per species.

Based on the effects of jackknifing, DLK1 proves to be by far the most
diverged target gene of ARID1B. This is also reflected in its expression
profile: it is upregulated much more in humans during the intermediate
stage of the neural differentiation than in gorilla or cynomolgus.

Some other target genes of ARID1B also show slight human-specific
differences in their expression patterns, e.g., TMSB15A and DTX4 are
both upregulated the most in humans towards the late pseudotime stage.

  
  

## Summary

CroCoNet provides a robust framework to quantitatively compare gene
regulatory networks across species, while accounting for unwanted
biological variation and technical noise. It defines joint co-expression
modules across all species based on a phylogeny-aware consensus network
and quantifies the conservation of these modules by contrasting
cross-species and within-species topological differences. In addition,
it can pinpoint individual target genes that contribute the most to
conservation/divergence with the help of jackknifing. The modules and
target genes detected as conserved or diverged are good candidates to
explain cross-species differences and similarities in gene regulation
and can be investigated in greater depth by follow-up analyses and
experiments.

  
  

## References

1\. Shumate A, Salzberg SL. Liftoff: accurate mapping of gene
annotations. Bioinformatics. 2021;37:1639–43.
<https://doi.org/10.1093/bioinformatics/btaa1016>.

2\. Rosen Y, Brbić M, Roohani Y, Swanson K, Li Z, Leskovec J. Toward
universal cell embeddings: integrating single-cell RNA-seq datasets
across species with SATURN. Nat Methods. 2024.
<https://doi.org/10.1038/s41592-024-02191-z>.

3\. Lun ATL, McCarthy DJ, Marioni JC. A step-by-step workflow for
low-level analysis of single-cell RNA-seq data with Bioconductor.
F1000Res. 2016;5:2122. <https://doi.org/10.12688/f1000research.9501.2>.

4\. Cannoodt R, Saelens W, Sichien D, Tavernier S, Janssens S, Guilliams
M, et al. SCORPIUS improves trajectory inference and identifies novel
modules in dendritic cell development. bioRxiv. 2016;079509.
<https://doi.org/10.1101/079509>.

5\. Aran D, Looney AP, Liu L, Wu E, Fong V, Hsu A, et al.
Reference-based analysis of lung single-cell sequencing reveals a
transitional profibrotic macrophage. Nat Immunol. 2019;20:163–72.
<https://doi.org/10.1038/s41590-018-0276-y>.

6\. Rhodes K, Barr KA, Popp JM, Strober BJ, Battle A, Gilad Y. Human
embryoid bodies as a novel system for genomic studies of functionally
diverse cell types. Elife. 2022;11.
<https://doi.org/10.7554/eLife.71361>.

7\. Moerman T, Santos SA, González-Blas CB, Simm J, Moreau Y, Aerts J,
et al. GRNBoost2 and Arboreto: efficient and scalable inference of gene
regulatory networks. Bioinformatics. 2019;35:2159–61.
<https://doi.org/10.1093/BIOINFORMATICS/BTY916>.

8\. Csardi G, Nepusz T. [The igraph software package for complex network
research](https://igraph.org). InterJournal. 2006;Complex Systems:1695.

9\. Agrawal H. Extreme Self-Organization in Networks Constructed from
Gene Expression Data. Phys Rev Lett. 2002;89:268702.
<https://doi.org/10.1103/PHYSREVLETT.89.268702/FIGURES/3/MEDIUM>.

10\. Bergmann S, Ihmels J, Barkai N. Similarities and Differences in
Genome-Wide Expression Data of Six Organisms. PLoS Biol. 2003;2:e9.
<https://doi.org/10.1371/JOURNAL.PBIO.0020009>.

11\. Barabási AL, Oltvai ZN. Network biology: understanding the cell’s
functional organization. Nature Reviews Genetics 2004 5:2.
2004;5:101–13. <https://doi.org/10.1038/nrg1272>.

12\. Broido AD, Clauset A. Scale-free networks are rare. Nature
Communications 2019 10:1. 2019;10:1–10.
<https://doi.org/10.1038/s41467-019-08746-5>.

13\. Bininda-Emonds ORP, Cardillo M, Jones KE, MacPhee RDE, Beck RMD,
Grenyer R, et al. The delayed rise of present-day mammals. Nature.
2007;446:507–12. <https://doi.org/10.1038/nature05634>.

14\. Rauluseviciute I, Riudavets-Puig R, Blanc-Mathieu R,
Castro-Mondragon JA, Ferenc K, Kumar V, et al. JASPAR 2024: 20th
anniversary of the open-access database of transcription factor binding
profiles. Nucleic Acids Res. 2023;52:D174–82.
<https://doi.org/10.1093/nar/gkad1059>.

15\. Madsen JGS, Rauch A, Van Hauwaert EL, Schmidt SF, Winnefeld M,
Mandrup S. Integrated analysis of motif activity and gene expression
changes of transcription factors. Genome Res. 2018;28:243–55.
<https://doi.org/10.1101/gr.227231.117>.

16\. Langfelder P, Horvath S. WGCNA: An R package for weighted
correlation network analysis. BMC Bioinformatics. 2008;9:1–13.
<https://doi.org/10.1186/1471-2105-9-559>.

17\. Ashburner M, Ball CA, Blake JA, Botstein D, Butler H, Cherry JM, et
al. Gene ontology: tool for the unification of biology. The Gene
Ontology Consortium. Nat Genet. 2000;25:25–9.
<https://doi.org/10.1038/75556>.

18\. Gene Ontology Consortium, Aleksander SA, Balhoff J, Carbon S,
Cherry JM, Drabkin HJ, et al. The Gene Ontology knowledgebase in 2023.
Genetics. 2023;224:iyad031. <https://doi.org/10.1093/genetics/iyad031>.

19\. Alexa A, Rahnenfuhrer J. topGO: enrichment analysis for gene
ontology. R package version. 2024;2.58.0.

20\. Milacic M, Beavers D, Conley P, Gong C, Gillespie M, Griss J, et
al. The Reactome Pathway Knowledgebase 2024. Nucleic Acids Res.
2024;52:D672–8. <https://doi.org/10.1093/nar/gkad1025>.

21\. Yu G, He Q-Y. ReactomePA: an R/Bioconductor package for reactome
pathway analysis and visualization. Mol Biosyst. 2016;12:477–9.
<https://doi.org/10.1039/c5mb00663e>.

22\. Szklarczyk D, Kirsch R, Koutrouli M, Nastou K, Mehryary F, Hachilif
R, et al. The STRING database in 2023: protein-protein association
networks and functional enrichment analyses for any sequenced genome of
interest. Nucleic Acids Res. 2023;51:D638–46.
<https://doi.org/10.1093/nar/gkac1000>.

23\. Yu G, Wang L-G, Han Y, He Q-Y. clusterProfiler: an R package for
comparing biological themes among gene clusters. OMICS. 2012;16:284–7.
<https://doi.org/10.1089/omi.2011.0118>.

24\. Langfelder P, Luo R, Oldham MC, Horvath S. Is my network module
preserved and reproducible? PLoS Comput Biol. 2011;7:1001057.
<https://doi.org/10.1371/journal.pcbi.1001057>.

25\. Wang J, Xie G, Singh M, Ghanbarian AT, Raskó T, Szvetnik A, et al.
Primate-specific endogenous retrovirus-driven transcription defines
naive-like stem cells. Nature. 2014;516:405–9.
<https://doi.org/10.1038/nature13804>.

26\. Hsieh F-K, Ji F, Damle M, Sadreyev RI, Kingston RE. HERVH-derived
lncRNAs negatively regulate chromatin targeting and remodeling mediated
by CHD7. Life Sci Alliance. 2022;5:e202101127.
<https://doi.org/10.26508/lsa.202101127>.

27\. Ito J, Sugimoto R, Nakaoka H, Yamada S, Kimura T, Hayano T, et al.
Systematic identification and characterization of regulatory elements
derived from human endogenous retroviruses. PLoS Genet.
2017;13:e1006883. <https://doi.org/10.1371/journal.pgen.1006883>.
