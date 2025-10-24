#' Transcriptional regulators
#'
#' 7 important transcriptional regulators during early neural differentiation that are used as cores to assemble modules around them.
#'
#' @format A character vector of 7 elements.
"regulators"


#' All genes
#'
#' The 300 genes that are used for the network inference.
#'
#' @format A character vector of 300 elements.
"genes"


#' JASPAR 2024 vertebrate core transcriptional regulators
#'
#' Transcriptional regulators that have at least 1 annotated motif in the JASPAR 2024 vertebrate core collection.
#'
#' @format A character vector of 784 elements.
"jaspar_core_TRs"


#' JASPAR 2024 unvalidated transcriptional regulators
#'
#' Transcriptional regulators that have at least 1 annotated motif in the JASPAR 2024 unvalidated collection.
#'
#' @format A character vector of 590 elements.
"jaspar_unvalidated_TRs"


#' IMAGE transcriptional regulators
#'
#' Transcriptional regulators that have at least 1 annotated motif in the IMAGE database (Madsen et al. 2018).
#'
#' @format A character vector of 1351 elements.
"image_TRs"


#' Replicate-species conversion
#'
#' A data frame that specifies which replicate belongs to which species.
#'
#' @format A data frame with 9 rows and 2 columns:
#' \describe{
#' \item{replicate}{Name of the replicate/cell line.}
#' \item{species}{Name of the species.}
#' }
"replicate2species"


#' SCE object of the primate neural differentiation dataset
#'
#' A subset of the primate neural differentiation scRNA-seq dataset in an SCE format. The data was collected during the early neural differentiation of human, gorilla and cynomolgus macaque iPS cells with 3 human, 2 gorilla and 4 cynomolgus cell lines (replicates). Using a directed differentiation protocol, cells were differentiated into neural progenitor cells (NPCs) over the course of 9 days, and scRNA-seq data was obtained at six time points (days 0, 1, 3, 5, 7 and 9) during this process. The SCE object contains the raw and log-normalized counts as well as the metadata for 300 genes and 900 cells (100 cells per replicate).
#'
#' @format An SCE object with 300 rows, 900 columns, 9 metadata columns and 2 assays.
#'
#' Metadata columns:
#' \describe{
#' \item{species}{Name of the species.}
#' \item{replicate}{Name of the replicate/cell line.}
#' \item{day}{The day when the cell was collected.}
#' \item{n_UMIs}{Number of UMIs detected.}
#' \item{n_genes}{Number of genes detected.}
#' \item{perc_mito}{Percent of mitochondrial reads.}
#' \item{sizeFactor}{Size factor for scaling normalization, calculated first per replicate using [scran::computeSumFactors] and [scran::quickCluster], then adjusted by [batchelor::multiBatchNorm] to remove systematic differences in covergae across replicates.}
#' \item{pseudotime}{Pseudotime inferred by [SCORPIUS::infer_trajectory].}
#' \item{cell_type}{Cell type labels predicted by [SingleR::classifySingleR] using the embryoid body dataset from Rhodes et al. 2022 as reference.}
#' }
#' Assays:
#' \describe{
#' \item{counts}{Raw counts.}
#' \item{logcounts}{Log-normalized counts created by first calculating size factors per replicate using [scran::computeSumFactors] and [scran::quickCluster], then adjusted them to remove systematic differences in covergae across replicates and log-normalizing by [batchelor::multiBatchNorm].}
#' }
"sce"


#' List of raw networks
#'
#' List of networks per replicate with raw edge weights. The networks were inferred using GRNBoost2 based on a subset of the primate neural differentiation scRNA-seq dataset. All 300 genes in the subsetted data were used as potential regulators. To circumvent the stochastic nature of the algorithm, GRNBoost2 was run 10 times on the same count matrices, then the results were averaged across runs, and rarely occurring edges were removed altogether. In addition, edges inferred between the same gene pair but in opposite directions were also averaged.
#'
#' @format A named list of 9 [igraph] objects. Each network contains 300 nodes, and has 1 node and 2 edge attributes:
#'
#' Node attributes:
#' \describe{
#' \item{name}{Name of the node (gene).}
#' }
#' Edge attributes:
#' \describe{
#' \item{weight}{Edge weight, the importance score calculate by GRNBoost2.}
#' }
"network_list_raw"


#' List of genomic annotations
#'
#' List of genomic annotations per species as GRanges objects. The human annotation is the GTF of Hg38 GENCODE release 32 primary assembly, while the gorilla and cynomolgus macaque annotations were created by transferring the human annotation onto the gorGor6 and macFas6 genomes via the tool Liftoff (https://github.com/agshumate/Liftoff). Each annotation was subsetted for the 300 genes that feature in this example dataset.
#'
#' @format A named list of 3 GRanges objects. Each object has 8 columns.
#'
#' Attributes:
#' \describe{
#' \item{seqnames}{Chromosome or contig name.}
#' \item{start}{Genomic start location.}
#' \item{end}{Genomic end location.}
#' \item{width}{Width of the feature in base pairs.}
#' \item{strand}{Genomic strand ("+" or "-").}
#' \item{source}{The prediction program or public database where the annotations came from.}
#' \item{type}{Feature type (gene, transcript, exon, CDS, UTR, start_codon, stop_codon or Selenocysteine).}
#' \item{score}{The degree of confidence in the feature's existence and coordinates.}
#' \item{phase}{One of '0', '1' or '2'. '0' means that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on.}
#' \item{gene_id}{Unique identifier of the gene.}
#' \item{gene_name}{Name of the gene.}
#' \item{transcript_id}{Unique identifier of the transcript.}
#' \item{transcript_name}{Name of the transcript.}
#' }
#' @source <ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz>, <https://hgdownload.soe.ucsc.edu/goldenPath/gorGor6/bigZips/gorGor6.fa.gz>, <https://ftp.ensembl.org/pub/release-109/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna_sm.toplevel.fa.gz>
"gtf_list"


#' List of networks
#'
#' List of networks per replicate with edge weights re-scaled between 0 and 1, and gene pairs (edges) with overlapping annotations removed. The networks were inferred using GRNBoost2 based on a subset of the primate neural differentiation scRNA-seq dataset. All 300 genes in the subsetted data were used as potential regulators. To circumvent the stochastic nature of the algorithm, GRNBoost2 was run 10 times on the same count matrices, then the results were averaged across runs, and rarely occurring edges were removed altogether. In addition, edges inferred between the same gene pair but in opposite directions were also averaged. Edge weights were scaled by the maximum edge weight across all replicates. Gene pairs that have overlapping annotations in any of the species' genomes were removed from all networks.
#'
#' @format A named list of 9 [igraph] objects. Each network contains 300 nodes, and has 1 node attribute and 3 edge attributes:
#'
#' Node attributes:
#' \describe{
#' \item{name}{Name of the node (gene).}
#' }
#' Edge attributes:
#' \describe{
#' \item{weight}{Edge weight, the importance score calculate by GRNBoost2 rescaled between 0 and 1.}
#' \item{genomic_dist}{Numeric, the genomic distance of the 2 genes that form the edge (Inf if the 2 genes are annotated on different chromosomes/contigs).}
#' }
"network_list"


#' Phylogenetic tree
#'
#' Rooted phylogenetic tree of 3 primate species: human (Homo sapiens), gorilla (Gorilla gorilla) and cynomolgus macaque (Macaca Fascicularis). The tree was created by subsetting the mammalian tree with the best estimates of branch lengths from Bininda-Edmons et al. 2007.
#'
#' @format A [phylo] object with 4 edges and 3 nodes.
#' @source Olaf R. P. Bininda-Emonds, Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman & Andy Purvis. "The delayed rise of present-day mammals" Nature 446, 507-512(29 March 2007). doi:10.1038/nature05634
"tree"


#' Consensus network
#'
#' Consensus network of the 9 primate replicates in the example dataset. For each edge, the consensus adjacency was calculated as the weighted average of replicate-wise adjacencies using weights that correct for 1) the phylogenetic distances between species and 2) the different numbers of replicates per species. If an edge was not detected in certain replicate, the adjacency of that replicate was regarded as 0 for the calculation of the consensus. The directionality of each edge was determined based on a modified Spearman's correlation between the corresponding 2 genes' expression profiles (positive expression correlation - activating interaction, negative expression correlation - repressing interaction). The correlations were calculated per replicate, then the mean correlation was taken across all replicates.
#'
#' @format An [igraph] object with 300 nodes, and has 1 node attribute and 3 edge attributes:
#'
#' Node attributes:
#' \describe{
#' \item{name}{Name of the node (gene).}
#' }
#' Edge attributes:
#' \describe{
#' \item{weight}{Consensus edge weight/adjacency, the weighted average of replicate-wise adjacencies.}
#' \item{rho}{Approximate Spearman's correlation coefficient of the 2 genes' expression profiles that form the edge.}
#' \item{p.adj}{BH-corrected approximate p-value of rho.}
#' \item{direction}{Direction of the interaction between the 2 genes that form the edge ("+" or "-").}
#' }
"consensus_network"


#' Initial modules
#'
#' Initial modules created by assigning the top 250 targets to each of 7 transcriptional regulators involved in the early neuronal differentiation of primates (\code{regulators}). For each regulator, the top 250 targets were selected based edge weights in the consensus network: all targets of the regulator were ranked based on their edge weight to the regulator (regulator-taregt adjacency) and the 250 targets with the highest regulator-target adjacencies were kept.
#'
#' @format A data frame with 1750 rows and 8 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Target gene of the transcriptional regulator (member of the regulator's initial module).}
#' \item{weight}{Consensus edge weight/adjacency, the weighted average of replicate-wise adjacencies.}
#' \item{rho}{Approximate Spearman's correlation coefficient of the 2 genes' expression profiles that form the edge.}
#' \item{p.adj}{BH-corrected approximate p-value of rho.}
#' \item{direction}{Direction of the interaction between the 2 genes that form the edge ("+" or "-").}
#' }
"initial_modules"


#' Pruned modules
#'
#' Pruned modules created by the dynamic pruning of the initial modules via the "UIK_adj_kIM" method. The initial module members were filtered in successive steps based on their adjacency to the regulator and their intramodular connectivitiy alternately. In each step, the cumulative sum curve based on one of these two characteristics was calculated per module, then the targets below the knee point of the curve were kept. This process was continued until the module sizes became as small as possible without the median falling below a pre-defined minimum of 20 genes.
#'
#' @format A data frame with 225 rows and 9 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted average of replicate-wise adjacencies.}
#' \item{rho}{Approximate Spearman's correlation coefficient of the 2 genes' expression profiles that form the edge.}
#' \item{p.adj}{BH-corrected approximate p-value of rho.}
#' \item{direction}{Character specifying the direction of regulation between the regulator and the target, either "+" or "-".}
#' \item{module_size}{Module size, the numer of target genes assigned to a regulator.}
#' }
"pruned_modules"


#' Random modules
#'
#' Random modules that match the actual (pruned) modules in size. These random modules have the same regulators and contain the same number of target genes as the actual modules, but the target genes were randomly drawn from all genes in the network.
#'
#' @format A data frame with 225 rows and 3 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Member of the regulator's random module.}
#' \item{module_size}{Module size, the numer of target genes assigned to a regulator.}
#' }
"random_modules"


#' Eigengenes
#'
#' Eigengenes of the pruned modules calculated using the activated targets in each module. An eigengene summarizes the expression profile of the module as a whole, mathematically it is the first principal component of the module expression data (i.e. the scaled and centered logcounts subsetted for the activated targets and the regulators of the given module). In this case, the principal component was taken across all cells irrespective of species.
#'
#' @format A data frame with 6300 rows and 8 columns:
#'\describe{
#' \item{cell}{Character, the cell barcode.}
#' \item{species}{Character, the name of the species.}
#' \item{pseudotime}{Numeric, inferred pseudotime.}
#' \item{cell_type}{Character, cell type annotation.}
#' \item{module}{Character, transcriptional regulator and direction of regulation (in this case always nameOfRegulator(+)).}
#' \item{eigengene}{Numeric, the eigengene (i.e. the first principal component of the scaled and centered logcounts) of the module.}
#' \item{mean_expr}{Numeric, the mean of the scaled and centered logcounts across all genes in the module.}
#' \item{regulator_expr}{Numeric, the scaled and centered logcounts of the regulator.}
#' }
"eigengenes"


#' Eigengenes per species
#'
#' Eigengenes of the pruned modules calculated per species using the activated targets in each module. An eigengene summarizes the expression profile of the module as a whole, mathematically it is the first principal component of the module expression data (i.e. the scaled and centered logcounts subsetted for the activated targets and the regulator of the given module). In this case, the principal component was calculated for each species separately.
#'
#' @format A data frame with 6300 rows and 8 columns:
#'\describe{
#' \item{cell}{Character, the cell barcode.}
#' \item{species}{Character, the name of the species.}
#' \item{pseudotime}{Numeric, inferred pseudotime.}
#' \item{cell_type}{Character, cell type annotation.}
#' \item{module}{Character, transcriptional regulator and direction of regulation (in this case always nameOfRegulator(+)).}
#' \item{eigengene}{Numeric, the eigengene (i.e. the first principal component of the scaled and centered logcounts) of the module.}
#' \item{mean_expr}{Numeric, the mean of the scaled and centered logcounts across all genes in the module.}
#' \item{regulator_expr}{Numeric, the scaled and centered logcounts of the regulator.}
#' }
"eigengenes_per_species"


#' Preservation statistics of the original and jackknifed pruned modules
#'
#' Correlation of intramodular connectivities (cor_kIM) per replicate pair for the original and all jackknifed versions of the pruned modules. The jackknifed versions of the modules were created by removing each target gene assigned to a module (the regulators were never excluded). The preservation statistic cor.kIM was then calculated for the original module as well as each jackknife module version by comparing each replicate to all others, both within and across species. cor.kIM quantifies how well the connectivity patterns are preserved between the networks of two replicates, mathematically it is the correlation of the intramodular connectivities per module member gene in the network of the 1st replicate VS the intramodular connectivities per module member gene in the network of the 2nd replicate.
#'
#' @format A data frame with 8352 rows and 8 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type (orig = original or jk = jackknifed).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig').}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type 'orig').}
#' \item{replicate1, replicate2}{Character, the names of the replicates compared.}
#' \item{species1, species2}{Character, the names of the species \code{replicate1} and \code{replicate2} belong to, respectively.}
#' \item{cor_adj}{Numeric, correlation of adjacencies.}
#' \item{cor_kIM}{Numeric, correlation of intramodular connectivities.}
#' }
"pres_stats_jk"


#' Preservation statistics of the the original and jackknifed random modules
#'
#' Correlation of intramodular connectivities (cor_kIM) per replicate pair for the original and all jackknifed versions of the random modules. The jackknifed versions of the modules were created by removing each target gene assigned to a module (the regulators were never excluded). The preservation statistic cor.kIM was then calculated for the original module as well as each jackknife module version by comparing each replicate to all others, both within and across species. cor.kIM quantifies how well the connectivity patterns are preserved between the networks of two replicates, mathematically it is the correlation of the intramodular connectivities per module member gene in the network of the 1st replicate VS the intramodular connectivities per module member gene in the network of the 2nd replicate.
#'
#' @format A data frame with 8352 rows and 8 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type (orig = original or jk = jackknifed).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig').}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type 'orig').}
#' \item{replicate1, replicate2}{Character the names of the replicates compared.}
#' \item{species1, species2}{Character, the names of the species \code{replicate1} and \code{replicate2} belong to, respectively.}
#' \item{cor_adj}{Numeric, correlation of adjacencies.}
#' \item{cor_kIM}{Numeric, correlation of intramodular connectivities.}
#' }
"random_pres_stats_jk"


#' Distance measures of the original and jackknifed pruned modules
#'
#' Distance measures per replicate pair for the original and all jackknifed versions of the pruned modules. The jackknifed versions of the modules were created by removing each target gene assigned to a module (the regulators were never excluded). The distance measures were calculated based on the correlation of intramodular connectivities: \eqn{dist = \frac{1 - cor.kIM}{2}} (a correlation of 1 corresponds to a distance of 0, whereas a correlation of -1 corresponds to a distance of 1). Each element of the list corresponds a jackknifed/original module version and contains the distance measures between all possible pairs of replicates for this module version.
#'
#' @format A named list with 232 elements containing the distance measures per (original or jackknifed) module version. Each element is a data frame with 72 rows and 10 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type (orig = original or jk = jackknifed).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig').}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type 'orig').}
#' \item{replicate1, replicate2}{Character the names of the replicates compared.}
#' \item{species1, species2}{Character, the names of the species \code{replicate1} and \code{replicate2} belong to, respectively.}
#' \item{dist}{Numeric, distance measure ranging from 0 to 1, calculated based on the correlation of intramodular connectivities.}
#' }
"dist_jk"


#' Distance measures of the pruned modules
#'
#' Distance measures per replicate pair for the original (i.e. not jackknifed) pruned modules.
#'
#' The distance measures were calculated based on the correlation of intramodular connectivities: \deqn{dist = \frac{1 - cor.kIM}{2}} (a correlation of 1 corresponds to a distance of 0, whereas a correlation of -1 corresponds to a distance of 1). Each element of the list corresponds to a module and contains the distance measures between all possible pairs of replicates for this module.
#'
#' @format A named list with 12 elements containing the distance measures per module. Each element is a data frame with 72 rows and 7 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{replicate1, replicate2}{Character the names of the replicates compared.}
#' \item{species1, species2}{Character, the names of the species \code{replicate1} and \code{replicate2} belong to, respectively.}
#' \item{dist}{Numeric, distance measure ranging from 0 to 1, calculated based on the correlation of intramodular connectivities.}
#' }
"dist"


#' Distance measures of the the original and jackknifed random modules
#'
#' Distance measures per replicate pair for the original and all jackknifed versions of the random modules.
#'
#' The jackknifed versions of the modules were created by removing each target gene assigned to a module (the regulators were never excluded). The distance measures were calculated based on the correlation of intramodular connectivities: \deqn{dist = \frac{1 - cor.kIM}{2}} (a correlation of 1 corresponds to a distance of 0, whereas a correlation of -1 corresponds to a distance of 1). Each element of the list corresponds a jackknife module version and contains the distance measures between all possible pairs of replicates for this module version.
#'
#' @format A named list with 232 elements containing the distance measures per (original or jackknifed) module version. Each element is a data frame with 72 rows and 10 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type (orig = original or jk = jackknifed).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig').}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type 'orig').}
#' \item{replicate1, replicate2}{Character the names of the replicates compared.}
#' \item{species1, species2}{Character, the names of the species \code{replicate1} and \code{replicate2} belong to, respectively.}
#' \item{dist}{Numeric, distance measure ranging from 0 to 1, calculated based on the correlation of intramodular connectivities.}
#' }
"random_dist_jk"


#' Trees of the original and jackknifed pruned modules
#'
#' Neighbor-joining trees representing the similarities of connectivity patterns across the 9 primate replicates for the original and all jackknifed versions of the pruned modules. The jackknifed versions of the modules were created by removing each target gene assigned to a module (the regulators were never excluded). For each of these module versions, the trees were inferred based on the preservation statistic cor.kIM (correlation of intramodular connectivities): first, the preservation scores were calculated between all possible replicate pairs, then they were converted into a distance matrix of replicates, and finally trees were reconstructed based on this distance matrix using the neighbor-joining algorithm. The result is a single tree for each original or jackknife module where the tips represent the replicates and the branch lengths represent the dissimilarity of connectivity patterns between these replicates.
#' @format A named list with 232 elements containing the neighbor-joining trees as [phylo] objects.
"trees_jk"


#' Trees of the pruned modules
#'
#' Neighbor-joining trees representing the similarities of connectivity patterns across the 9 primate replicates for the original (i.e. not jackknifed) pruned modules. For each of the modules, the trees were inferred based on the preservation statistic cor.kIM (correlation of intramodular connectivities): first, the preservation scores were calculated between all possible replicate pairs, then they were converted into a distance matrix of replicates, and finally trees were reconstructed based on this distance matrix using the neighbor-joining algorithm. The result is a single tree per module where the tips represent the replicates and the branch lengths represent the dissimilarity of connectivity patterns between these replicates.
#' @format A named list with 12 elements containing the neighbor-joining trees as [phylo] objects.
"trees"


#' Trees of the original and jackknifed random modules
#'
#' Neighbor-joining trees representing the similarities of connectivity patterns across the 9 primate replicates for the original and all jackknifed versions of the random modules. The jackknifed versions of the modules were created by removing each target gene assigned to a module (the regulators were never excluded). For each of these module versions, the trees were inferred based on the preservation statistic cor.kIM (correlation of intramodular connectivities): first, the preservation scores were calculated between all possible replicate pairs, then they were converted into a distance matrix of replicates, and finally trees were reconstructed based on this distance matrix using the neighbor-joining algorithm. The result is a single tree for each original or jackknife module where the tips represent the replicates and the branch lengths represent the dissimilarity of connectivity patterns between these replicates.
#' @format A named list with 232 elements containing the neighbor-joining trees as [phylo] objects.
"random_trees_jk"


#' Tree-based statistics of the original and jackknifed pruned modules
#'
#' Tree-based statistics calculated for the original and all jackknifed versions of the pruned modules. For each module version, a tree was reconstructed based on the module preservation scores within and across species. The tips of resulting tree represent the replicates and the branch lengths represent the dissimilarity of connectivity patterns between the replicates. Branches between replicates of different species carry information about cross-species differences, while the branches between replicates of the same species carry information about the within-species diversity.
#' @format A data frame with 232 rows and 16 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Module size, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type (orig = original or jk = jackknifed).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig').}
#' \item{gene_removed}{The name of the gene removed by jackknifing (NA in case of module type 'orig').}
#' \item{within_species_diversity}{Numeric, the sum of the human, gorilla and cynomolgus diversity.}
#' \item{total_tree_length}{Numeric, the sum of all branch lengths in the tree.}
#' \item{human_diversity}{Numeric, the sum of the branch lengths in the subtree that contains only the human tips.}
#' \item{gorilla_diversity}{Numeric, the sum of the branch lengths in the subtree that contains only the gorilla tips.}
#' \item{cynomolgus_diversity}{Numeric, the sum of the branch lengths in the subtree that contains only the cynomolgus tips.}
#' \item{human_monophyl}{Logical indicating whether the tree is monophyletic for the human replicates.}
#' \item{gorilla_monophyl}{Logical indicating whether the tree is monophyletic for the gorilla replicates.}
#' \item{cynomolgus_monophyl}{Logical indicating whether the tree is monophyletic for the cynomolgus replicates.}
#' \item{human_subtree_length}{Numeric, the sum of the branch lengths in the subtree that is defined by the human replicates and includes the internal branch connecting the human replicates to the rest of the tree. NA if the tree is not monophyletic for the human replicates.}
#' \item{gorilla_subtree_length}{Numeric, the sum of the branch lengths in the subtree that is defined by the gorilla replicates and includes the internal branch connecting the gorilla replicates to the rest of the tree. NA if the tree is not monophyletic for the gorilla replicates.}
#' \item{cynomolgus_subtree_length}{Numeric, the sum of the branch lengths in the subtree that is defined by the cynomolgus replicates and includes the internal branch connecting the cynomolgus replicates to the rest of the tree. NA if the tree is not monophyletic for the cynomolgus replicates.}
#' }
"tree_stats_jk"


#' Tree-based statistics of the original and jackknifed random modules
#'
#' Tree-based statistics calculated for the original and all jackknifed versions of the random modules. For each module version, a tree was reconstructed based on the module preservation scores within and across species. The tips of resulting tree represent the replicates and the branch lengths represent the dissimilarity of connectivity patterns between the replicates. Branches between replicates of different species carry information about cross-species differences, while the branches between replicates of the same species carry information about the within-species diversity.
#' @format A data frame with 232 rows and 16 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Module size, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type (orig = original or jk = jackknifed).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig').}
#' \item{gene_removed}{The name of the gene removed by jackknifing (NA in case of module type 'orig').}
#' \item{within_species_diversity}{Numeric, the sum of the human, gorilla and cynomolgus diversity.}
#' \item{total_tree_length}{Numeric, the sum of all branch lengths in the tree.}
#' \item{human_diversity}{Numeric, the sum of the branch lengths in the subtree that contains only the human tips.}
#' \item{gorilla_diversity}{Numeric, the sum of the branch lengths in the subtree that contains only the gorilla tips.}
#' \item{cynomolgus_diversity}{Numeric, the sum of the branch lengths in the subtree that contains only the cynomolgus tips.}
#' \item{human_monophyl}{Logical indicating whether the tree is monophyletic for the human replicates.}
#' \item{gorilla_monophyl}{Logical indicating whether the tree is monophyletic for the gorilla replicates.}
#' \item{cynomolgus_monophyl}{Logical indicating whether the tree is monophyletic for the cynomolgus replicates.}
#' \item{human_subtree_length}{Numeric, the sum of the branch lengths in the subtree that is defined by the human replicates and includes the internal branch connecting the human replicates to the rest of the tree. NA if the tree is not monophyletic for the human replicates.}
#' \item{gorilla_subtree_length}{Numeric, the sum of the branch lengths in the subtree that is defined by the gorilla replicates and includes the internal branch connecting the gorilla replicates to the rest of the tree. NA if the tree is not monophyletic for the gorilla replicates.}
#' \item{cynomolgus_subtree_length}{Numeric, the sum of the branch lengths in the subtree that is defined by the cynomolgus replicates and includes the internal branch connecting the cynomolgus replicates to the rest of the tree. NA if the tree is not monophyletic for the cynomolgus replicates.}
#' }
"random_tree_stats_jk"


#' Preservation statistics of the pruned modules
#'
#' Correlation of intramodular connectivities (cor_kIM) per replicate pair and module. The preservation statistic cor.kIM quantifies how well the connectivity patterns are preserved between the networks of two replicates, mathematically it is the correlation of the intramodular connectivities per module member gene in the network of the 1st replicate VS the intramodular connectivities per module member gene in the network of the 2nd replicate.
#' This statistic was calculated for all possible jackknifed versions of the modules, each of which was created by removing a target gene assigned to the given module (the regulators were never excluded). Each jackknifed module was compared between all posible pairs of replicates, both within and across species, resulting in a cor.kIM value per jackknifed module version and replicate pair. Finally, the cor.kIM values were summarized per module and replicate pair by taking the median and its 95\% confidence interval across all jackknifed module versions.
#'
#' @format A data frame with 252 rows and 10 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Module size, the numer of target genes assigned to a regulator.}
#' \item{replicate1, replicate2}{The names of the replicates compared.}
#' \item{species1, species2}{The names of the species 'replicate1' and 'replicate2' belongs to, respectively.}
#' \item{cor_kIM}{The median of cor.kIM across all jackknifed versions of the module.}
#' \item{var_cor_kIM}{The variance of cor.kIM across all jackknifed versions of the module.}
#' \item{lwr_cor_kIM}{The lower bound of the 95\% confidence interval of cor.kIM calculated by jackknifing.}
#' \item{upr_cor_kIM}{The upper bound of the 95\% confidence interval of cor.kIM calculated by jackknifing.}
#' \item{cor_adj}{The median of cor.adj across all jackknifed versions of the module.}
#' \item{var_cor_adj}{The variance of cor.adj across all jackknifed versions of the module.}
#' \item{lwr_cor_adj}{The lower bound of the 95\% confidence interval of cor.adj calculated by jackknifing.}
#' \item{upr_cor_adj}{The upper bound of the 95\% confidence interval of cor.adj calculated by jackknifing.}
#' }
"pres_stats"


#' Preservation statistics of the random modules
#'
#' Correlation of intramodular connectivities (cor_kIM) per replicate pair and random module. The preservation statistic cor.kIM quantifies how well the connectivity patterns are preserved between the networks of two replicates, mathematically it is the correlation of the intramodular connectivities per module member gene in the network of the 1st replicate VS the intramodular connectivities per module member gene in the network of the 2nd replicate.
#' This statistic was calculated for all possible jackknifed versions of the modules, each of which was created by removing a target gene assigned to the given module (the regulators were never excluded). Each jackknifed module was compared between all posible pairs of replicates, both within and across species, resulting in a cor.kIM value per jackknifed module version and replicate pair. Finally, the cor.kIM values were summarized per module and replicate pair by taking the median and its 95\% confidence interval across all jackknifed module versions.
#'
#' @format A data frame with 252 rows and 10 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Module size, the numer of target genes assigned to a regulator.}
#' \item{replicate1, replicate2}{The names of the replicates compared.}
#' \item{species1, species2}{The names of the species 'replicate1' and 'replicate2' belongs to, respectively.}
#' \item{cor_kIM}{The median of cor.kIM across all jackknifed versions of the module.}
#' \item{var_cor_kIM}{The variance of cor.kIM across all jackknifed versions of the module.}
#' \item{lwr_cor_kIM}{The lower bound of the 95\% confidence interval of cor.kIM calculated by jackknifing.}
#' \item{upr_cor_kIM}{The upper bound of the 95\% confidence interval of cor.kIM calculated by jackknifing.}
#' \item{cor_adj}{The median of cor.adj across all jackknifed versions of the module.}
#' \item{var_cor_adj}{The variance of cor.adj across all jackknifed versions of the module.}
#' \item{lwr_cor_adj}{The lower bound of the 95\% confidence interval of cor.adj calculated by jackknifing.}
#' \item{upr_cor_adj}{The upper bound of the 95\% confidence interval of cor.adj calculated by jackknifing.}
#' }
"random_pres_stats"


#' Tree-based statistics of the pruned modules
#'
#' Tree-based statistics per module. The trees were reconstructed based on module preservation scores (cor.kIM) within and across species. The tips of resulting tree represent the replicates and the branch lengths represent the dissimilarity of connectivity patterns between the replicates. The total tree length is the sum of all branch lengths in the tree and it measures module variability both within and across species, while the within-species_diversity is the sum of the human, gorilla and cynomolgus within-species branch lengths and it measures module variability within species only.
#' These statistics were calculated for all possible jackknifed versions of the modules, each of which was created by removing a target gene assigned to the given module (the regulators were never excluded). Then the statistics were summarized per module by taking the median and its 95\% confidence interval across all jackknifed module versions.
#'
#' @format A data frame with 7 rows and 10 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Module size, the numer of target genes assigned to a regulator.}
#' \item{total_tree_length}{The median of total tree lengths across all jackknifed versions of the module.}
#' \item{var_total_tree_length}{The variance of total tree lengths across all jackknifed versions of the module.}
#' \item{lwr_total_tree_length}{The lower bound of the 95\% confidence interval of the total tree length calculated by jackknifing.}
#' \item{upr_total_tree_length}{The upper bound of the 95\% confidence interval of the total tree length calculated by jackknifing.}
#' \item{within_species_diversity}{The median of within-species_diversity values across all jackknifed versions of the module.}
#' \item{var_within_species_diversity}{The variance of within-species_diversity values across all jackknifed versions of the module.}
#' \item{lwr_within_species_diversity}{The lower bound of the 95\% confidence interval of the within-species_diversity calculated by jackknifing.}
#' \item{upr_within_species_diversity}{The upper bound of the 95\% confidence interval of the within-species_diversity calculated by jackknifing.}
#' }
"tree_stats"


#' Tree-based statistics of the random modules
#'
#' Tree-based statistics per random module. The trees were reconstructed based on module preservation scores (cor.kIM) within and across species. The tips of resulting tree represent the replicates and the branch lengths represent the dissimilarity of connectivity patterns between the replicates. The total tree length is the sum of all branch lengths in the tree and it measures module variability both within and across species, while the within-species_diversity is the sum of the human, gorilla and cynomolgus within-species branch lengths and it measures module variability within species only.
#' These statistics were calculated for all possible jackknifed versions of the modules, each of which was created by removing a target gene assigned to the given module (the regulators were never excluded). Then the statistics were summarized per module by taking the median and its 95\% confidence interval across all jackknifed module versions.
#'
#' @format A data frame with 7 rows and 10 columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Module size, the numer of target genes assigned to a regulator.}
#' \item{total_tree_length}{The median of total tree lengths across all jackknifed versions of the module.}
#' \item{var_total_tree_length}{The variance of total tree lengths across all jackknifed versions of the module.}
#' \item{lwr_total_tree_length}{The lower bound of the 95\% confidence interval of the total tree length calculated by jackknifing.}
#' \item{upr_total_tree_length}{The upper bound of the 95\% confidence interval of the total tree length calculated by jackknifing.}
#' \item{within_species_diversity}{The median of within-species_diversity values across all jackknifed versions of the module.}
#' \item{var_within_species_diversity}{The variance of within-species_diversity values across all jackknifed versions of the module.}
#' \item{lwr_within_species_diversity}{The lower bound of the 95\% confidence interval of the within-species_diversity calculated by jackknifing.}
#' \item{upr_within_species_diversity}{The upper bound of the 95\% confidence interval of the within-species_diversity calculated by jackknifing.}
#' }
"random_tree_stats"


#' Linear model for the characterization of overall module conservation
#'
#' A weighted linear model between the total tree length and within-species diversity of the 12 modules in the subsetted early neuronal differentiation dataset. It captures the general relationship between the two tree-based statistics in this dataset and thus also describes the expected total tree length for any given within-species diversity. By identifying modules that differ the most from this expectation, i.e. identifying the outlier data points, it can be used to pinpoint conserved and overall diverged modules. In combination with jackknifing, it can also be used to identify target genes within these modules that contribute the most to the conservation/divergence. The linear model was fit by calling \code{\link{lm}} with weights inversely proportional to the variance of the total tree lengths.
#'
#' @format An object of class \code{\link{lm}} fitted on the object \code{tree_stats} with the formula "total_tree_length ~ within_species_diversity".
"lm_overall"


#' Cross-species conservation measures per module
#'
#' Cross-species conservation measures of the 12 modules in the subsetted early neuronal differentiation dataset. Modules that were found to be conserved or diverged overall (across all species) are labelled as "conserved" or "diverged" in the column \code{conservation}.
#'
#' To determine whether a module as whole is conserved or diverged overall, module trees were reconstructed from pairwise preservation scores between replicates and based on these trees 2 useful statistics were calculated for each module: the total tree length and the within-species diversity (for details please see \code{\link{calculatePresStats}}, \code{\link{reconstructTrees}}, \code{\link{calculateTreeStats}}). After fitting a weighted linear model between the total tree length and within-species diversity values of all modules, a module was considered diverged if it fell above the prediction interval of the regression line, while a module was considered conserved if it fell below the prediction interval of the regression line (for details please see \code{\link{fitTreeStatsLm}} and \code{\link{findConservedDivergedModules}}. The degree of conservation/divergence can be further compared between the modules categorized as conserved/diverged using 2 measures, the residual and the t-score.
#'
#' @format A data frame with 12 rows and 16 columns:
#' \describe{
#' \item{focus}{Character, the focus of interest in terms of cross-species conservation ("overall").}
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{total_tree_length}{Numeric, total tree length per module (the median across all jackknife versions of the module).}
#' \item{lwr_total_tree_length}{Numeric, the lower bound of the confidence interval of the total tree length calculated based on the jackknifed versions of the module.}
#' \item{upr_total_tree_length}{Numeric, the upper bound of the confidence interval of the total tree length calculated based on the jackknifed versions of the module.}
#' \item{within_species_diversity}{Numeric, within-species diveristy per module (the median across all jackknife versions of the module).}
#' \item{lwr_within_species_diversity}{Numeric, the lower bound of the confidence interval of the within-species diversity calculated based on the jackknifed versions of the module.}
#' \item{upr_within_species_diversity}{Numeric, the upper bound of the confidence interval of the within-species diversity calculated based on the jackknifed versions of the module.}
#' \item{fit}{Numeric, the fitted total tree length at the within-species diversity value of the module.}
#' \item{lwr_fit}{Numeric, the lower bound of the prediction interval of the fit.}
#' \item{upr_fit}{Numeric, the upper bound of the prediction interval of the fit.}
#' \item{residual}{Numeric, the residual of the module in the linear model. It is calculated as the difference between the observed and expected (fitted) total tree lengths.}
#' \item{weight}{Numeric, the weight of the module in the linear regression, inversely proportional to the variance of total tree lengths.}
#' \item{t_score}{Numeric, the t-score of the module. It is calculated as the residual normalized by the standard error of the total tree length prediction at the given within-species diversity value.}
#' \item{conservation}{Character, "not_significant" if the module falls inside the prediction interval of the fit, "diverged" if a module has a higher total tree length than the upper boundary of the prediction interval, and "conserved" if a module has a lower total tree length than the lower boundary of the prediction interval.}
#' }
"module_conservation_overall"


#' Cross-species conservation measures of the target genes in the POU5F1 module
#'
#' Cross-species conservation measures of the 21 target genes in the POU5F1 module in the subsetted early neuronal differentiation dataset.
#'
#' To determine which target genes are the most responsible for the divergence of the POU5F1 module, the same statistics were be used as for the identification of conserved/diverged modules but in combination with jackknifing. This means that each target gene in the module was removed and all statistics, including the tree-based statistics that carry information about cross-species conservation, were re-calculated (please see \code{\link{calculatePresStats}}, \code{\link{reconstructTrees}} and \code{\link{calculateTreeStats}}). Our working hypothesis is that if removing a target gene from the diverged POU5F1 module makes the module more conserved, then this gene contributed to the divergence of the original module in the first place. This can be quantified using the residuals compared to the regression line between the total tree lengths and within-species diversities of the tree representations across all modules (please see \code{\link{fitTreeStatsLm}} and \code{\link{findConservedDivergedTargets}}. The jackknifed versions of the POU5F1 module that have the lowest residuals are the most conserved, therefore the corresponding target genes are considered the most diverged.
#'
#' @format A data frame with 23 rows and 12 columns:
#' \describe{
#' \item{focus}{Character, the focus of interest in terms of cross-species conservation ("overall").}
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type (orig = original, jk = jackknifed or summary = the summary across all jackknife module versions).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk', nameOfRegulator_orig in case of module type 'orig' and nameOfRegulator_summary in case of module type 'summary').}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module types 'orig' and 'summary').}
#' \item{total_tree_length}{Numeric, total tree length per jackknife module version.}
#' \item{within_species_diversity}{Numeric, within-species diversity per jackknife module version.}
#' \item{fit}{Numeric, the fitted total tree length at the within-species diversity value of a jackknife module version.}
#' \item{lwr_fit}{Numeric, the lower bound of the prediction interval of the fit.}
#' \item{upr_fit}{Numeric, the upper bound of the prediction interval of the fit.}
#' \item{residual}{Numeric, the residual of the jackknife module version in the linear model. It is calculated as the difference between the observed and expected (fitted) total tree lengths.}
#' \item{contribution}{Numeric, the contribution score of the removed target gene.}
#' }
"POU5F1_target_conservation"
