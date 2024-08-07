library(CroCoNet)
library(foreach)
library(tidyverse)
library(scran)
library(transformGamPoi)
library(ape)
library(SingleCellExperiment)
library(foreach)


library(igraph)

library(DescTools)
library(foreach)
library(igraph)
library(data.table)

library(plyranges)
library(BiocParallel)
library(doParallel)
library(scran)
library(inflection)
library(stats)
library(Matrix)
library(caper)
library(ggtree)
library(tidytree)
library(ggrepel)



## Genes ----------------------------------------------------------

# list of regulators
regulators <- sort(c("SOX2", "POU5F1", "ONECUT2", "TCF4", "BNC2", "TEAD1", "SMAD2", "ZNF701", "SOX15", "CREB3", "PROX1", "FOS"))
usethis::use_data(regulators, compress = "bzip2", overwrite = T)

# genes in the modules of 7 famous TFs (and all other TFs excluded)
genes <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_scripts/toy_data_playaround4/gene9_cell15/RDS/genes.rds")
usethis::use_data(genes, compress = "bzip2", overwrite = T)



## Clone 2 species --------------------------------------------------------

# read and save object
clone2species <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/RDS/clone2species.rds") %>%
  dplyr::mutate(clone = gsub(".i", "", clone))
usethis::use_data(clone2species, compress = "bzip2", overwrite = T)

# clone names
clone_names <- clone2species$clone


## SCE object -----------------------------------------------------------

# randomly select 100 cells per clone
sce <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_scripts/toy_data_playaround4/gene9_cell15/RDS/sce.rds")
counts(sce)["FOS", c( "TAGGCATGGGCAGC", "TAGGCATGATTAAG", "CAGAGAGGCTTTAA", "GGAGCTACTTATAC", "GGAGCTACGTATAA", "TAGCGCTCAGACCT",  "ATCTCAGGTACTAT", "CGGAGCCTCAGGCC", "AAGAGGCAGCCAGG", "AAGAGGCAGGCCAC", "TACGCTGCGGCACG", "TACGCTGCTAACAA", "CAGAGAGGACAATA", "TAGGCATGCACGCG", "GTAGAGGACTATAT",  "CTCTCTACGAGGGC")] <- c(2, 3, 3, 1, 1, 0, 1,  1, 1, 0,  1, 1, 1, 1, 1, 1)
sce <- computeSumFactors(sce)
preclusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = preclusters)
sce <- logNormCounts(sce)
usethis::use_data(sce, compress = "xz", overwrite = T)


## Networks -------------------------------------------------------------

# normalise each clone separately using transformGamPoi
sce_list <- foreach(clone_name = clone_names,
                    .final = function(x) setNames(x, clone_names)) %do% {

                      # subset genes and cells
                      sce_clone <- sce[, sce$clone == clone_name]
                      print(dim(sce_clone))

                      # find genes that are only expressed in 1 cell and make that 1 cell 0 as well to avoid spurious correlations
                      gene_idx <- rowSums(counts(sce_clone) > 0) == 1 | rowSums(counts(sce_clone) > 0) == 2
                      print(clone_name)
                      print("Number of genes expressed in only 1 or 2 cell:")
                      print(sum(gene_idx))
                      counts(sce_clone)[gene_idx, ] <- 0

                      # log normalisation
                      sce_clone <- computeSumFactors(sce_clone)
                      preclusters <- quickCluster(sce_clone)
                      sce_clone <- computeSumFactors(sce_clone, clusters = preclusters)
                      sce_clone <- logNormCounts(sce_clone)

                      mat <- read_csv(paste0("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_scripts/toy_data_playaround4/gene9_cell15/GRNBoost2_files/count_matrices/", clone_name, ".csv"))
                      coln <- mat[[1]]
                      mat <- mat[, 2:ncol(mat)] %>% t()
                      colnames(mat) <- coln
                      assay(sce_clone, "rand_quantile_res") <- mat

                      sce_clone

                    }
usethis::use_data(sce_list, compress = "xz", overwrite = T)

# write input for GRNBoost2
write.table(rownames(sce_list[[1]]), "data-raw/intermediate_files/regulators.txt", quote = F, row.names = F, col.names = F)
source("/data/share/htp/hack_GRN/NPC_diff_network_analysis/functions/writeTableForGRNBoost2.R")
writeTableForGRNBoost2(sce_list, "rand_quantile_res", "data-raw/intermediate_files/count_matrices/")

# # run GRNBoost2
# system('sbatch data-raw/grnboost2_per_clone.sh "data-raw/intermediate_files/count_matrices" "data-raw/intermediate_files/regulators.txt" "inst/extdata"', wait = T)
system('cp "/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_scripts/toy_data_playaround4/gene9_cell15/GRNBoost2_files/output/"/* "inst/extdata"')

# load networks
network_list_raw <- loadGRNBoost2output("inst/extdata", clone_names)
usethis::use_data(network_list_raw, compress = "xz", overwrite = T)

# rescale interaction scores
network_list_scaled <- rescaleEdgeWeights(network_list_raw)
usethis::use_data(network_list_scaled, compress = "xz", overwrite = T)

# load gtfs
gtf_list <- list(human = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/hg38/genes.gtf"),
                 gorilla = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/gorGor6/genes.gtf"),
                 cynomolgus = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping/genomes/macFas6/genes.gtf"))

# filter
gtf_list <- foreach(gtf = gtf_list,
                    .final = function(x) {setNames(x, names(gtf_list))}) %do% {

                      gtf %>%
                        plyranges::filter(gene_name %in% genes) %>%
                        as_tibble() %>%
                        dplyr::select(seqnames, start, end, width, strand, source, type, score, phase, gene_id, gene_name, transcript_id, transcript_name) %>%
                        plyranges::as_granges()

                    }
usethis::use_data(gtf_list, compress = "xz", overwrite = T)

# remove gene pairs that overlap in any of the genomes
network_list_scaled_filt <- removeOverlappingGenePairs(network_list_scaled, gtf_list, clone2species, "gene_name")
usethis::use_data(network_list_scaled_filt, compress = "xz", overwrite = T)

# phylogenetic tree
tree <- read.tree("/data/share/htp/TRNP1/paper_data/Co-evolution-TRNP1-and-GI/protein/trees/mammaltree.txt") %>%
  drop.tip(.$tip.label[!.$tip.label %in% c("Homo_sapiens","Gorilla_gorilla", "Macaca_fascicularis")])
tree$tip.label <- c("cynomolgus", "gorilla", "human")
usethis::use_data(tree, compress = "bzip2", overwrite = T)

# consensus network
consensus_network <- createConsensus(network_list_scaled_filt, clone2species, tree)

# add consensus to network list
network_list_scaled_filt_withCons <- network_list_scaled_filt
network_list_scaled_filt_withCons[["consensus"]] <- consensus_network
usethis::use_data(network_list_scaled_filt_withCons, compress = "xz", overwrite = T)

# add directionality
network_list_scaled_filt_withCons_withDir <- addDirectionality(network_list_scaled_filt_withCons, sce_list, assay = "logcounts")
usethis::use_data(network_list_scaled_filt_withCons_withDir, compress = "xz", overwrite = T)

# save consensus network
consensus_network <- network_list_scaled_filt_withCons_withDir[["consensus"]]
usethis::use_data(consensus_network, compress = "xz", overwrite = T)

# save network list without consensus
network_list <- network_list_scaled_filt_withCons_withDir[names(network_list_scaled_filt_withCons_withDir) != "consensus"]
usethis::use_data(network_list, compress = "xz", overwrite = T)


## Module assignment ----------------------------------------------------

# initial modules
initial_modules <- assignInitialModules(consensus_network, regulators, N = 250)
usethis::use_data(initial_modules, compress = "xz", overwrite = T)

# pruned modules
pruned_modules <- pruneModules(initial_modules, consensus_network, "UIK_adj_kIM")
usethis::use_data(pruned_modules, compress = "xz", overwrite = T)

# random modules
random_modules <- createRandomModules(pruned_modules, genes)
usethis::use_data(random_modules, compress = "bzip2", overwrite = T)


## Eigengenes ----------------------------------------------------

# calculate eigengenes across all cells
eigengenes <- calculateEigengenes(regulators, pruned_modules, sce)
usethis::use_data(eigengenes, compress = "xz", overwrite = T)
plotEigengeneHeatmap(eigengenes)


# calculate eigengenes per species
eigengenes_per_species <- calculateEigengenes(regulators, pruned_modules, sce, per_species = T)
usethis::use_data(eigengenes_per_species, compress = "xz", overwrite = T)
plotEigengenesAlongPseudotime(eigengenes_per_species)

eigengenes_per_dir <- calculateEigengenes(regulators, pruned_modules, sce, "+-_separately")
plotEigengeneHeatmap(eigengenes_per_dir)

## Preservation statistics & tree reconstruction ----------------------------------------------

# calculate cor.kIM
pres_stats_jk <- calculatePresStats(pruned_modules, network_list, "cor_kIM", clone2species)
usethis::use_data(pres_stats_jk, compress = "xz", overwrite = T)
random_pres_stats_jk <- calculatePresStats(random_modules, network_list, "cor_kIM", clone2species)
usethis::use_data(random_pres_stats_jk, compress = "xz", overwrite = T)

# summarise jackknife values per module
pres_stats <- summarizeJackknifeStats(pres_stats_jk)
usethis::use_data(pres_stats, compress = "xz", overwrite = T)
random_pres_stats <- summarizeJackknifeStats(random_pres_stats_jk)
usethis::use_data(random_pres_stats, compress = "xz", overwrite = T)

# calculate distances
dist_jk <- convertPresToDist(pres_stats_jk, "cor_kIM")
usethis::use_data(dist_jk, compress = "xz", overwrite = T)
random_dist_jk <- convertPresToDist(random_pres_stats_jk, "cor_kIM")
usethis::use_data(random_dist_jk, compress = "xz", overwrite = T)

# reconstruct neighbour-joining trees
trees_jk <- reconstructTrees(dist_jk)
usethis::use_data(trees_jk, compress = "xz", overwrite = T)
random_trees_jk <- reconstructTrees(random_dist_jk)
usethis::use_data(random_trees_jk, compress = "xz", overwrite = T)

# calculate tree-based statistics
tree_stats_jk <- calculateTreeStats(trees_jk)
usethis::use_data(tree_stats_jk, compress = "xz", overwrite = T)
random_tree_stats_jk <- calculateTreeStats(random_trees_jk)
usethis::use_data(random_tree_stats_jk, compress = "xz", overwrite = T)

# summarise jackknife values per module
tree_stats <- summarizeJackknifeStats(tree_stats_jk, c("total_tree_length", "within_species_diversity"))
usethis::use_data(tree_stats, compress = "xz", overwrite = T)
random_tree_stats <- summarizeJackknifeStats(random_tree_stats_jk, c("total_tree_length", "within_species_diversity"))
usethis::use_data(random_tree_stats, compress = "xz", overwrite = T)


## Cross-species conservation -------------------------------------------


lm_overall <- fitTreeStatsLm(tree_stats, focus = "overall")
usethis::use_data(lm_overall, compress = "xz", overwrite = T)
module_conservation_overall <- findConservedDivergedModules(tree_stats, lm_overall)
usethis::use_data(module_conservation_overall, compress = "xz", overwrite = T)
plotConservedDivergedModules(module_conservation_overall)

target_conservation_POU5F1 <- findConservedDivergedTargets(tree_stats_jk, lm_overall, "POU5F1")
usethis::use_data(target_conservation_POU5F1, compress = "xz", overwrite = T)
plotConservedDivergedTargets(target_conservation_POU5F1)
