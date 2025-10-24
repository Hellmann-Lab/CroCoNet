library(CroCoNet)
library(foreach)
library(tidyverse)
library(scran)
library(transformGamPoi)
library(ape)
library(SingleCellExperiment)
library(foreach)



## Genes ----------------------------------------------------------

# list of regulators
regulators <- sort(c("SOX2", "POU5F1", "ONECUT2", "TCF4", "BNC2", "TEAD1", "SMAD2", "ZNF701", "SOX15", "CREB3", "PROX1", "FOS"))
usethis::use_data(regulators, compress = "bzip2", overwrite = T)

# genes in the modules of 7 famous TFs (and all other TFs excluded)
genes <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_scripts/toy_data_playaround2/gene8_cell1/RDS/genes.rds")
usethis::use_data(genes, compress = "bzip2", overwrite = T)

# list of TRs
jaspar_core_TRs <- data.table::fread("grep '>' data-raw/intermediate_files/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt",
                                     col.names = c("motif_id","SYMBOL")) %>%
  separate_rows(SYMBOL, sep="::") %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort() %>%
  toupper()
usethis::use_data(jaspar_core_TRs, compress = "bzip2", overwrite = T)

jaspar_unvalidated_TRs <- data.table::fread("grep '>' data-raw/intermediate_files/JASPAR2024_UNVALIDATED_non-redundant_pfms_jaspar.txt",
                                            col.names = c("motif_id","SYMBOL")) %>%
  separate_rows(SYMBOL, sep="::") %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort() %>%
  toupper()
usethis::use_data(jaspar_unvalidated_TRs, compress = "bzip2", overwrite = T)

image_TRs <- data.table::fread("data-raw/intermediate_files/IMAGE_motifs.txt",
                               col.names = c("SYMBOL", "motif_id", "Evidence")) %>%
  pull(SYMBOL) %>%
  unique() %>%
  sort()
usethis::use_data(image_TRs, compress = "bzip2", overwrite = T)

## Clone 2 species --------------------------------------------------------

# read and save object
replicate2species <- readRDS("/data/share/htp/hack_GRN/CroCoNet_scripts_and_data/data/neural_differentiation_dataset/processed_data/replicate2species.rds")
usethis::use_data(replicate2species, compress = "bzip2", overwrite = T)

# replicate names
replicate_names <- replicate2species$replicate


## SCE object -----------------------------------------------------------

# randomly select 100 cells per replicate
sce <- readRDS("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_scripts/toy_data_playaround2/gene8_cell1/RDS/sce.rds")

sce_list <- foreach(replicate_name = replicate_names,
                    .final = function(x) setNames(x, replicate_names)) %do% {

                      # subset genes and cells
                      sce[, sce$replicate == replicate_name]

                    }

cells_sox2 <- lapply(replicate2species$replicate[replicate2species$species %in% c("human", "cynomolgus")], function(replicate) {

  set.seed(0)
  sample(colData(sce_list[[replicate]]) %>% as.data.frame() %>% arrange(desc(pseudotime)) %>% rownames() %>% .[1:15], 13)

}) %>% unlist()
counts(sce)["SOX2", cells_sox2] <- counts(sce)["SOX2", cells_sox2] + 1


cells_sox2_big <- lapply(replicate2species$replicate[replicate2species$species %in% c("human", "cynomolgus")], function(replicate) {

  sce_replicate <- sce[,sce$replicate == replicate]
  colnames(sce_replicate)[counts(sce_replicate)["SOX2", ] >= 8]

}) %>% unlist()

counts(sce)["SOX2", cells_sox2_big] <- counts(sce)["SOX2", cells_sox2_big] - round((counts(sce)["SOX2", cells_sox2_big] - 2) / 2)

cells_FOS <- lapply(replicate2species$replicate[replicate2species$species %in% c("gorilla", "human")], function(replicate) {

  set.seed(0)
  c <- colData(sce_list[[replicate]]) %>% as.data.frame() %>% filter(pseudotime < 0.5) %>% rownames()
  sample(c, length(c) - 10)

}) %>% unlist()
counts(sce)["FOS", cells_FOS] <- counts(sce)["FOS", cells_FOS] - round((counts(sce)["FOS", cells_FOS]) / 3 * 2)

cells_FOS_small <- lapply(replicate2species$replicate[replicate2species$species %in% c( "human")], function(replicate) {

  set.seed(0)
  c <- colData(sce_list[[replicate]]) %>% as.data.frame() %>% filter(pseudotime > 0.75) %>% rownames()
  sample(c, round(length(c) /8 ))

}) %>% unlist()
counts(sce)["FOS", cells_FOS_small] <- counts(sce)["FOS", cells_FOS_small] + 1

cells_SMAD2_small <- lapply(replicate2species$replicate[replicate2species$species %in% c( "human", "gorilla")], function(replicate) {

  set.seed(0)
  c <- colData(sce_list[[replicate]]) %>% as.data.frame() %>% filter(pseudotime > 0.6) %>% rownames()
  sample(c, round(length(c) /3 ))

}) %>% unlist()
cells_SMAD2_small2 <- lapply(replicate2species$replicate[replicate2species$species %in% c("cynomolgus")], function(replicate) {

  set.seed(0)
  c <- colData(sce_list[[replicate]]) %>% as.data.frame() %>% filter(pseudotime > 0.6) %>% rownames()
  sample(c, round(length(c) /7 ))

}) %>% unlist()
counts(sce)["SMAD2", c(cells_SMAD2_small, cells_SMAD2_small2)] <- counts(sce)["SMAD2", c(cells_SMAD2_small,cells_SMAD2_small2)] + 1

cells_TEAD1_small <- lapply(replicate2species$replicate[replicate2species$species %in% c( "gorilla")], function(replicate) {

  set.seed(0)
  c <- colData(sce_list[[replicate]]) %>% as.data.frame() %>% filter(pseudotime > 0.65) %>% rownames()
  sample(c, round(length(c) /4 ))

}) %>% unlist()
counts(sce)["TEAD1", cells_TEAD1_small] <- counts(sce)["TEAD1", cells_TEAD1_small] + 1

cells_ZNF701 <- lapply(replicate2species$replicate[replicate2species$species %in% c("gorilla", "human")], function(replicate) {

  set.seed(0)
  c <- colData(sce_list[[replicate]]) %>% as.data.frame() %>% filter(pseudotime < 0.5) %>% rownames()
  sample(c, length(c) / 4)

}) %>% unlist()
counts(sce)["ZNF701", cells_ZNF701] <- counts(sce)["ZNF701", cells_ZNF701] + 1

cells_ZNF701 <- lapply(replicate2species$replicate[replicate2species$species %in% c("cynomolgus")], function(replicate) {

  set.seed(0)
  c <- colData(sce_list[[replicate]]) %>% as.data.frame() %>% filter(pseudotime > 0.5) %>% rownames()
  sample(c, length(c) - 10)

}) %>% unlist()
counts(sce)["ZNF701", cells_ZNF701] <- counts(sce)["ZNF701", cells_ZNF701] - round((counts(sce)["ZNF701", cells_ZNF701]) / 3 * 2)

sce <- computeSumFactors(sce)
preclusters <- quickCluster(sce, min.size = 50)
sce <- computeSumFactors(sce, clusters = preclusters)
sce <- logNormCounts(sce)
usethis::use_data(sce, compress = "xz", overwrite = T)


## Networks -------------------------------------------------------------

# normalise each replicate separately using transformGamPoi
sce_list <- foreach(replicate_name = replicate_names,
                    .final = function(x) setNames(x, replicate_names)) %do% {

                      # subset genes and cells
                      sce_replicate <- sce[, sce$replicate == replicate_name]
                      print(dim(sce_replicate))

                      # find genes that are only expressed in 1 cell and make that 1 cell 0 as well to avoid spurious correlations
                      gene_idx <- rowSums(counts(sce_replicate) > 0) == 1 | rowSums(counts(sce_replicate) > 0) == 2
                      print(replicate_name)
                      print("Number of genes expressed in only 1 or 2 cell:")
                      print(sum(gene_idx))
                      counts(sce_replicate)[gene_idx, ] <- 0

                      # log normalisation
                      sce_replicate <- computeSumFactors(sce_replicate)
                      preclusters <- quickCluster(sce_replicate, min.size = 50)
                      sce_replicate <- computeSumFactors(sce_replicate, clusters = preclusters)
                      sce_replicate <- logNormCounts(sce_replicate)

                      mat <- read_csv(paste0("/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_scripts/toy_data_playaround2/gene8_cell1/GRNBoost2_files/count_matrices/", replicate_name, ".csv"))
                      coln <- mat[[1]]
                      mat <- mat[, 2:ncol(mat)] %>% t()
                      colnames(mat) <- coln
                      assay(sce_replicate, "rand_quantile_res") <- mat

                      sce_replicate

                    }
# usethis::use_data(sce_list, compress = "xz", overwrite = T)

# write input for GRNBoost2
write.table(rownames(sce_list[[1]]), "data-raw/intermediate_files/regulators.txt", quote = F, row.names = F, col.names = F)
source("/data/share/htp/hack_GRN/NPC_diff_network_analysis/functions/writeTableForGRNBoost2.R")
writeTableForGRNBoost2(sce_list, "rand_quantile_res", "data-raw/intermediate_files/count_matrices/")

# # run GRNBoost2
# system('sbatch data-raw/grnboost2_per_replicate.sh "data-raw/intermediate_files/count_matrices" "data-raw/intermediate_files/regulators.txt" "inst/extdata"', wait = T)
system('cp "/data/share/htp/hack_GRN/NPC_diff_network_analysis/additional_scripts/toy_data_playaround2/gene8_cell1/GRNBoost2_files/output/"/* "inst/extdata"')

# load networks
network_list_raw <- loadNetworks("~/CroCoNet_backup/inst/extdata", replicate_names, rep = 10)
usethis::use_data(network_list_raw, compress = "xz", overwrite = T)

# load gtfs
gtf_list <- list(human = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping_and_QC/genomes/hg38/genes.gtf"),
                 gorilla = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping_and_QC/genomes/gorGor6/genes.gtf"),
                 cynomolgus = plyranges::read_gff("/data/share/htp/hack_GRN/NPC_diff_network_analysis/01.mapping_and_QC/genomes/macFas6/genes.gtf"))

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

# rescale interaction scores
# remove gene pairs that overlap in any of the genomes
network_list <- network_list_raw %>%
  removeOverlappingGenePairs(gtf_list, replicate2species, "gene_name") %>%
  normalizeEdgeWeights()
usethis::use_data(network_list, compress = "xz", overwrite = T)

# phylogenetic tree
tree <- read.tree("/data/share/htp/TRNP1/paper_data/Co-evolution-TRNP1-and-GI/protein/trees/mammaltree.txt") %>%
  drop.tip(.$tip.label[!.$tip.label %in% c("Homo_sapiens","Gorilla_gorilla", "Macaca_fascicularis")])
tree$tip.label <- c("cynomolgus", "gorilla", "human")
usethis::use_data(tree, compress = "bzip2", overwrite = T)

# consensus network
consensus_network <- network_list %>%
  createConsensus(replicate2species, tree) %>%
  addDirectionality(sce)
usethis::use_data(consensus_network, compress = "xz", overwrite = T)


## Module assignment ----------------------------------------------------

# initial modules
initial_modules <- assignInitialModules(consensus_network, regulators, N = 250)
usethis::use_data(initial_modules, compress = "xz", overwrite = T)

# pruned modules
pruned_modules <- pruneModules(initial_modules, "UIK_adj_kIM", consensus_network)
usethis::use_data(pruned_modules, compress = "xz", overwrite = T)

# random modules
random_modules <- createRandomModules(pruned_modules, genes)
usethis::use_data(random_modules, compress = "bzip2", overwrite = T)


## Eigengenes ----------------------------------------------------

# calculate eigengenes across all cells
eigengenes <- calculateEigengenes(pruned_modules, sce)
usethis::use_data(eigengenes, compress = "xz", overwrite = T)
plotEigengeneHeatmap(eigengenes)

# calculate eigengenes per species
eigengenes_per_species <- calculateEigengenes(pruned_modules, sce, per_species = T)
usethis::use_data(eigengenes_per_species, compress = "xz", overwrite = T)
plotEigengenesAlongPseudotime(eigengenes_per_species)

eigengenes_per_dir <- calculateEigengenes(pruned_modules, sce, "+-_separately")
plotEigengeneHeatmap(eigengenes_per_dir)

## Preservation statistics & tree reconstruction ----------------------------------------------

# calculate cor.kIM
pres_stats_jk <- calculatePresStats(pruned_modules, network_list, c("cor_adj", "cor_kIM"), replicate2species)
usethis::use_data(pres_stats_jk, compress = "xz", overwrite = T)
random_pres_stats_jk <- calculatePresStats(random_modules, network_list, c("cor_adj", "cor_kIM"), replicate2species)
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

# distances of the original modules
dist <- dist_jk[paste0(regulators, "_orig")]
names(dist) <- regulators
dist <- lapply(dist, function(df) {

  df$type = df$id = df$gene_removed <- NULL
  df

})
usethis::use_data(dist, compress = "xz", overwrite = T)

# reconstruct neighbour-joining trees
trees_jk <- reconstructTrees(dist_jk)
usethis::use_data(trees_jk, compress = "xz", overwrite = T)
random_trees_jk <- reconstructTrees(random_dist_jk)
usethis::use_data(random_trees_jk, compress = "xz", overwrite = T)

# trees of the original modules
trees <- trees_jk[paste0(regulators, "_orig")]
names(trees) <- regulators
trees <- lapply(trees, function(tree) {

  tree$info$type = tree$info$id = tree$info$gene_removed <- NULL
  tree

})
usethis::use_data(trees, compress = "xz", overwrite = T)

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

POU5F1_target_conservation <- findConservedDivergedTargets("POU5F1", tree_stats_jk, lm_overall)
usethis::use_data(POU5F1_target_conservation, compress = "xz", overwrite = T)
plotConservedDivergedTargets(POU5F1_target_conservation)
