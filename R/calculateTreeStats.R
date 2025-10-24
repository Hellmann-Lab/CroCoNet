#' Calculate tree statistics across all modules
#'
#' Calculates various tree statistics (total tree length, within-species diversity, diversity, monophyleticity and subtree length of each species) for all modules/jackknifed module versions based on a list of tree representations.
#'
#' As part of the CroCoNet approach, pairwise module preservation scores are calculated between replicates, both within and across species (see \code{\link{calculatePresStats}}) and neighbor-joining trees are reconstructed based on these preservation scores per module (see \code{\link{convertPresToDist}} and \code{\link{reconstructTrees}}). The tips of the resulting tree represent the replicates and the branch lengths represent the dissimilarity of module connectivity patterns between the networks of 2 replicates.
#'
#' Various useful statistics can be defined based on these module trees:
#' \itemize{
#'  \item{Total tree length: The sum of all branch lengths in the tree; measures module variability both within and across species.}
#'  \item{Diversity of a species: The total length of the branches connecting the replicates of the given species to each other; measures module variability within this particular species.}
#'  \item{Within-species diversity: The sum of the diversity values across all species; measures module variability within species in general.}
#'  \item{Monophyleticity of a species: Indicates whether the tree is monophyletic for the replicates of the given species. Only if a module tree is monophyletic for a species of interest can the module be tested for divergence between this species and all others.}
#'  \item{Subtree length of a species: The sum of the branch lengths in the subtree that is defined by the replicates of the species and includes the internal branch connecting these replicates to the rest of the tree. Undefined if the tree is not monophyletic for the species of interest.}
#' }
#'
#' In the later steps of the pipeline, these tree-based statistics can be used to 1) identify modules that are conserved, diverged overall or diverged between a species and all others (see \code{\link{findConservedDivergedModules}}), and 2) pinpoint individual target genes within these modules that contribute the most to the conservation/divergence (see \code{\link{findConservedDivergedTargets}}).
#'
#' @param tree_list A named list of \code{\link{phylo}} objects containing the tree representations of all modules/jackknifed module versions, with the tips of the trees corresponding to replicates. All trees are expected to have a component \code{species} that specifies which species each tip belongs to. The trees can also contain a component \code{info} that stores metadata of the module/jackknifed module version in a data frame format.
#' @param n_cores Integer, the number of cores (default: 1).
#'
#' @return A data frame of tree statistics with the following columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Module size, the number of target genes assigned to a regulator (only present if the column is also present in the component \code{info} of the trees).}
#' \item{type}{Character, module type (orig = original or jk = jackknifed, only present if the input trees were reconstructed with jackknifing).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig', only present if the input trees were reconstructed with jackknifing).}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type 'orig', only present if the input trees were reconstructed with jackknifing).}
#' \item{within_species_diversity}{Numeric, the sum of the species-wise diversities.}
#' \item{total_tree_length}{Numeric, the total length of all branches in the tree.}
#' \item{\{\{species\}\}_diversity}{Numeric, as many columns as there are species, each of them containing the the total length of the branches connecting the replicates of the given species.}
#' \item{\{\{species\}\}_monophyl}{Logical, as many columns as there are species, each of them indicating whether the tree is monophyletic for the replicates of the given species.}
#' \item{\{\{species\}\}_subtree_length}{Numeric, as many columns as there are species, each of them containing the sum of the branch lengths in the subtree that is defined by the replicates of the species and includes the internal branch connecting these replicates to the rest of the tree. NA if the tree is not monophyletic for the replicates of the given species.}
#' }
#' @export
#' @examples tree_stats_jk <- calculateTreeStats(trees_jk)
calculateTreeStats <- function(tree_list, n_cores = 1L) {

  if (!inherits(tree_list, "list"))
    stop("The argument \"tree_list\" should be a named list.")

  if (is.null(names(tree_list)))
    stop("The argument \"tree_list\" should be a named list.")

  if (any(!sapply(tree_list, function(tree) {inherits(tree, "phylo")})))
    stop("All elements of \"tree_list\" should be of class \"phylo\".")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") && !(inherits(n_cores, "numeric") && n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  warning_single_replicate <- character()

  for (id in names(tree_list)) {

    tree <- tree_list[[id]]

    if (is.null(tree$species))
      stop(paste0("No species information found for tree \"", id, "\". Please add a component \"species\" to all trees specifying which species each tip belongs to."))

    if (length(unique(tree$species)) < 2)
      stop(paste0("Only 1 species found in tree \"", id, "\". For the calculation of tree statistics, each tree has to consist of at least 2 species."))

    spec_count <- table(tree$species)
    if (all(spec_count == 1))
      stop(paste0("Tree \"", id, "\" contains only 1 replicate per species, therefore diversity cannot be estimated."))

    if (any(spec_count == 1))
      warning_single_replicate <- c(warning_single_replicate, paste0("One or more trees contain only 1 replicate of species \"", paste(names(spec_count)[spec_count == 1], collapse = '\" and \"'), "\". The diversity of this/these species cannot be estimated and the within-species diversity will be estimated using the other species only."))

    if (any(tree$edge_lengths < 0))
      warning(paste0("Negative branch lengths found in tree \"", id, "\"."))

  }

  if (length(warning_single_replicate) > 0)
    warning(paste(unique(warning_single_replicate), collapse = "\n"))

  id = NULL

  doParallel::registerDoParallel(n_cores)

  tree_stats <-

    foreach::foreach(id = names(tree_list),
                     .combine = rbind) %dopar%

    {

      tree <- tree_list[[id]]
      if (is.null(tree$info))
        tree$info <- data.frame(id = id)
      getTreeStats(tree)

    }

  doParallel::stopImplicitCluster()

  tree_stats %>%
    tidyr::pivot_wider(names_from = "spec", names_glue = c("{spec}_{.value}"), values_from = c("diversity", "monophyl", "subtree_length"))

}




#' Calculate diversity
#'
#' @description Calculates the diversity within a species given a module tree where the tips correspond to replicates. This diversity is defined as the total length of the branches connecting the replicates of that species.
#' @param tree Object of class \code{\link{phylo}} containing the tree representation of a module/jackknifed module version, with the tips of the tree corresponding to replicates. The tree is expected to have a component \code{species} that specifies which species each tip belongs to.
#' @param spec Character, the name of the species.
#'
#' @return Numeric, the diversity within the specified species.
#' @noRd
calculateDiversity <- function(tree, spec) {

  # replicates belonging to the species
  replicates_spec <- tree$tip.label[tree$species == spec]

  if (length(replicates_spec) == 1)
    return(NA)

  tree %>%
    ape::keep.tip(replicates_spec) %>%
    tidytree::as_tibble() %>%
    dplyr::pull(.data[["branch.length"]]) %>%
    sum(na.rm = TRUE)

}


#' Get monophyleticity
#'
#' @description Determines whether a module tree is monophyletic for a species.
#' @param tree Object of class \code{\link{phylo}} containing the tree representation of a module/jackknifed module version, with the tips of the tree corresponding to replicates. The tree is expected to have a component \code{species} that specifies which species each tip belongs to.
#' @param spec Character, the name of the species.
#'
#' @return Logical indicating whether the input tree is monophyletic for the specified species.
#' @noRd
isMonophyletic <- function(tree, spec) {

  # replicates belonging to the species
  replicates_spec <- tree$tip.label[tree$species == spec]

  if (length(replicates_spec) == 1 || length(replicates_spec) == length(tree$tip.label) - 1)
    return(NA)

  ape::is.monophyletic(tree, tips = replicates_spec)

}


#' Calculate the subtree length of a species
#'
#' @description Calculates the subtree length of a species given a module tree where the tips correspond to replicates. The subtree length of a species is defined as the sum of the branch lengths in the subtree that is defined by the replicates of the species and includes the internal branch connecting these replicates to the rest of the tree.
#' @param tree Object of class \code{\link{phylo}} containing the tree representation of a module/jackknifed module version, with the tips of the tree corresponding to replicates. The tree is expected to have a component \code{species} that specifies which species each tip belongs to.
#' @param spec Character, the name of the species.
#'
#' @return Numeric, the subtree length of the specified species.
#' @noRd
calculateSpeciesSubtreeLength <- function(tree, spec) {

  if (!ape::is.monophyletic(tree, tips = tree$tip.label[tree$species == spec])) {

    return(NA)

  } else {

    # replicates belonging to the species
    replicates_spec <- tree$tip.label[tree$species == spec]

    if (length(replicates_spec) == 1)
      return(NA)

    # replicates NOT belonging to the species
    replicates_rest <- tree$tip.label[tree$species != spec]

    if (length(replicates_rest) == 1)
      return(NA)

    # reroot tree
    tree_rooted <- tree %>%
      ape::root(outgroup = replicates_rest, resolve.root = TRUE)

    # get the node numbers of the replicates belonging to the species
    nodes_spec <- which(tree_rooted$tip.label %in% replicates_spec)

    # get the paths from the root to each replicate
    paths <- ape::nodepath(tree_rooted)[nodes_spec]

    # collect all unique edges along those paths
    keep_nodes <- unique(unlist(paths))
    keep_edges <- which(tree_rooted$edge[,2] %in% keep_nodes)

    # sum the branch lengths of those edges
    sum(tree_rooted$edge.length[keep_edges])

  }

}


#' Calculate tree statistics
#'
#' @description Calculates various tree statistics (total tree length, within-species diversity, diversity, monophyleticity and subtree length of each species) for a single module based on a tree where the tips correspond to replicates and the branch lengths correspond to dissimilarity of module topology between the replicates.
#' @param tree Object of class \code{\link{phylo}} containing the tree representation of a module/jackknifed module version, with the tips of the tree corresponding to replicates. The tree is expected to have a component \code{species} that specifies which species each tip belongs to.
#'
#' @return A data frame of tree statistics with the following columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Module size, the number of target genes assigned to a regulator (only present if the column is also present in the component \code{info} of the trees).}
#' \item{type}{Character, module type (orig = original or jk = jackknifed, only present if the input trees were reconstructed with jackknifing).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig', only present if the input trees were reconstructed with jackknifing).}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type 'orig', only present if the input trees were reconstructed with jackknifing).}
#' \item{within_species_diversity}{Numeric, the sum of the species-wise diversities.}
#' \item{total_tree_length}{Numeric, the sum of all branch lengths in the tree.}
#' \item{species}{Character, the name of the species.}
#' \item{diversity}{Numeric, the sum of the branch lengths in the subtree that contains only the replicates of the given species.}
#' \item{monophyl}{Logical indicating whether the tree is monophyletic for the replicates of the given species.}
#' \item{subtree_length}{Numeric, the sum of the branch lengths in the subtree that is defined by the replicates of the species and includes the internal branch connecting these replicates to the rest of the tree.}
#' }
#' @noRd
getTreeStats <- function (tree) {

  tree$info %>%
    tidyr::expand_grid(spec = unique(tree$species)) %>%
    dplyr::group_by(.data[["spec"]]) %>%
    dplyr::mutate(diversity = calculateDiversity(tree, .data[["spec"]]),
                  monophyl = isMonophyletic(tree, .data[["spec"]]),
                  subtree_length = calculateSpeciesSubtreeLength(tree, .data[["spec"]])) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(total_tree_length = tree %>%
                    tidytree::as_tibble() %>%
                    dplyr::pull(.data[["branch.length"]]) %>%
                    sum(na.rm = TRUE),
                  within_species_diversity = sum(.data[["diversity"]]),
                  .after = ncol(tree$info))

}

