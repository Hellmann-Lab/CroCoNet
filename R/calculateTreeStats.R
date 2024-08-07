#' Calculate tree statistics across all modules
#'
#' Calculates various tree statistics (total tree length, total diversity, diversity within each species, monophyleticity for each species and branch length between each species and all others) for all modules based on a list of trees. The tips of the trees should correspond to clones and the branch lengths of the trees should correspond to dissimilarity of module topology between the clones.
#' @param tree_list A list of trees per module as \code{\link{phylo}} objects with clones as tips. All trees are expected to have a component 'species' that specifies which species each tip belongs to.
#' @param n_cores Number of cores.
#'
#' @return A data frame of tree statistics.
#' @export
#' @examples tree_stats_jk <- calculateTreeStats(trees_jk)
calculateTreeStats <- function(tree_list, n_cores = 1L) {

  if (!inherits(tree_list, "list"))
    stop("The argument \"tree_list\" should be a named list.")

  if (is.null(names(tree_list)))
    stop("The argument \"tree_list\" should be a named list.")

  if (any(!sapply(tree_list, function(tree) {inherits(tree, "phylo")})))
    stop("All elements of \"tree_list\" should be of class \"phylo\".")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") & !(inherits(n_cores, "numeric") & n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  warning_single_clone <- character()

  for (id in names(tree_list)) {

    tree <- tree_list[[id]]

    if (is.null(tree$species))
      stop(paste0("No species information found for tree \"", id, "\". Please add a component \"species\" to all trees specifying which species each tip belongs to."))

    if (length(unique(tree$species)) < 2)
      stop(paste0("Only 1 species found in tree \"", id, "\". For the calculation of tree statistics, each tree has to consist of at least 2 species."))

    spec_count <- table(tree$species)
    if (all(spec_count == 1))
      stop(paste0("Tree \"", id, "\" contains only 1 clone per species, therefore diversity cannot be estimated."))

    if (any(spec_count == 1))
      warning_single_clone <- c(warning_single_clone, paste0("One or more trees contain only 1 clone of species \"", paste(names(spec_count)[spec_count == 1], collapse = '\" and \"'), "\". The diversity of this/these species cannot be estimated and the within-species diversity will be estimated using the other species only."))

    if (any(tree$edge_lengths < 0))
      warning(paste0("Negative branch lengths found in tree \"", id, "\"."))

  }

  if (length(warning_single_clone) > 0)
    warning(paste(unique(warning_single_clone), collapse = "\n"))

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
    tidyr::pivot_wider(names_from = "spec", names_glue = c("{spec}_{.value}"), values_from = c("diversity", "monophyl", "to_other_branch_length"))

}




#' Calculate diversity
#'
#' @description Calculates the diversity within a species given a module tree where the tips correspond to clones. This diversity is defined as the total lengths of the branches connecting the clones of that species.
#' @param tree Object of class \code{\link{phylo}} with clones as tips. It is expected to have a component 'species' that specifies which species each tip belongs to.
#' @param spec Character, the name of the species.
#'
#' @return Numeric, the diversity within a species.
#' @noRd
calculateDiversity <- function(tree, spec) {

  # clones belonging to the species
  clones_spec <- tree$tip.label[tree$species == spec]

  if (length(clones_spec) == 1)
    return(NA)

  tree %>%
    ape::keep.tip(clones_spec) %>%
    tidytree::as_tibble() %>%
    dplyr::pull(.data[["branch.length"]]) %>%
    sum(na.rm = TRUE)

}


#' Get monophyleticity
#'
#' @description Determines whether a module tree is monophyletic for a certain species.
#' @param tree Object of class \code{\link{phylo}} with clones as tips. It is expected to have a component 'species' that specifies which species each tip belongs to.
#' @param spec Character, the name of the species.
#'
#' @return Logical indicating whether the input tree is monophyletic for the specified species.
#' @noRd
isMonophyletic <- function(tree, spec) {

  # clones belonging to the species
  clones_spec <- tree$tip.label[tree$species == spec]

  if (length(clones_spec) == 1 || length(clones_spec) == length(tree$tip.label) - 1)
    return(NA)

  ape::is.monophyletic(tree, tips = clones_spec)

}


#' Calculate the branch length between a species and all others
#'
#' @description Calculates the species-to-other branch length given a module tree where the tips correspond to clones. The species-to-other branch length is defined as the length of the branch connecting the subtree of the clones belonging to the species of interest and the subtree of all other clones.
#' @param tree Object of class \code{\link{phylo}} with clones as tips. It is expected to have a component 'species' that specifies which species each tip belongs to.
#' @param spec Character, the name of the species.
#'
#' @return Numeric, the branch length between the species of interest and the rest of the tree.
#' @noRd
calculateSpeciesToOtherBranchLength <- function(tree, spec) {

  if (!ape::is.monophyletic(tree, tips = tree$tip.label[tree$species == spec])) {

    return(NA)

  } else {

    # clones belonging to the species
    clones_spec <- tree$tip.label[tree$species == spec]

    if (length(clones_spec) == 1)
      return(NA)

    # clones NOT belonging to the species
    clones_rest <- tree$tip.label[tree$species != spec]

    if (length(clones_rest) == 1)
      return(NA)

    # total number of clones
    n_clones <- length(tree$tip.label)

    # reroot tree
    tree_rooted <- tree %>%
      ape::root(outgroup = clones_spec, resolve.root = TRUE)

    # convert tree to data frame while preserving the order of edges (will be handy for plotting)
    tree_df_rooted <- getTreeDf(tree_rooted)

    # get the most recent common ancestor of the clones (it'll have 2 node numbers because it's the root, but these 2 nodes have a distance of 0, 1st node number: the number right after the tips, 2nd node number: the number after the last node (note: n_nodes = 2*n_tips - 2))
    mrca_spec <- c(n_clones + 1, 2*n_clones - 1)

    # get the most recent common ancestor of the gorilla and cynomolgus clones
    mrca_rest <- ape::getMRCA(tree_rooted, clones_rest)

    # flag the internal branch that connect the human to the gorilla+cyno group
    tree_df_rooted <- tree_df_rooted %>%
      dplyr::mutate(is_internal = (.data[["parent"]] %in% mrca_spec & .data[["node"]] == mrca_rest) | (.data[["node"]] %in% mrca_spec & .data[["parent"]] == mrca_rest))

    # get internal branch length
    tree_df_rooted %>%
      dplyr::filter(.data[["is_internal"]]) %>%
      dplyr::pull(.data[["branch.length"]])

  }

}


#' Calculate tree statistics
#'
#' @description Calculates various tree statistics (total tree length, total diversity, diversity within each species, monophyleticity for each species and branch length between each species and all others) for a single module based on a tree where the tips correspond to clones and the branch lengths correspond to dissimilarity of module topology between the clones.
#' @param tree Object of class \code{\link{phylo}} with clones as tips. It is expected to have a component 'species' that specifies which species each tip belongs to.
#'
#' @return A data frame of tree statistics.
#' @noRd
getTreeStats <- function (tree) {

  tree$info %>%
    tidyr::expand_grid(spec = unique(tree$species)) %>%
    dplyr::group_by(.data[["spec"]]) %>%
    dplyr::mutate(diversity = calculateDiversity(tree, .data[["spec"]]),
                  monophyl = isMonophyletic(tree, .data[["spec"]]),
                  to_other_branch_length = calculateSpeciesToOtherBranchLength(tree, .data[["spec"]])) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(total_tree_length = tree %>%
                    tidytree::as_tibble() %>%
                    dplyr::pull(.data[["branch.length"]]) %>%
                    sum(na.rm = TRUE),
                  within_species_diversity = sum(.data[["diversity"]]),
                  .after = ncol(tree$info))

}

