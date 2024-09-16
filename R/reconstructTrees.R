#' Convert preservation statistics to distances across all modules
#'
#' Converts preservation statistics between clones to distance measures ranging from 0 to 1.
#'
#' As part of the CroCoNet approach, pairwise module preservation scores are calculated between clones, both within and across species (see \code{\link{calculatePresStats}}) to gain information about the cross-species differences but also about the within-species diversity of the modules. These correlation-based preservation statistics quantify how well the module topology is preserved between the networks of two clones.
#'
#' This function converts a chosen preservation statistic (\emph{p}) specified by the argument \code{stat} into a distance measure (\emph{d}) using the following formula: \deqn{d = \frac{max(p) - p}{max(p) - min(p)}}
#' If the theoretical minimum and maximum of the preservation statistic are known, these can be provided using the parameters \code{min_stat} and \code{max_stat}. As the preservation statistics implemented by CroCoNet are correlation-based, they all range between -1 and 1, and thus  \code{min_stat} and \code{max_stat} should be set to -1 and 1, respectively (default). If a custom preservation statistic is used, \code{min_stat} and \code{max_stat} might have to be set to different values. If they are set to NULL or ±Inf, the minimum and maximum of the preservation statistic are calculated empirically using the data.
#'
#' The function also splits up the distance measures into modules/jackknife module versions (depending on whether the preservation statistics were calculated with or without jackknifing, see see the parameter \code{jackknife} in the function \code{\link{calculatePresStats}}) and outputs a list of data frames per module/jackknife module version. Modules/jackknifed module versions where the preservation statistic for any of the clone pairs is NA are removed.
#'
#' In the next step of the pipeline, these distance measures are used to reconstruct a neighbor-joining tree per module/jackknife module version that represent the dissimilarity of module topology across all clones (see \code{\link{reconstructTrees}}).
#'
#' @param pres_stats Data frame of the preservation statistics, required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type ("orig" = original or "jk" = jackknifed, only needed if the preservation statistics were calculated with jackknifing).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type "jk" and nameOfRegulator_orig in case of module type "orig", only needed if the preservation statistics were calculated with jackknifing).}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type "orig", only needed if the preservation statistics were calculated with jackknifing).}
#' \item{clone1, clone2}{Character the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species \code{clone1} and \code{clone2} belong to, respectively.}
#' \item{\{\{stat\}\}}{Numeric, the preservation statistic specified by the parameter \code{stat}.}
#' }
#' @param stat Character, the name of the column containing the preservation statistic that is to be converted into a distance measure.
#' @param min_stat,max_stat Numeric, the theoretical minimum and maximum value of \code{stat} (default: -1 and 1, respectively). For the preservation statistics implemented by CroConet and all other correlation-based statistics, please leave \code{min_stat} and \code{max_stat} at -1 and 1, respectively. For custom preservation statistics, a different value might have to be used. If set to NULL or ±Inf, the conversion to a distance measure is performed using the empirical minimum and maximum.
#' @param n_cores Integer, the number of cores (default: 1).
#'
#' @return A named list containing the distance measures as data frames for all modules or jackknifed module versions in the input \code{pres_stats}. The data frames contain 1 new column in addition to the relevant columns of \code{pres_stats}:
#' \describe{
#' \item{dist}{Numeric, the distance measure ranging from 0 to 1 calculated based on \code{stat}.}
#' }
#' While \code{pres_stats} contains only non-redundant clone pairs, the data frames in the output contain all possible clone pairs (i.e. cloneA-cloneB and cloneB-cloneA are 2 separate entries with the same distance).
#' @export
#'
#' @examples dist_jk <- convertPresToDist(pres_stats_jk, "cor_kIM")
convertPresToDist <- function(pres_stats, stat, min_stat = -1, max_stat = 1, n_cores = 1L) {

  # check input data
  if (!is.data.frame(pres_stats))
    stop("The argument \"pres_stats\" should be a data frame.")

  if (any(!c("regulator", "module_size", "clone1", "clone2", "species1", "species2") %in% colnames(pres_stats)))
    stop("The argument \"pres_stats\" should contain the columns \"regulator\", \"module_size\", \"clone1\", \"clone2\", \"species1\" and \"species2\".")

  if (!inherits(stat, "character") || length(stat) != 1)
    stop("The argument \"stat\" should be a string specifying the statistic that is to be converted into a distance measure.")

  if (!stat %in% colnames(pres_stats))
    stop("The statistic specified by \"stat\" was not found in \"pres_stats\".")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") && !(inherits(n_cores, "numeric") && n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  module_id = NULL

  # number of clones
  n_clones <- length(unique(c(pres_stats$clone1, pres_stats$clone2)))

  # ID column: "regulator" if preservation stats were calculated without jackknifing, "id" if preservation stats were calculated with jackknifing
  id_column <- ifelse("id" %in% colnames(pres_stats), "id", "regulator")

  # remove modules/jackknifed module versions where the preservation statistic for any of the clone pairs is NA
  pres_stats <- pres_stats %>%
    dplyr::rename(stat = paste(stat)) %>%
    tidyr::drop_na(.data[["stat"]]) %>%
    dplyr::group_by(.data[[id_column]]) %>%
    dplyr::filter(length(.data[["stat"]]) == n_clones*(n_clones - 1)/2) %>%
    dplyr::ungroup()

  # define the minimum and maximum value of the preservation statistic that correspond to a distance of 1 and 0, respectively
  # if min_stat or max_stat are NULL or ±Inf, take the observed minimum or maximum
  if (is.null(min_stat) || min_stat == -Inf) min_stat <- min(pres_stats$stat)
  if (is.null(max_stat) || max_stat == Inf) max_stat <- max(pres_stats$stat)

  # all modules/jackknifed module versions
  module_ids <- unique(pres_stats[[id_column]])

  doParallel::registerDoParallel(n_cores)

  # convert preservation statistics into distance measures and split data frame into a list of data frames
  dist_list <- foreach::foreach(module_id = module_ids) %dopar%

    {

      pres_stats_id <- pres_stats %>%
        dplyr::filter(.data[[id_column]] == module_id)
      getDist(pres_stats_id, min_stat, max_stat)

    }

  doParallel::stopImplicitCluster()

  names(dist_list) <- module_ids

  dist_list

}


#' Reconstruct trees across all modules
#'
#' Reconstructs neighbor-joining trees based on pairwise distance measures between clones.
#'
#' As part of the CroCoNet approach, pairwise module preservation scores are calculated between clones, both within and across species, to quantify how similar module connectivity patterns are between the networks of two clones (see \code{\link{calculatePresStats}}). These preservation scores are then converted into distance measures (see \code{\link{convertPresToDist}}).
#'
#' This function first sorts the distance measures of each module/jackknifed module version into a distance matrix of all clones, then based on this distance matrix reconstructs a tree using the neighbor-joining algorithm.
#'
#' The procedure results in a single tree per module/jackknifed module version (depending on whether the preservation statistics were calculated with or without jackknifing, see see the parameter \code{jackknife} in the function \code{\link{calculatePresStats}}). The tips of the trees represent the clones and the branch lengths represent the dissimilarity of module connectivity patterns between the networks of 2 clones. The trees are output as a list of \code{\link{phylo}} objects.
#'
#' In the next steps of the pipeline, statistics based on these trees can be used to identify conserved and diverged modules and pinpoint target genes within these modules that contribute the most to conservation/divergence (see \code{\link{calculateTreeStats}}, \code{\link{fitTreeStatsLm}}, \code{\link{findConservedDivergedModules}} and \code{\link{findConservedDivergedTargets}}).
#'
#' @param dist_list A named list containing the distance measures as data frames for all modules or jackknifed module versions, the output of \code{\link{convertPresToDist}}. Required columns for the data frames:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type ("orig" = original or "jk" = jackknifed, only needed if the distances were calculated with jackknifing).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type "jk" and nameOfRegulator_orig in case of module type "orig", only needed if the distances were calculated with jackknifing).}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type "orig", only needed if the distances were calculated with jackknifing).}
#' \item{clone1, clone2}{Character the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species \code{clone1} and \code{clone2} belong to, respectively.}
#' \item{dist}{Numeric, the distance measure to be used for the tree reconstruction.}
#' }
#' @param n_cores Number of cores.
#'
#' @return A named list containing the neighbor-joining trees as \code{\link{phylo}} objects for all modules or jackknifed module versions in the input \code{dist_list}. All trees contain a component \code{species} that specifies which species each tip belongs to and a component \code{info} that stores metadata of the module/jackknifed module version in a data frame format.
#' @export
#' @examples trees_jk <- reconstructTrees(dist_jk)
reconstructTrees <- function(dist_list, n_cores = 1L) {

  # check input data
  if (!inherits(dist_list, "list"))
    stop("The argument \"dist_list\" should be a named list.")

  if (is.null(names(dist_list)))
    stop("The argument \"dist_list\" should be a named list.")

  if (any(!sapply(dist_list, is.data.frame)))
    stop("One or more list elements of \"dist_list\" are not data frames.")

  if (!all(sapply(dist_list, function(df) {

    all(c("regulator", "module_size", "clone1", "clone2", "species1", "species2", "dist") %in% colnames(df))

  })))
    stop("All data frames in \"dist_list\" are expected to contain the columns \"regulator\", \"module_size\", \"clone1\", \"clone2\", \"species1\", \"species2\" and \"dist\".")

  if (!all(sapply(dist_list, function(df) {

    all(df$dist >= 0 & df$dist <= 1)

  })))
    warning("Not all distance measures in \"dist_list\" are between 0 and 1.")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") && !(inherits(n_cores, "numeric") && n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  id = NULL

  # convert each list element of dist_list into a tree, then create a list of these trees
  doParallel::registerDoParallel(n_cores)

  tree_list <- foreach::foreach(id = names(dist_list)) %dopar%

    {

      # reconstruct tree
      dist <- dist_list[[id]]
      tree <- getNjTree(dist)

      # add metadata
      tree$info <- dplyr::distinct(dist[, colnames(dist) %in% c("regulator", "module_size", "type", "id", "gene_removed")])

      tree

    }

  doParallel::stopImplicitCluster()

  names(tree_list) <- names(dist_list)

  tree_list

}


#' Convert preservation statistics to distances
#'
#' @description Converts preservation statistics of a single module to distance measures ranging from 0 to 1.
#' @param pres_stats Data frame of the preservation statistics for a single module or jackknifed module version, required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type (orig = original or jk = jackknifed, only needed if the preservation statistics were calculated with jackknifing).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig', only needed if the preservation statistics were calculated with jackknifing).}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type 'orig', only needed if the preservation statistics were calculated with jackknifing).}
#' \item{clone1, clone2}{Character the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species \code{clone1} and \code{clone2} belong to, respectively.}
#' \item{stat}{Numeric, the preservation statistic that is to be converted into a distance measure.}
#' }
#' @param min_stat Numeric, the minimum value of the preservation statistics that corresponds to a distance of 1.
#' @param max_stat Numeric, the maximum value of the preservation statistics that corresponds to a distance of 0.
#'
#' @return A data frame of distances containing 1 new column in addition to the relevant columns of \code{pres_stats}:
#' \describe{
#' \item{dist}{Numeric, the distance measure ranging from 0 to 1, calculated based on \code{stat}.}
#' }
#' @noRd
getDist <- function(pres_stats, min_stat, max_stat) {

  # convert the preservation statistics to distance measures
  dist <- pres_stats %>%
    dplyr::mutate(dist = (max_stat - .data[["stat"]]) / (max_stat - min_stat))

  # expand the data frame so that it contains all clone pairs (i.e. cloneA-cloneB and cloneB-cloneA are 2 separate entries with the same distance)
  dist_full <- dplyr::bind_rows(dist,
                   dist %>%
                     dplyr::rename(clone2_new = .data[["clone1"]], clone1_new = .data[["clone2"]], species2_new = .data[["species1"]], species1_new = .data[["species2"]]) %>%
                     dplyr::rename(clone1 = .data[["clone1_new"]], clone2 = .data[["clone2_new"]], species1 = .data[["species1_new"]], species2 = .data[["species2_new"]])) %>%
    dplyr::arrange(.data[["clone2"]]) %>%
    dplyr::arrange(.data[["clone1"]])

  dist_full[, colnames(dist_full) %in% c("regulator", "module_size", "type", "id", "gene_removed", "clone1", "clone2", "species1", "species2", "dist")]

}


#' Reconstruct trees
#'
#' @description Reconstructs a neighbor-joining tree based on pairwise distance measures between clones for a single module.
#' @param dist Data frame of the distance measures for a single module or jackknifed module version, required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of target genes assigned to a regulator.}
#' \item{type}{Character, module type (orig = original or jk = jackknifed, only needed if the distances were calculated with jackknifing).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig', only needed if the distances were calculated with jackknifing).}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type 'orig', only needed if the distances were calculated with jackknifing).}
#' \item{clone1, clone2}{Character the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species \code{clone1} and \code{clone2} belong to, respectively.}
#' \item{dist}{Numeric, the distance measure to be used for the tree reconstruction.}
#' }
#'
#' @return The neighbor-joining tree as a [phylo] obejct.
#' @noRd
getNjTree <- function(dist) {

  # convert the data frame of distance measures into a distance matrxi
  distMat <- dist %>%
    dplyr::mutate(clone1 = paste0(.data[["clone1"]], "_", .data[["species1"]]),
                  clone2 = paste0(.data[["clone2"]], "_", .data[["species2"]])) %>%
    dplyr::select(.data[["clone1"]], .data[["clone2"]], .data[["dist"]]) %>%
    tidyr::pivot_wider(names_from = .data[["clone2"]], values_from = .data[["dist"]], values_fill = 0) %>%
    tibble::column_to_rownames("clone1") %>%
    as.matrix()
  distMat <- distMat[, rownames(distMat)]

  # set seed and schedule the restoration of the old seed
  old_seed <- globalenv()$.Random.seed
  on.exit(suspendInterrupts({
    if (is.null(old_seed)) {
      rm(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      assign(".Random.seed", value = old_seed, envir = globalenv(), inherits = FALSE)
    }
  }), add = TRUE)
  set.seed(0)

  # reconstruct tree using the neighbor-joining algorithm
  tree <- ape::nj(distMat)

  tree$species <- factor(unlist(lapply(tree$tip.label, function(lab) {strsplit(lab, "\\_")[[1]][2]})), levels(dist$species1))
  tree$tip.label <- unlist(lapply(tree$tip.label, function(lab) {strsplit(lab, "\\_")[[1]][1]}))

  # get rid of negative branch lengths
  if(sum(tree$edge.length < 0) > 0) {

    id <- ifelse(is.null(dist$id), unique(as.character(dist$regulator)), unique(dist$id))

    warning(paste0("Negative branch lengths inferred for ", id, ", setting these to 0."))
    tree$edge.length[tree$edge.length < 0] <- 0

  }

  tree

}
