#' Prune modules
#'
#' Prunes the initial modules by keeping only the best targets of each transcriptional regulator.
#'
#' 3 methods are implemented to choose the best targets:
#' \itemize{
#'  \item{topN: Takes a fixed number of targets per regulator with the highest regulator-target adjacencies (for details see \code{\link{pruneModules_topN}}).}
#'  \item{UIK_adj: Applies a dynamic stepwise pruning based on the regulator-target adjacencies (for details see \code{\link{pruneModules_UIK_adj}}).}
#'  \item{UIK_adj_kIM: Applies a dynamic stepwise pruning based on the regulator-target adjacencies and intramodular connectivities (for details see \code{\link{pruneModules_UIK_adj_kIM}}).}
#'  }
#'
#' @param initial_modules Data frame of initial modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, member gene of the regulator's initial module.}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#' @param method Character, the pruning method, one of "UIK_adj", "UIK_adj_kIM", "topN".
#' @param consensus_network \code{\link{igraph}} object, the consensus network across all species and clones.
#' @param N An integer or a named integer vector specifying the desired pruned module size(s) in case of the method "topN" (default: 50).
#' @param min_module_size Integer, the lower threshold of module size in case of the methods "UIK_adj" and "UIK_adj_kIM" (default: 20).
#' @param max_frac_modules_lost Numeric, the upper threshold of the fraction of removed modules in case of the methods "UIK_adj" and "UIK_adj_kIM" (default: 0.02).
#' @param exponent Integer, the exponent the regulator-target adjacency and intramodular connectivity is raised to the power of during the cumulative sum curve calculation in case of the methods "UIK_adj" and "UIK_adj_kIM" (default: 1, i.e. the regulator-target adjacencies and intramodular connectivities stay unchanged).
#'
#' @return Data frame of the pruned modules with the following columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of genes assigned to a regulator.}
#' \item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#' Additional columns present in \code{initial_modules} will also be preserved in \code{pruned_modules}.
#' @export
#'
#' @examples
#' pruned_modules <- pruneModules(initial_modules, "topN", N = 30)
#' pruned_modules <- pruneModules(initial_modules, "UIK_adj")
#' pruned_modules <- pruneModules(initial_modules, "UIK_adj_kIM", consensus_network)
pruneModules <- function(initial_modules, method = c("UIK_adj", "UIK_adj_kIM", "topN"), consensus_network = NULL, min_module_size = 20L, max_frac_modules_lost = 0.02, exponent = 1L, N = 50L) {

  if (method == "UIK_adj") {

    pruned_modules <- pruneModules_UIK_adj(initial_modules = initial_modules, min_module_size = min_module_size, max_frac_modules_lost = max_frac_modules_lost, exponent = exponent)

  } else if (method == "UIK_adj_kIM") {

    pruned_modules <- pruneModules_UIK_adj_kIM(initial_modules = initial_modules, consensus_network = consensus_network, min_module_size = min_module_size, max_frac_modules_lost = max_frac_modules_lost, exponent = exponent)

  } else if (method == "topN") {

    pruned_modules <- pruneModules_topN(initial_modules = initial_modules, N = N)

  } else {

    stop("Please set method to one of 'UIK_adj, 'UIK_adj_kIM', 'topN'.")

  }

  return(pruned_modules)

}


#' Prune modules based on the regulator-target adjacencies using dynamic filtering
#'
#' Prunes the initial modules by applying a dynamic stepwise pruning based on the regulator-target adjacencies.
#'
#' For each module, the initial module members are filtered in successive steps based on their regulator-target edge weight/adjacency, which quantifies how strongly a target is connected to the regulator. In each step, the cumulative sum curve of the regulator-target adjacencies is calculated per module, the knee point of the curve is identified using the Unit Invariant Knee (UIK) method (see \code{\link{uik}}), then only the targets that rank higher than the knee point are kept.
#'
#' The modules containing less target genes than \code{min_module_size} are removed after each pruning step. The steps continue until the fraction of removed modules becomes higher than \code{max_frac_modules_lost}.
#'
#' It is recommended to set \code{min_module_size} to at least 20, because the correlation-based preservation statistics in the next steps might be coupled with high uncertainty for modules smaller than this (see \code{\link{calculatePresStats}}). If you would like to keep all modules for further analysis, please set \code{max_frac_modules_lost} to 0.
#'
#' If \code{exponent} is set higher than 1, the adjacencies are raised to the equivalent power when calculating the cumulative sum curves and their knee points.
#'
#' While setting the parameter \code{min_module_size} prevents the modules from becoming too small, the exact number of target genes per regulator does not have to be pre-defined, in line with the notion that different regulators can have an effect on different numbers of genes. There is also no hard cutoff applied to the regulator-target adjacencies, but by using knee point detection the target genes are filtered in a data-driven way.
#'
#' The modules are allowed to overlap, and in addition to having its own module, a regulator can be assigned to another regulator's module as well, in line with the notion that genes can be multifunctional and gene regulation can be combinatorial.
#'
#' @param initial_modules Data frame of initial modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, member gene of the regulator's initial module.}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#' @param min_module_size Integer, the lower threshold of module size. Modules with a smaller size than this threshold are removed after each pruning step (default: 20).
#' @param max_frac_modules_lost Numeric, the upper threshold of the fraction of removed modules (default: 0.02).
#' @param exponent Integer, the exponent the regulator-target adjacency and intramodular connectivity is raised to the power of during the cumulative sum curve calculation (default: 1, i.e. the regulator-target adjacencies and intramodular connectivities stay unchanged).
#'
#' @return Data frame of the pruned modules with the following columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of genes assigned to a regulator.}
#' \item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#' Additional columns present in \code{initial_modules} will also be preserved in \code{pruned_modules}.
#' @export
#'
#' @examples pruned_modules_UIK_adj <- pruneModules_UIK_adj(initial_modules)
#' @family methods to prune modules
#' @references
#' Christopoulos, D. (2016). Introducing Unit Invariant Knee (UIK) As an Objective Choice for Elbow Point in Multivariate Data Analysis Techniques. SSRN Electronic Journal. https://doi.org/10.2139/SSRN.3043076
pruneModules_UIK_adj <- function(initial_modules, min_module_size = 20L, max_frac_modules_lost = 0.02,  exponent = 1L) {

  # check input data
  if (!is.data.frame(initial_modules))
    stop("The argument \"initial_modules\" should be a data frame.")

  if (any(!c("regulator", "target", "weight") %in% colnames(initial_modules)))
    stop("The argument \"initial_modules\" should contain the columns \"regulator\", \"target\" and \"weight\".")

  if (length(min_module_size) != 1 || (!inherits(min_module_size, "integer") && !(inherits(min_module_size, "numeric") && min_module_size == round(min_module_size))) || min_module_size < 1)
    stop("The argument \"min_module_size\" should be a positive integer.")

  if (!inherits(max_frac_modules_lost, "numeric") || length(max_frac_modules_lost) != 1 || max_frac_modules_lost < 0 || max_frac_modules_lost >= 1)
    stop("The argument \"max_frac_modules_lost\" should be a numeric value between 0 and 1.")

  if (length(exponent) != 1 || (!inherits(exponent, "integer") && !(inherits(exponent, "numeric") && exponent == round(exponent))) || exponent < 1)
    stop("The argument \"exponent\" should be a positive integer.")

  # initialize list of module data frames
  module_data <- list(initial_modules)

  # repeat steps iteratively until too many modules fall below the specified minimum module size
  repeat{

    # filter targets based on their regulator-target adjacency
    intermediate_modules <- filterBasedOnAdj(module_data[[1]], exponent)

    # remove modules that fall below the specified minimum module size after the filtering
    intermediate_modules_filt <- removeSmallModules(intermediate_modules, min_module_size)

    # if too many modules had to be removed, break and revert to the previous version
    if (getNumberOfModules(intermediate_modules_filt) / getNumberOfModules(intermediate_modules) < 1 - max_frac_modules_lost) break

    # get median module size
    median_module_size <- getMedianModuleSize(intermediate_modules_filt)

    # print stuff
    message(paste0("Step ", length(module_data), ": filtering targets based on their adjacencies to the regulator"))
    message(paste0("Median module size after filtering: ", median_module_size))

    # add new module data frame to the list
    module_data <- c(list(intermediate_modules_filt), module_data)

  }

  # print message if any modules were lost due to low module size
  modules_lost <- setdiff(unique(initial_modules$regulator), unique(module_data[[1]]$regulator))
  if (length(modules_lost) > 0)
    message(paste0("The following modules were removed due to too low module size after pruning: ", paste(modules_lost, collapse = ", "), ". If you would like to keep them, consider lowering 'min_module_size'."))

  # take the last module data frame that still passed the size criterion
  module_data[[1]] %>%
    # sort targets
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::arrange(.data[["target"]], .by_group = TRUE) %>%
    dplyr::ungroup()

}


#' Prune modules based on the regulator-target adjacencies and intramodular connectivities using dynamic filtering
#'
#' Prunes the initial modules by applying a dynamic stepwise pruning based on the regulator-target adjacencies and intramodular connectivities.
#'
#' For each module, the initial module members are filtered in successive steps based on 2 metrics alternately: 1) the regulator-target edge weight/adjacency, which quantifies how strongly a target is connected to the regulator, and 2) the intramodular connectivity, which quantifies how strongly a target is connected to all other genes in the module. In each step, the cumulative sum curve based on one of these two characteristics is calculated per module, the knee point of the curve is identified using the Unit Invariant Knee (UIK) method (see \code{\link{uik}}), then only the targets that rank higher than the knee point are kept. The first step is based on the regulator-target adjacencies, the second step is based on the intramodular connectivities, the third step is again based on the regulator-target adjacencies and so on. The intramodular connectivity is recalculated in each relevant filtering step based on the then-current module assignment.
#'
#' The modules containing less target genes than \code{min_module_size} are removed after each pruning step. The steps continue until the fraction of removed modules becomes higher than \code{max_frac_modules_lost}.
#'
#' It is recommended to set \code{min_module_size} to at least 20, because the correlation-based preservation statistics in the next steps might be coupled with high uncertainty for modules smaller than this (see \code{\link{calculatePresStats}}). If you would like to keep all modules for further analysis, please set \code{max_frac_modules_lost} to 0.
#'
#' If \code{exponent} is set higher than 1, the regulator-target adjacencies and intramodular connectivities are raised to the equivalent power when calculating the cumulative sum curves and their knee points.
#'
#' Pruning based on intramodular connectivities in addition to the regulator-target adjacencies ensures that the chosen targets co-vary not just with the main regulator but also with the rest of the module. These intramodular connections between targets can carry important information about combinatorial regulation, feedback loops and co-functionality.
#'
#' While setting the parameter \code{min_module_size} prevents the modules from becoming too small, the exact number of target genes per regulator does not have to be pre-defined, in line with the notion that different regulators can have an effect on different numbers of genes. There are also no hard cutoffs applied to the regulator-target adjacencies or intramodular connectivities, but by using knee point detection the target genes are filtered in a data-driven way.
#'
#' The modules are allowed to overlap, and in addition to having its own module, a regulator can be assigned to another regulator's module as well, in line with the notion that genes can be multifunctional and gene regulation can be combinatorial.
#'
#' @param initial_modules Data frame of initial modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, member gene of the regulator's initial module.}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#' @param consensus_network \code{\link{igraph}} object, the consensus network across all species and clones.
#' @param min_module_size Integer, the lower treshold of module size. Modules with a smaller size than this threshold are removed after each pruning step (default: 20).
#' @param max_frac_modules_lost Numeric, the upper threshold of the fraction of removed modules (default: 0.02).
#' @param exponent Integer, the exponent the regulator-target adjacency and intramodular connectivity is raised to the power of during the cumulative sum curve calculation (default: 1, i.e. the regulator-target adjacencies and intramodular connectivities stay unchanged).
#'
#' @return Data frame of the pruned modules with the following columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of genes assigned to a regulator.}
#' \item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#' Additional columns present in \code{initial_modules} will also be preserved in \code{pruned_modules}.
#' @export
#'
#' @examples pruned_modules_UIK_adj_kIM <- pruneModules_UIK_adj_kIM(initial_modules, consensus_network)
#' @family methods to prune modules
#' @references
#' Christopoulos, D. (2016). Introducing Unit Invariant Knee (UIK) As an Objective Choice for Elbow Point in Multivariate Data Analysis Techniques. SSRN Electronic Journal. https://doi.org/10.2139/SSRN.3043076
pruneModules_UIK_adj_kIM <- function(initial_modules, consensus_network, min_module_size = 20L, max_frac_modules_lost = 0.02, exponent = 1L) {

  # check input data
  if (!is.data.frame(initial_modules))
    stop("The argument \"initial_modules\" should be a data frame.")

  if (any(!c("regulator", "target", "weight") %in% colnames(initial_modules)))
    stop("The argument \"initial_modules\" should contain the columns \"regulator\", \"target\" and \"weight\".")

  if (!inherits(consensus_network, "igraph"))
    stop("The argument \"consensus_network\" should be of class \"igraph\".")

  genes <- sort(unique(c(as.character(initial_modules$regulator), initial_modules$target)))
  genes_not_found <- genes[!genes %in% V(consensus_network)$name]

  if (length(genes_not_found) > 0)
    stop(paste0("The following genes in the initial modules were not found in the provided network: ", paste(genes_not_found, ", "), "."))

  if (length(min_module_size) != 1 || (!inherits(min_module_size, "integer") && !(inherits(min_module_size, "numeric") && min_module_size == round(min_module_size))) || min_module_size < 1)
    stop("The argument \"min_module_size\" should be a positive integer.")

  if (!inherits(max_frac_modules_lost, "numeric") || length(max_frac_modules_lost) != 1 || max_frac_modules_lost < 0 || max_frac_modules_lost >= 1)
    stop("The argument \"max_frac_modules_lost\" should be a numeric value between 0 and 1.")

  if (length(exponent) != 1 || (!inherits(exponent, "integer") && !(inherits(exponent, "numeric") && exponent == round(exponent))) || exponent < 1)
    stop("The argument \"exponent\" should be a positive integer.")

  # convert igraph to adjacency matrix
  adj_mat <- igraph::as_adjacency_matrix(consensus_network, attr = "weight")

  # initialize list of module data frames
  module_data <- list(initial_modules)

  # repeat steps iteratively until too many modules fall below the specified minimum module size
  repeat{

    # filter targets based on their regulator-target adjacency
    intermediate_modules <- filterBasedOnAdj(module_data[[1]], exponent)

    # remove modules that fall below the specified minimum module size after the filtering
    intermediate_modules_filt <- removeSmallModules(intermediate_modules, min_module_size)

    # if too many modules had to be removed, break and revert to the previous version
    if (getNumberOfModules(intermediate_modules_filt) / getNumberOfModules(intermediate_modules) < 1 - max_frac_modules_lost) break

    # get median module size
    median_module_size <- getMedianModuleSize(intermediate_modules_filt)

    # print stuff
    message(paste0("Step ", length(module_data), ": filtering targets based on their adjacencies to the regulator"))
    message(paste0("Median module size after filtering: ", median_module_size))

    # add new module data frame to the list
    module_data <- c(list(intermediate_modules_filt), module_data)


    # filter targets based on their intramodular connectivity
    intermediate_modules <- filterBasedOnkIM(module_data[[1]], adj_mat, exponent)

    # remove modules that fall below the specified minimum module size after the filtering
    intermediate_modules_filt <- removeSmallModules(intermediate_modules, min_module_size)

    # if too many modules had to be removed, break and revert to the previous version
    if (getNumberOfModules(intermediate_modules_filt) / getNumberOfModules(intermediate_modules) < 1 - max_frac_modules_lost) break

    # get median module size
    median_module_size <- getMedianModuleSize(intermediate_modules_filt)

    # print stuff
    message(paste0("Step ", length(module_data), ": filtering targets based on their intramodular connectivities"))
    message(paste0("Median module size after filtering: ", median_module_size))

    # add new module data frame to the list
    module_data <- c(list(intermediate_modules_filt), module_data)

  }

  # print message if any modules were lost due to low module size
  modules_lost <- setdiff(unique(initial_modules$regulator), unique(module_data[[1]]$regulator))
  if (length(modules_lost) > 0)
    message(paste0("The following modules were removed due to too low module size after pruning: ", paste(modules_lost, collapse = ", "), ". If you would like to keep them, consider lowering 'min_module_size'."))

  # take the last module data frame that still passed the size criterion
  module_data[[1]] %>%
    # sort targets
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::arrange(.data[["target"]], .by_group = TRUE) %>%
    dplyr::ungroup()

}


#' Prune modules based on the regulator-target adjacencies by keeping a fixed number of top targets
#'
#' Prunes the initial modules by keeping a fixed number of targets per transcriptional regulator with the highest regulator-target adjacencies.
#'
#' Each pruned module output by the function contains the regulator and its \code{N} best target genes. When choosing the best targets, the genes are ranked based on how strongly they are connected to the regulator (regulator-target edge weight/adjacency).
#'
#' Based on prior biological knowledge, \code{N} can be set to a different value for different regulators, however, in most cases it will just be the same desired module size for all modules. Fixing the sizes of all modules to the same number is a simple but widespread approach.
#'
#' The modules are allowed to overlap, and in addition to having its own module, a regulator can be assigned to another regulator's module as well, in line with the notion that genes can be multifunctional and gene regulation can be combinatorial.
#'
#' @param initial_modules Data frame of initial modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, member gene of the regulator's initial module.}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#' @param N Either an integer specifying a single desired pruned module size for all modules or a named integer vector specifying the desired pruned module size for each regulator (default: 50).
#'
#' @return Data frame of the pruned modules with the following columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the numer of genes assigned to a regulator.}
#' \item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#' Additional columns present in \code{initial_modules} will also be preserved in \code{pruned_modules}.
#' @export
#'
#' @examples pruned_modules_top50 <- pruneModules_topN(initial_modules, 50)
#' @family methods to prune modules
pruneModules_topN <- function(initial_modules, N = 50L) {

  # check input data
  if (!is.data.frame(initial_modules))
    stop("The argument \"initial_modules\" should be a data frame.")

  if (any(!c("regulator", "target", "weight") %in% colnames(initial_modules)))
    stop("The argument \"initial_modules\" should contain the columns \"regulator\", \"target\" and \"weight\".")

  if (length(N) == 1) {

    if ((!inherits(N, "integer") && !(inherits(N, "numeric") && N == round(N))) || N < 1)
      stop("The argument \"N\" should be a positive integer or a named vector of positive integers.")

    # if N is a single integer, keep the N best targets based on regulator-target adjacency in each module
    pruned_modules <- initial_modules %>%
      dplyr::group_by(.data[["regulator"]]) %>%
      dplyr::slice_max(order_by = .data[["weight"]], n = N) %>%
      dplyr::arrange(.data[["target"]], .by_group = TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(module_size = N)

  } else {

    if ((!inherits(N, "integer") && !(inherits(N, "numeric") && all(N == round(N)))) || any(N < 1) || is.null(names(N)))
      stop("The argument \"N\" should be a positive integer or a named vector of positive integers.")

    if (length(N) != length(unique(initial_modules$regulator)))
      stop("If \"N\" is a named vector of integers, it should have the same length as there are regulators in the column \"regulator\" of \"initial_modules\".")

    if (any(!names(N) %in% initial_modules$regulator))
      stop("The names of \"N\" should be regulators present in the column \"regulator\" of \"initial_modules\".")

    # if N is a named vector of integers, keep the individually specified number of best targets based on regulator-target adjacency in each module
    pruned_modules <- initial_modules %>%
      dplyr::group_by(.data[["regulator"]]) %>%
      dplyr::mutate(module_size = N[unique(as.character(.data[["regulator"]]))]) %>%
      dplyr::filter(.data[["weight"]] %in% sort(.data[["weight"]], decreasing = TRUE)[1:unique(.data[["module_size"]])]) %>%
      dplyr::ungroup()

  }

  pruned_modules

}


#' Filter modules based on the regulator-target adjacencies using knee-point detection
#'
#' Calculates the cumulative sum curves of the regulator-target adjacencies per module, identifies the knee point of the curve using the Unit Invariant Knee (UIK) method (see \code{\link{uik}}), and keeps only the target genes that rank higher than the knee point.
#'
#' @param modules Data frame of the initial or intermediate modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, member gene of the regulator's initial/intermediate module.}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#' @param exponent Integer, the exponent the regulator-target adjacency is raised to the power of during the cumulative sum curve calculation.
#'
#' @return Data frame of the filtered modules.
#' @noRd
filterBasedOnAdj <- function(modules, exponent = 1L) {

  modules %>%
    # order and rank targets by regulator-target adjacency
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::arrange(dplyr::desc(.data[["weight"]]), .by_group = TRUE) %>%
    # calculate the cumulative sum of regulator-target adjacencies at each rank and determine the knee point of the cumulative sum curve
    dplyr::mutate(rank_IS = 1:length(.data[["target"]]),
                  cumsum_IS = cumsum(.data[["weight"]]^exponent),
                  rank_cut_IS = inflection::uik(.data[["rank_IS"]], .data[["cumsum_IS"]])) %>%
    # remove targets that rank lower than the knee point
    dplyr::filter(.data[["rank_IS"]] <= .data[["rank_cut_IS"]]) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data[["rank_IS"]], -.data[["cumsum_IS"]], -.data[["rank_cut_IS"]])

}


#' Filter modules based on the intramodular connectivities using knee-point detection
#'
#' Calculates the cumulative sum curves of the intramodular connectivities per module, identifies the knee point of the curve using the Unit Invariant Knee (UIK) method (see \code{\link{uik}}), and keeps only the target genes that rank higher than the knee point.
#'
#' @param modules Data frame of the initial or intermediate modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, member gene of the regulator's initial/intermediate module.}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#' @param adj_mat Adjacency matrix of the consensus network in a sparse or dense matrix format.
#' @param exponent Integer, the exponent the regulator-target adjacency is raised to the power of during the cumulative sum curve calculation.
#'
#' @return Data frame of the filtered modules.
#' @noRd
filterBasedOnkIM <- function(modules, adj_mat, exponent = 1L) {

  modules %>%
    # calculate the intramodular connectivity of each target
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::mutate(kIM = Matrix::rowSums(adj_mat[.data[["target"]], c(.data[["target"]], unique(as.character(.data[["regulator"]])))])) %>%
    # order and rank targets by intarmodular connectivity
    dplyr::arrange(dplyr::desc(.data[["kIM"]]), .by_group = TRUE) %>%
    # calculate the cumulative sum of intramodular connectivities at each rank and determine the knee point of the cumulative sum curve
    dplyr::mutate(rank_kIM = 1:length(.data[["target"]]),
                  cumsum_kIM = cumsum(.data[["kIM"]]^exponent),
                  rank_cut_kIM = inflection::uik(.data[["rank_kIM"]], .data[["cumsum_kIM"]])) %>%
    # remove targets that rank lower than the knee point
    dplyr::filter(.data[["rank_kIM"]] <= .data[["rank_cut_kIM"]]) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data[["kIM"]], -.data[["rank_kIM"]], -.data[["cumsum_kIM"]], -.data[["rank_cut_kIM"]])

}


#' Calculate the median module size
#'
#' Calculates the median module size (the median number of target genes assigned to a regulator) across all modules.
#'
#' @param modules Data frame of the initial, intermediate or pruned modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, member gene of the regulator's initial/intermediate/pruned module.}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#'
#' @return Numeric, the median size of the modules.
#' @noRd
getMedianModuleSize <- function(modules) {

  modules %>%
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::summarize(module_size = length(.data[["target"]])) %>%
    dplyr::pull(.data[["module_size"]]) %>%
    stats::median()

}


#' Remove too small modules
#'
#' Removes the modules that have a lower module size (fewer target gene assigned to a regulator) than \code{min_module_size}.
#'
#' @param modules Data frame of the initial, intermediate or pruned modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, member gene of the regulator's initial/intermediate/pruned module.}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#' @param min_module_size Integer, the size (number of genes) a module needs to reach to be kept.
#'
#' @return Data frame of filtered modules.
#' @noRd
removeSmallModules <- function(modules, min_module_size) {

  modules %>%
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::mutate(module_size = length(.data[["target"]])) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data[["module_size"]] >= min_module_size)

}


#' Get the number of modules
#'
#' @param modules Data frame of the initial, intermediate or pruned modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, member gene of the regulator's initial/intermediate/pruned module.}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted mean of clonewise edge weights.}
#' }
#'
#' @return Integer, the number of modules.
#' @noRd
getNumberOfModules <- function(modules) {

  length(unique(modules$regulator))

}
