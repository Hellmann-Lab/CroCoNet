#' Prune modules
#'
#' Prunes the initial modules by keeping only the best targets of each transcriptional regulator. 3 methods are implemented to choose the best targets: 1) topN: takes a fixed number of targets per regulator with the highest regulator-target adjacencies. 2) UIK_adj: applies a dynamic stepwise pruning based on the regulator-target adjacencies, and 3) applies a dynamic stepwise pruning based on the regulator-target adjacencies and intramodular connectivities.
#' @param initial_modules Data frame of the initial modules with columns 'regulator', 'target' and 'weight' containing the transcriptional regulators, their target genes and the edge weights between each regulator-target pair.
#' @param network Consensus network in an igraph format.
#' @param method Pruning method, one of "UIK_adj", "UIK_adj_kIM" or "topN".
#' @param N An integer or a named integer vector specifying the desired pruned module size(s) in case of the method 'topN'.
#' @param min_module_size Integer, the lower threshold of module size in case of the methods 'UIK_adj' and UIK_adj_kIM'.
#' @param min_median_module_size Integer, the lower threshold of median module size in case of the methods 'UIK_adj' and UIK_adj_kIM'.
#' @param max_frac_modules_lost Numeric, the upper threshold of the fraction of removed modules in case of the methods 'UIK_adj' and UIK_adj_kIM'.
#' @param exponent Integer, the exponent the regulator-target adjacency and intramodular connectivity is raised to the power of during the cumulative sum curve calculation in case of the methods 'UIK_adj' and UIK_adj_kIM'.
#'
#' @return Data frame of the pruned modules.
#' @export
#'
#' @examples pruned_modules <- pruneModules(initial_modules, consensus_network, "UIK_adj_kIM")
pruneModules <- function(initial_modules, network = NULL, method = c("UIK_adj", "UIK_adj_kIM", "topN"), N = 50L, min_module_size = 20L, min_median_module_size = 20L, max_frac_modules_lost = 0.02, exponent = 1L) {

  if (method == "UIK_adj") {

    pruned_modules <- pruneModules_UIK_adj(initial_modules = initial_modules, min_module_size = min_module_size, min_median_module_size = min_median_module_size, max_frac_modules_lost = max_frac_modules_lost, exponent = exponent)

  } else if (method == "UIK_adj_kIM") {

    pruned_modules <- pruneModules_UIK_adj_kIM(initial_modules = initial_modules, network = network, min_module_size = min_module_size, min_median_module_size = min_median_module_size, max_frac_modules_lost = max_frac_modules_lost, exponent = exponent)

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
#' @param initial_modules Data frame of the initial modules with columns 'regulator', 'target' and 'weight' containing the transcriptional regulators, their target genes and the edge weights between each regulator-target pair.
#' @param min_module_size Integer, the lower treshold of module size. Modules with a smaller size than this threshold are removed after each pruning step.
#' @param min_median_module_size Integer, the lower threshold of median module size.
#' @param max_frac_modules_lost Numeric, the upper threshold of the fraction of removed modules.
#' @param exponent Integer, the exponent the regulator-target adjacency and intramodular connectivity is raised to the power of during the cumulative sum curve calculation.
#' The initial module members are filtered in successive steps based on their adjacency to the regulator. In each step, the cumulative sum curves of the regulator-target adjacencies are calculated per module, the knee point of the curve is identified using the Unit Invariant Knee (UIK) method, then only the targets that rank higher than the knee point are kept. If 'exponent' is set higher than 1, the adjacencies are raised to the equivalent power when calculating the cumulative sum curves and its knee point. The modules containing less target genes than 'min_module_size' are removed after each pruning step. The steps continue until the median module size becomes as small as possible without falling below 'min_median_module_size' OR until the fraction of removed modules becomes higher than 'max_frac_modules_lost'.
#' @return Data frame of the pruned modules.
#' @export
#'
#' @examples pruned_modules_UIK_adj <- pruneModules_UIK_adj(initial_modules)
pruneModules_UIK_adj <- function(initial_modules, min_module_size = 20L, min_median_module_size = 20L, max_frac_modules_lost = 0.02,  exponent = 1L) {

  if (!is.data.frame(initial_modules))
    stop("The argument \"initial_modules\" should be a data frame.")

  if (any(!c("regulator", "target", "weight") %in% colnames(initial_modules)))
    stop("The argument \"initial_modules\" should contain the columns \"regulator\", \"target\" and \"weight\".")

  if (length(min_module_size) != 1 || (!inherits(min_module_size, "integer") & !(inherits(min_module_size, "numeric") & min_module_size == round(min_module_size))) || min_module_size < 1)
    stop("The argument \"min_module_size\" should be a positive integer.")

  if (length(min_median_module_size) != 1 || (!inherits(min_median_module_size, "integer") & !(inherits(min_median_module_size, "numeric") & min_median_module_size == round(min_median_module_size))) || min_median_module_size < 1)
    stop("The argument \"min_median_module_size\" should be a positive integer.")

  if (!inherits(max_frac_modules_lost, "numeric") || length(max_frac_modules_lost) != 1 || max_frac_modules_lost < 0 || max_frac_modules_lost >= 1)
    stop("The argument \"max_frac_modules_lost\" should be a numeric value between 0 and 1.")

  if (length(exponent) != 1 || (!inherits(exponent, "integer") & !(inherits(exponent, "numeric") & exponent == round(exponent))) || exponent < 1)
    stop("The argument \"exponent\" should be a positive integer.")

  module_data <- list(initial_modules)

  repeat{

    intermediate_modules <- filterBasedOnAdj(module_data[[1]], exponent)

    intermediate_modules_filt <- removeSmallModules(intermediate_modules, min_module_size)

    median_module_size <- getMedianModuleSize(intermediate_modules)

    if ((getNumberOfModules(intermediate_modules_filt) / getNumberOfModules(intermediate_modules) < 1 - max_frac_modules_lost) || (median_module_size < min_median_module_size)) break

    message(paste0("Step ", length(module_data), ": filtering targets based on their adjacencies to the regulator"))

    message(paste0("Median module size after filtering: ", median_module_size))

    module_data <- c(list(intermediate_modules_filt), module_data)

  }

  module_data[[1]] %>%
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::arrange(.data[["target"]], .by_group = TRUE) %>%
    dplyr::ungroup()

}


#' Prune modules based on the regulator-target adjacencies and intramodular connectivities using dynamic filtering
#'
#' Prunes the initial modules by applying a dynamic stepwise pruning based on the regulator-target adjacencies and intramodular connectivties.
#' @param initial_modules Data frame of the initial modules with columns 'regulator', 'target' and 'weight' containing the transcriptional regulators, their target genes and the edge weights between each regulator-target pair.
#' @param network Consensus network in an igraph format.
#' @param min_module_size Integer, the lower treshold of module size. Modules with a smaller size than this threshold are removed after each pruning step.
#' @param min_median_module_size Integer, the lower threshold of median module size.
#' @param max_frac_modules_lost Numeric, the upper threshold of the fraction of removed modules.
#' @param exponent Integer, the exponent the regulator-target adjacency and intramodular connectivity is raised to the power of during the cumulative sum curve calculation.
#' The initial module members are filtered in successive steps based on their adjacency to the regulator and their intramodular connectivitiy alternately. In each step, the cumulative sum curves based on one of these two characteristics are calculated per module, the knee point of the curve is identified using the Unit Invariant Knee (UIK) method, then only the targets that rank higher than the knee point are kept. The intramodular connectivity is recalculated in each relevant filtering step based on the then-current module assignment. If 'exponent' is set higher than 1, the regulator-target adjacencies and intramodular connectivities are raised to the equivalent power when calculating the cumulative sum curves and its knee point. The modules containing less target genes than 'min_module_size' are removed after each pruning step. The steps continue until the median module size becomes as small as possible without falling below 'min_median_module_size' OR until the fraction of removed modules becomes higher than 'max_frac_modules_lost'.
#' @return Data frame of the pruned modules.
#' @export
#'
#' @examples pruned_modules_UIK_adj_kIM <- pruneModules_UIK_adj_kIM(initial_modules, consensus_network)
pruneModules_UIK_adj_kIM <- function(initial_modules, network, min_module_size = 20L, min_median_module_size = 20L, max_frac_modules_lost = 0.02,  exponent = 1L) {

  if (!is.data.frame(initial_modules))
    stop("The argument \"initial_modules\" should be a data frame.")

  if (any(!c("regulator", "target", "weight") %in% colnames(initial_modules)))
    stop("The argument \"initial_modules\" should contain the columns \"regulator\", \"target\" and \"weight\".")

  if (!inherits(network, "igraph"))
    stop("The argument \"network\" should be of class \"igraph\".")

  genes <- sort(unique(c(as.character(initial_modules$regulator), initial_modules$target)))
  genes_not_found <- genes[!genes %in% V(network)$name]

  if (length(genes_not_found) > 0)
    stop(paste0("The following genes in the initial modules were not found in the provided network: ", paste(genes_not_found, ", "), "."))

  if (length(min_module_size) != 1 || (!inherits(min_module_size, "integer") & !(inherits(min_module_size, "numeric") & min_module_size == round(min_module_size))) || min_module_size < 1)
    stop("The argument \"min_module_size\" should be a positive integer.")

  if (length(min_median_module_size) != 1 || (!inherits(min_median_module_size, "integer") & !(inherits(min_median_module_size, "numeric") & min_median_module_size == round(min_median_module_size))) || min_median_module_size < 1)
    stop("The argument \"min_median_module_size\" should be a positive integer.")

  if (!inherits(max_frac_modules_lost, "numeric") || length(max_frac_modules_lost) != 1 || max_frac_modules_lost < 0 || max_frac_modules_lost >= 1)
    stop("The argument \"max_frac_modules_lost\" should be a numeric value between 0 and 1.")

  if (length(exponent) != 1 || (!inherits(exponent, "integer") & !(inherits(exponent, "numeric") & exponent == round(exponent))) || exponent < 1)
    stop("The argument \"exponent\" should be a positive integer.")

  adj_mat <- igraph::as_adjacency_matrix(network, attr = "weight")

  module_data <- list(initial_modules)

  repeat{

    intermediate_modules <- filterBasedOnAdj(module_data[[1]], exponent)

    intermediate_modules_filt <- removeSmallModules(intermediate_modules, min_module_size)

    median_module_size <- getMedianModuleSize(intermediate_modules)

    if ((getNumberOfModules(intermediate_modules_filt) / getNumberOfModules(intermediate_modules) < 1 - max_frac_modules_lost) || (median_module_size < min_median_module_size)) break

    message(paste0("Step ", length(module_data), ": filtering targets based on their adjacencies to the regulator"))

    message(paste0("Median module size after filtering: ", median_module_size))

    module_data <- c(list(intermediate_modules_filt), module_data)



    intermediate_modules <- filterBasedOnkIM(module_data[[1]], adj_mat, exponent)

    intermediate_modules_filt <- removeSmallModules(intermediate_modules, min_module_size)

    median_module_size <- getMedianModuleSize(intermediate_modules)

    if ((getNumberOfModules(intermediate_modules_filt) / getNumberOfModules(intermediate_modules) < 1 - max_frac_modules_lost) || (median_module_size < min_median_module_size)) break

    message(paste0("Step ", length(module_data), ": filtering targets based on their intramodular connectivities"))

    message(paste0("Median module size after filtering: ", median_module_size))

    module_data <- c(list(intermediate_modules_filt), module_data)

  }

  module_data[[1]] %>%
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::arrange(.data[["target"]], .by_group = TRUE) %>%
    dplyr::ungroup()

}


#' Prune modules based on the regulator-target adjacencies by keeping a fixed number of top targets
#'
#' Prunes the initial modules by keeping a fixed number of targets per transcriptional regulator with the highest regulator-target adjacencies
#' @param initial_modules Data frame of the initial modules with columns 'regulator', 'target' and 'weight' containing the transcriptional regulators, their target genes and the edge weights between each regulator-target pair.
#' @param N Either an integer specifying a single desired pruned module size for all modules or a named integer vector specifying the desired pruned module size for each regulator.
#'
#' @return Data frame of the pruned modules.
#' @export
#'
#' @examples pruned_modules_top50 <- pruneModules_topN(initial_modules, 50)
pruneModules_topN <- function(initial_modules, N = 50L) {

  if (!is.data.frame(initial_modules))
    stop("The argument \"initial_modules\" should be a data frame.")

  if (any(!c("regulator", "target", "weight") %in% colnames(initial_modules)))
    stop("The argument \"initial_modules\" should contain the columns \"regulator\", \"target\" and \"weight\".")

  if (length(N) == 1) {

    if ((!inherits(N, "integer") & !(inherits(N, "numeric") & N == round(N))) || N < 1)
      stop("The argument \"N\" should be a positive integer or a named vector of positive integers.")

    pruned_modules <- initial_modules %>%
      dplyr::group_by(.data[["regulator"]]) %>%
      dplyr::slice_max(order_by = .data[["weight"]], n = N) %>%
      dplyr::arrange(.data[["target"]], .by_group = TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(module_size = N)

  } else {

    if ((!inherits(N, "integer") & !(inherits(N, "numeric") & all(N == round(N)))) || any(N < 1) || is.null(names(N)))
      stop("The argument \"N\" should be a positive integer or a named vector of positive integers.")

    if (length(N) != length(unique(initial_modules$regulator)))
      stop("If \"N\" is a named vector of integers, it should have the same length as there are regulators in the column \"regulator\" of \"initial_modules\".")

    if (any(!names(N) %in% initial_modules$regulator))
      stop("The names of \"N\" should be regulators present in the column \"regulator\" of \"initial_modules\".")

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
#' @param modules Data frame of the initial or intermediate modules with columns 'regulator', 'target' and 'weight' containing the transcriptional regulators, their target genes and the edge weights between each regulator-target pair.
#' @param exponent Integer, the exponent the regulator-target adjacency is raised to the power of during the cumulative sum curve calculation.
#'
#' @return Data frame of the filtered modules.
#' @noRd
filterBasedOnAdj <- function(modules, exponent = 1L) {

  modules %>%
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::arrange(dplyr::desc(.data[["weight"]]), .by_group = TRUE) %>%
    dplyr::mutate(rank_IS = 1:length(.data[["target"]]),
                  cumsum_IS = cumsum(.data[["weight"]]^exponent),
                  rank_cut_IS = inflection::uik(.data[["rank_IS"]], .data[["cumsum_IS"]])) %>%
    dplyr::filter(.data[["rank_IS"]] <= .data[["rank_cut_IS"]]) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data[["rank_IS"]], -.data[["cumsum_IS"]], -.data[["rank_cut_IS"]])

}


#' Filter modules based on the intramodular connectivities using knee-point detection
#'
#' @param modules Data frame of the initial or intermediate modules with columns 'regulator', 'target' and 'weight' containing the transcriptional regulators, their target genes and the edge weights between each regulator-target pair.
#' @param adj_mat Adjacency matrix of the consensus network in a sparse or dense matrix format.
#' @param exponent Integer, the exponent the regulator-target adjacency is raised to the power of during the cumulative sum curve calculation.
#'
#' @return Data frame of the filtered modules.
#' @noRd
filterBasedOnkIM <- function(modules, adj_mat, exponent = 1L) {

  modules %>%
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::mutate(kIM = Matrix::rowSums(adj_mat[.data[["target"]], c(.data[["target"]], unique(as.character(.data[["regulator"]])))])) %>%
    dplyr::arrange(dplyr::desc(.data[["kIM"]]), .by_group = TRUE) %>%
    dplyr::mutate(rank_kIM = 1:length(.data[["target"]]),
                  cumsum_kIM = cumsum(.data[["kIM"]]^exponent),
                  rank_cut_kIM = inflection::uik(.data[["rank_kIM"]], .data[["cumsum_kIM"]])) %>%
    dplyr::filter(.data[["rank_kIM"]] <= .data[["rank_cut_kIM"]]) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data[["kIM"]], -.data[["rank_kIM"]], -.data[["cumsum_kIM"]], -.data[["rank_cut_kIM"]])

}


#' Calculate the median module size
#'
#' @param modules Data frame of the initial, intermediate or pruned modules with columns 'regulator', 'target' and 'weight' containing the transcriptional regulators, their target genes and the edge weights between each regulator-target pair.
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
#' @param modules Data frame of the initial, intermediate or pruned modules with columns 'regulator', 'target' and 'weight' containing the transcriptional regulator, their target genes and the edge weights between each regulator-target pair.
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
#' @param modules Data frame of the initial, intermediate or pruned modules with columns 'regulator', 'target' and 'weight' containing the transcriptional regulator, their target genes and the edge weights between each regulator-target pair.
#'
#' @return Integer, the number of modules.
#' @noRd
getNumberOfModules <- function(modules) {

  length(unique(modules$regulator))

}
