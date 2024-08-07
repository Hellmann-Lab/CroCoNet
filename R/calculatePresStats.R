#' Calculate preservation statistics
#'
#' Calculates one or more preservation statistics - correlation of adjacencies (cor_adj), correlation of regulator-target adjacencies (cor_adj_regulator), and/or correlation of intramodular connectivities (cor_kIM) - for all jackknifed versions of the modules and all clone pairs.
#'
#' The function calculates one or more preservation statistics adapted from Langfelder et al. 2011 (cor_adj, cor_adj_regulator and cor_kIM) for module and each clone pair. All three statistics quantify how well the module topology is preserved between the networks of two clones. The statistic cor_adj is the correlation of all edge weights within the module in the network of clone1 VS all edge weights within the module in the network of clone2. The statistic cor_adj_regulator is the correlation of all edge weights between the regulator and its module members in the network of clone1 VS all edge weights between the regulator and its module members in the network of clone2. Finally, the statistic cor_kIM is the correlation of the intramodular connectivities per module member gene in the network of clone1 VS the intramodular connectivities per module member gene in the network of clone2.
#'
#' All statistics assume a joint module assignment (typically derived from the consensus network) but compared topological properties directly between the clonewise networks. In this approach, a module is always defined as the same set of genes (nodes), but the interaction strengths (adjacencies) among these genes could differ from clone to clone; poorly preserved modules are expected to have many, while well-preserved modules are expected to have few such topological differences.
#'
#' If \code{jackknife} is set to TRUE, the function creates all possible jackknifed versions of each input module by removing each target gene assigned to that module (the regulator is never excluded), then it calculates the preservation statistics for all of these jackknifed module versions in addition to the original module. This way, a confidence interval can be calculated for each module and statistic (see \code{\link{summarizeJackknifeStats}}). Later on in the pipeline, jackknifing can also provide information about which target genes within a conserved/diverged module are particularly responsible for the conservation/divergence (see \code{\link{findConservedDivergedTargets}}). If jackknifing is not desired, please set \code{jackknife} to FALSE. This will also substantially reduce running times.
#'
#' @param pruned_modules Data frame of the pruned modules with columns 'regulator', 'target' and 'weight' containing the transcriptional regulators, their target genes and the edge weights between each regulator-target pair.
#' @param network_list A named list that contains the networks of each clone in an igraph format.
#' @param stats Character or character vector specifying which preservation statistics to calculate (one or more of "cor_adj", "cor_adj_regulator", and "cor_kIM").
#' @param clone2species A data frame with columns 'clone' and 'species' that specifies which species each clone belongs to. The names of clones should match the names of 'network_list'.
#' @param jackknife Logical specifying whether jackknifing should be performed or not (default: TRUE).
#' @param n_cores Number of cores.
#' @param corr_method The method for the calculation of correlation, one of "spearman" or "pearson".
#'
#' @return Data frame of the preservation statistics.
#' @export
#' @examples pres_stats_jk <- calculatePresStats(pruned_modules, network_list, "cor_kIM", clone2species)
calculatePresStats <- function(pruned_modules, network_list, stats = c("cor_adj", "cor_adj_regulator", "cor_kIM"), clone2species = NULL, jackknife = TRUE, n_cores = 1L, corr_method = "spearman") {

  if (!is.data.frame(pruned_modules))
    stop("The argument \"pruned_modules\" should be a data frame.")

  if (any(!(c("regulator", "module_size", "target") %in% colnames(pruned_modules))))
    stop("The argument \"pruned_modules\" should contain the columns \"regulator\", \"module_size\" and \"target\".")

  if (!inherits(network_list, "list"))
    stop("The argument \"network_list\" should be a named list.")

  if (is.null(names(network_list)))
    stop("The argument \"network_list\" should be a named list.")

  if (any(!sapply(network_list, function(net) {inherits(net, "igraph")})))
    stop("All elements of \"network_list\" should be of class \"igraph\".")

  if (is.null(stats) || !all(stats %in% c("cor_adj", "cor_adj_regulator", "cor_kIM")))
    stop("The argument \"stats\" should be one or more of \"cor_adj\", \"cor_adj_regulator\", \"cor_kIM\".")

  if (is.null(clone2species)) {

    warning("The argument \"clone2species\" is missing, no species information can be added to the preservation statistics. If you would like to proceed to the tree reconstruction step of the pipeline, please specify \"clone2species\" or add the columns \"species1\" and \"species2\" to the output manually.")

  } else {

    if (!is.data.frame(clone2species))
      stop("The argument \"clone2species\" should be a data frame.")

    if (any(!(c("clone", "species") %in% colnames(clone2species))))
      stop("The argument \"clone2species\" should contain the columns \"clone\" and \"species\".")

    if (any(!(names(network_list) %in% unique(clone2species$clone))))
      stop("The names of \"network_list\" should match the clone names in the column \"clone\" of \"clone2species\".")

  }

  if (!inherits(jackknife, "logical"))
    stop("The argument \"jackknife\" should be a logical value.")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") & !(inherits(n_cores, "numeric") & n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  if (is.null(corr_method) || !corr_method %in% c("spearman", "pearson", "kendall"))
    stop("The argument \"corr_method\" should be one of \"spearman\", \"pearson\", \"kendall\".")

  clone1 = clone2 = NULL

  # clone names
  clones <- names(network_list)

  # convert igraphs to adjacency matrices
  adjMat_list <- lapply(network_list, function(network) {

    igraph::as_adjacency_matrix(network, attr = "weight")

  })

  doParallel::registerDoParallel(n_cores)

  pres_stats <-

    foreach::foreach(clone1 = clones[1:(length(clones) - 1)],
                     .combine = rbind) %:%

    foreach::foreach(clone2 = clones[(which(clones == clone1) + 1):length(clones)],
                     .combine = rbind) %dopar%

    {

      message(paste0("Comparing ", clone1, " and ", clone2))

      adjMat1 <- adjMat_list[[clone1]]
      adjMat2 <- adjMat_list[[clone2]]

      pruned_modules_sum <- pruned_modules %>%
        dplyr::group_by(.data[["regulator"]], .data[["module_size"]]) %>%
        dplyr::summarise(genes = list(c(.data[["target"]], unique(as.character(.data[["regulator"]]))))) %>%
        dplyr::ungroup()

      if (jackknife) {

        pres_stats_clone1_clone2 <- dplyr::bind_rows(pruned_modules_sum %>%
                                                       # original modules
                                                       dplyr::mutate(type = "orig",
                                                                     id = paste0(.data[["regulator"]], "_orig")),
                                                     pruned_modules_sum %>%
                                                       # jackknifed module versions
                                                       dplyr::mutate(jk = mapply(addJk, r = .data[["regulator"]], g = .data[["genes"]], SIMPLIFY = F)) %>%
                                                       tidyr::unnest(.data[["jk"]]) %>%
                                                       dplyr::transmute(.data[["regulator"]],
                                                                        .data[["module_size"]],
                                                                        .data[["type"]],
                                                                        .data[["id"]],
                                                                        genes = .data[["genes_jk"]],
                                                                        .data[["gene_removed"]])) %>%
          dplyr::arrange(.data[["regulator"]]) %>%
          # add clone info
          dplyr::mutate(clone1 = clone1,
                        clone2 = clone2)

      } else {

        pres_stats_clone1_clone2 <- pruned_modules_sum %>%
          # add clone info
          dplyr::mutate(clone1 = clone1,
                        clone2 = clone2)

      }

      if(!is.null(clone2species)) {

        species1 <- clone2species$species[clone2species$clone == clone1]
        species2 <- clone2species$species[clone2species$clone == clone2]

        pres_stats_clone1_clone2 <- pres_stats_clone1_clone2 %>%
          dplyr::mutate(species1 = species1,
                        species2 = species2)

      }

      # calculate preservation statistics
      if ("cor_adj" %in% stats) {

        pres_stats_clone1_clone2 <- pres_stats_clone1_clone2 %>%
          dplyr::mutate(cor_adj = sapply(.data[["genes"]], calculateCorAdj, mat1 = adjMat1, mat2 = adjMat2, corr_method = corr_method))

      }

      if ("cor_adj_regulator" %in% stats) {

        pres_stats_clone1_clone2 <- pres_stats_clone1_clone2 %>%
          dplyr::mutate(cor_adj_regulator = sapply(.data[["genes"]], calculateCorAdjRegulator, mat1 = adjMat1, mat2 = adjMat2, corr_method = corr_method))

      }

      if ("cor_kIM" %in% stats) {

        pres_stats_clone1_clone2 <- pres_stats_clone1_clone2 %>%
          dplyr::mutate(cor_kIM = sapply(.data[["genes"]], calculateCorKIM, mat1 = adjMat1, mat2 = adjMat2, corr_method = corr_method))

      }

      pres_stats_clone1_clone2 %>%
        dplyr::select(-.data[["genes"]])

    }

  doParallel::stopImplicitCluster()

  pres_stats$clone1 <- factor(pres_stats$clone1, clones)
  pres_stats$clone2 <- factor(pres_stats$clone2, clones)

  if(!is.null(clone2species)) {

    species <- unique(c(pres_stats$species1, pres_stats$species2))
    pres_stats$species1 <- factor(pres_stats$species1, species)
    pres_stats$species2 <- factor(pres_stats$species2, species)

  }

  pres_stats

}


#' Vectorize matrix
#'
#' @description Vectorizes the non-redundant (i.e. upper-triangular) elements of a symmetrix matrix.
#' @param mat A symmetric matrix.
#'
#' @return A vector of non-redundnt (i.e. upper-triangular) matrix elements.
#' @noRd
vectorize <- function(mat) {

  mat[upper.tri(mat)]

}


#' Create jackknife modules
#'
#' @description Creates all possible jackknife versions of a module by removing each target gene (all jackknife versions still contain the regulator).
#' @param r Character, the name of the transcriptional regulator.
#' @param g Character vector, the names of the all genes in the module.
#'
#' @return A data frame with columns 'type', 'id', 'gene_removed' and 'target' containing the module type (orig = original or jk = jackknifed), the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type 'jk' and nameOfRegulator_orig in case of module type 'orig'), the name of the gene removed (NA in case of module type 'orig') and the remaining target genes.
#' @noRd
addJk <- function(r, g) {

  t <- g[g != r]

  data.frame(type = "jk",
             id = paste0(r, "_jk_", t),
             gene_removed = t) %>%
    dplyr::mutate(genes_jk = lapply(.data[["gene_removed"]], function(x){g[g != x]}))

}


#' Calculate correlation of adjacencies
#'
#' @param g Character vector, the genes of the module (regulator and targets) for which the correlation needs to be calculated.
#' @param mat1 Adjacency matrix of clone1.
#' @param mat2 Adjacency matrix of clone2.
#' @param corr_method The method for the calculation of correlation, one of "spearman" or "pearson".
#'
#' @return Numeric, the correlation coefficient of the intramodular adjacencies between clone1 and clone2.
#' @noRd
calculateCorAdj <- function(g, mat1, mat2, corr_method = corr_method) {

  suppressMessages(stats::cor(vectorize(mat1[g, g]),
                              vectorize(mat2[g, g]),
                              method = corr_method))

}


#' Calculate correlation of regulator-target adjacencies
#'
#' @param g Character vector, the genes of the module (regulator and targets) for which the correlation needs to be calculated.
#' @param mat1 Adjacency matrix of clone1.
#' @param mat2 Adjacency matrix of clone2.
#' @param corr_method The method for the calculation of correlation, one of "spearman" or "pearson".
#'
#' @return Numeric, the correlation coefficient of the intramodular regulator-target adjacencies between clone1 and clone2.
#' @noRd
calculateCorAdjRegulator <- function(g, mat1, mat2, corr_method = corr_method) {

  h <- g[length(g)]
  t <- g[-length(g)]

  suppressMessages(stats::cor(mat1[t, h],
                              mat2[t, h],
                              method = corr_method))

}


#' Calculate correlation of intramodular connectivities
#'
#' @param g Character vector, the genes of the module (regulator and targets) for which the correlation needs to be calculated.
#' @param mat1 Adjacency matrix of clone1.
#' @param mat2 Adjacency matrix of clone2.
#' @param corr_method The method for the calculation of correlation, one of "spearman" or "pearson".
#'
#' @return Numeric, the correlation coefficient of the intramodular connectivites between clone1 and clone2.
#' @noRd
calculateCorKIM <- function(g, mat1, mat2, corr_method = corr_method) {

  suppressMessages(stats::cor(Matrix::rowSums(mat1[g, g]),
                              Matrix::rowSums(mat2[g, g]),
                              method = corr_method))

}
