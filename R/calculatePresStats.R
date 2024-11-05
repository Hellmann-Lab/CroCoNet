#' Calculate preservation statistics
#'
#' Calculates one or more preservation statistics - correlation of adjacencies (cor_adj), correlation of regulator-target adjacencies (cor_adj_regulator), correlation of intramodular connectivities (cor_kIM) - for all modules and clone pairs with or without jackknifing.
#'
#' The function calculates one or more preservation statistics adapted from WGCNA (cor_adj, cor_adj_regulator and cor_kIM) for each module and each clone pair. All three statistics quantify how well the topology of a module is preserved between the networks of two clones. The statistic cor_adj is the correlation of all edge weights within the module in the network of clone1 VS in the network of clone2. The statistic cor_adj_regulator is the correlation of all edge weights between the regulator and its module members in the network of clone1 VS in the network of clone2. Finally, the statistic cor_kIM is the correlation of the intramodular connectivities per module member gene in the network of clone1 VS in the network of clone2.
#'
#' All statistics assume a joint module assignment (typically derived from the consensus network) but compare topological properties directly between the clonewise networks. In this approach, a module is always defined as the same set of genes, but the adjacencies/edge weights among these genes could differ from clone to clone; poorly preserved modules are expected to have many, while well-preserved modules are expected to have few such differences.
#'
#' If \code{jackknife} is set to TRUE (the default), the function creates all possible jackknifed versions of each input module by removing each target gene assigned to that module (the regulator is never excluded), then it calculates the preservation statistics for all of these jackknifed module versions in addition to the original module. This way, a confidence interval can be calculated for each module and statistic (see \code{\link{summarizeJackknifeStats}}). Later on in the pipeline, jackknifing can also provide information about which target genes within a conserved/diverged module are particularly responsible for the conservation/divergence (see \code{\link{findConservedDivergedTargets}}). If jackknifing is not desired, please set \code{jackknife} to FALSE. This will also substantially reduce running times.
#'
#' @param pruned_modules Data frame of pruned modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' }
#' @param network_list A named list of \code{\link{igraph}} objects containing the networks of all clones.
#' @param stats Character or character vector specifying which preservation statistics to calculate (one or more of "cor_adj", "cor_adj_regulator", "cor_kIM", default: c("cor_adj", "cor_adj_regulator", "cor_kIM")).
#' @param clone2species A data frame specifying which species each clone belongs to, required columns:
#' \describe{
#' \item{clone}{Character, name of the clone.}
#' \item{species}{Character, name of the species.}
#' }
#' If NULL (default), the output will contain no species information.
#' @param jackknife Logical specifying whether jackknifing should be performed or not (default: TRUE).
#' @param n_cores Integer, the number of cores (default: 1).
#' @param corr_method Character, the method for the calculation of correlation, one of "spearman", "pearson", "kendall" (default: "spearman").
#'
#' @return Data frame of the preservation statistics with the following columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the number of target genes assigned to a regulator (only present if the column is also present in \code{pruned_modules}).}
#' \item{type}{Character, module type ("orig" = original or "jk" = jackknifed, only present if parameter \code{jackknife} is set to TRUE).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type "jk" and nameOfRegulator_orig in case of module type "orig", only present if parameter \code{jackknife} is set to TRUE).}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type "orig", only present if parameter \code{jackknife} is set to TRUE).}
#' \item{clone1, clone2}{Character, the names of the clones compared.}
#' \item{species1, species2}{Character, the names of the species \code{clone1} and \code{clone2} belong to, respectively (only present if \code{clone2species} is not NULL).}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the preservation statistics specified by the parameter \code{stats}.}
#' }
#' @export
#' @examples pres_stats_jk <- calculatePresStats(pruned_modules, network_list, "cor_kIM", clone2species)
#' @references
#' Langfelder, P., Luo, R., Oldham, M. C., & Horvath, S. (2011). Is my network module preserved and reproducible? PLoS Computational Biology, 7(1), 1001057.
calculatePresStats <- function(pruned_modules, network_list, stats = c("cor_adj", "cor_adj_regulator", "cor_kIM"), clone2species = NULL, jackknife = TRUE, n_cores = 1L, corr_method = "spearman") {

  # check input data
  if (!is.data.frame(pruned_modules))
    stop("The argument \"pruned_modules\" should be a data frame.")

  if (any(!(c("regulator", "target") %in% colnames(pruned_modules))))
    stop("The argument \"pruned_modules\" should contain the columns \"regulator\" and \"target\".")

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

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") && !(inherits(n_cores, "numeric") && n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  if (is.null(corr_method) || !corr_method %in% c("spearman", "pearson", "kendall"))
    stop("The argument \"corr_method\" should be one of \"spearman\", \"pearson\", \"kendall\".")

  # avoid NSE notes in R CMD check
  clone = module = clone1 = clone2 = NULL

  # clone names
  clone_names <- names(network_list)

  # # if direction is present in the network objects, correlate direction*weight instead of weight
  # if ("direction" %in% names(igraph::edge_attr(network_list[[1]]))) {
  #
  #   network_list <- lapply(network_list, function(network) {
  #
  #     E(network)$weight = E(network)$weight * ifelse(E(network)$direction == "+", 1, -1)
  #     network
  #
  #   })
  #
  # }

  # convert igraphs to adjacency matrices
  adjMat_list <- lapply(network_list, function(network) {

    igraph::as_adjacency_matrix(network, attr = "weight")

  })

  # summarize module data frames (1 row - 1 module)
  pruned_modules_sum <- pruned_modules %>%
    dplyr::group_by(dplyr::across(dplyr::any_of(c("regulator", "module_size")))) %>%
    dplyr::summarise(genes = list(c(.data[["target"]], unique(as.character(.data[["regulator"]]))))) %>%
    dplyr::ungroup()

  # if jackknife is TRUE, add all possible jackknifed versions of each module
  if (jackknife) {

    pruned_modules_sum <- dplyr::bind_rows(pruned_modules_sum %>%
                                             # original modules
                                             dplyr::mutate(type = "orig",
                                                           id = paste0(.data[["regulator"]], "_orig")),
                                           pruned_modules_sum %>%
                                             # jackknifed module versions
                                             dplyr::mutate(jk = mapply(addJk, r = .data[["regulator"]], g = .data[["genes"]], SIMPLIFY = F)) %>%
                                             tidyr::unnest(.data[["jk"]]) %>%
                                             dplyr::mutate(genes = .data[["genes_jk"]]) %>%
                                             dplyr::select(dplyr::any_of(c("regulator", "module_size", "type", "id", "genes", "gene_removed")))) %>%
      dplyr::arrange(.data[["regulator"]])

  }

  if ("id" %in% colnames(pruned_modules_sum)) {

    id_column <- "id"

  } else {

    id_column <- "regulator"

  }

  module_names <- pruned_modules_sum[[id_column]]

  module_gene_list <- pruned_modules_sum$genes
  names(module_gene_list) <- module_names

  pruned_modules_sum$genes <- NULL

  if ("cor_adj" %in% stats) {

    doParallel::registerDoParallel(n_cores)

    adj_list <-

      foreach::foreach(clone = clone_names) %:%

      foreach::foreach(module = module_names,
                       .combine = c) %dopar% {

        adjMat <- adjMat_list[[clone]]
        genes <- module_gene_list[[module]]
        adj <- list(vectorize(adjMat[genes, genes]))
        names(adj) <- module
        adj

                       }

    names(adj_list) <- clone_names

    doParallel::stopImplicitCluster()

  }

  if ("cor_adj_regulator" %in% stats) {

    doParallel::registerDoParallel(n_cores)

    adj_regulator_list <-

      foreach::foreach(clone = clone_names) %:%

      foreach::foreach(module = module_names,
                       .combine = c) %dopar% {

                         adjMat <- adjMat_list[[clone]]
                         genes <- module_gene_list[[module]]
                         targets <- genes[-length(genes)]
                         regulator <- genes[length(genes)]
                         adj_regulator <- list(adjMat[targets, regulator])
                         names(adj_regulator) <- module
                         adj_regulator

                       }

    names(adj_regulator_list) <- clone_names

    doParallel::stopImplicitCluster()

  }

  if ("cor_kIM" %in% stats) {

    doParallel::registerDoParallel(n_cores)

    kIM_list <-

      foreach::foreach(clone = clone_names) %:%

      foreach::foreach(module = module_names,
                       .combine = c) %dopar% {

                         adjMat <- adjMat_list[[clone]]
                         genes <- module_gene_list[[module]]
                         kIM <- list(Matrix::rowSums(adjMat[genes, genes]))
                         names(kIM) <- module
                         kIM

                       }

    names(kIM_list) <- clone_names

    doParallel::stopImplicitCluster()

  }


  # calculate preservation statistics between all possible pairs of clones for all modules
  doParallel::registerDoParallel(n_cores)

  pres_stats <-

    foreach::foreach(clone1 = clone_names[1:(length(clone_names) - 1)],
                     .combine = rbind) %:%

    foreach::foreach(clone2 = clone_names[(which(clone_names == clone1) + 1):length(clone_names)],
                     .combine = rbind) %dopar%

    {

      # add clone info
      pres_stats_clone1_clone2 <- pruned_modules_sum %>%
          dplyr::mutate(clone1 = clone1,
                        clone2 = clone2)

      # calculate preservation statistics
      if ("cor_adj" %in% stats) {

        adj_clone1 <- adj_list[[clone1]]
        adj_clone2 <- adj_list[[clone2]]

        pres_stats_clone1_clone2 <- pres_stats_clone1_clone2 %>%
          dplyr::mutate(cor_adj = unname(sapply(.data[[id_column]], function(id) {

            stats::cor(adj_clone1[[id]],
                       adj_clone2[[id]],
                       method = corr_method)

          })))

      }

      if ("cor_adj_regulator" %in% stats) {

        adj_regulator_clone1 <- adj_regulator_list[[clone1]]
        adj_regulator_clone2 <- adj_regulator_list[[clone2]]

        pres_stats_clone1_clone2 <- pres_stats_clone1_clone2 %>%
          dplyr::mutate(cor_adj_regulator = unname(sapply(.data[[id_column]], function(id) {

            stats::cor(adj_regulator_clone1[[id]],
                       adj_regulator_clone2[[id]],
                       method = corr_method)

          })))

      }

      if ("cor_kIM" %in% stats) {

        kIM_clone1 <- kIM_list[[clone1]]
        kIM_clone2 <- kIM_list[[clone2]]

        pres_stats_clone1_clone2 <- pres_stats_clone1_clone2 %>%
          dplyr::mutate(cor_kIM = unname(sapply(.data[[id_column]], function(id) {

            stats::cor(kIM_clone1[[id]],
                       kIM_clone2[[id]],
                       method = corr_method)

          })))

      }

      pres_stats_clone1_clone2

    }

  doParallel::stopImplicitCluster()

  # if clone2species is provided, add species info
  if(!is.null(clone2species)) {

    suppressMessages(
      pres_stats <- pres_stats %>%
        dplyr::left_join(clone2species %>% dplyr::rename(clone1 = clone, species1 = species)) %>%
        dplyr::left_join(clone2species %>% dplyr::rename(clone2 = clone, species2 = species)) %>%
        dplyr::relocate(.data[["species1"]], .data[["species2"]], .after = "clone2")
    )

    # factor species to keep the order fixed in later functions
    species <- unique(c(pres_stats$species1, pres_stats$species2))
    pres_stats$species1 <- factor(pres_stats$species1, species)
    pres_stats$species2 <- factor(pres_stats$species2, species)

  }

  # factor clones to keep the order fixed in later functions
  pres_stats$clone1 <- factor(pres_stats$clone1, clone_names)
  pres_stats$clone2 <- factor(pres_stats$clone2, clone_names)

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
#' @return A data frame with the following columns:
#' \describe{
#' \item{type}{Character, module type ("orig" = original or "jk" = jackknifed, only present if parameter \code{jackknife} is set to TRUE).}
#' \item{id}{Character, the unique ID of the module version (format: nameOfRegulator_jk_nameOfGeneRemoved in case of module type "jk" and nameOfRegulator_orig in case of module type "orig", only present if parameter \code{jackknife} is set to TRUE).}
#' \item{gene_removed}{Character, the name of the gene removed by jackknifing (NA in case of module type "orig", only present if parameter \code{jackknife} is set to TRUE).}
#' \item{genes_jk}{A list containing a character vector, the remaining target genes of the transcriptional regulator.}
#' }
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
#' @param corr_method Character, the method for the calculation of correlation, one of "spearman", "pearson", "kendall".
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
#' @param corr_method Character, the method for the calculation of correlation, one of "spearman", "pearson", "kendall".
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
#' @param corr_method Character, the method for the calculation of correlation, one of "spearman", "pearson", "kendall".
#'
#' @return Numeric, the correlation coefficient of the intramodular connectivites between clone1 and clone2.
#' @noRd
calculateCorKIM <- function(g, mat1, mat2, corr_method = corr_method) {

  suppressMessages(stats::cor(Matrix::rowSums(mat1[g, g]),
                              Matrix::rowSums(mat2[g, g]),
                              method = corr_method))

}
