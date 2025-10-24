#' Calculate preservation statistics
#'
#' Calculates one or both of two preservation statistics - correlation of adjacencies (cor_adj) and correlation of intramodular connectivities (cor_kIM) - for all modules and replicate pairs with or without jackknifing.
#'
#' The function calculates two preservation statistics adapted from WGCNA (cor_adj and cor_kIM) for each module and each replicate pair. Both statistics quantify how well the topology of a module is preserved between the networks of two replicates. The statistic cor_adj is the correlation of all edge weights within the module in the network of replicate1 VS in the network of replicate2, while the statistic cor_kIM is the correlation of the intramodular connectivities per module member gene in the network of replicate1 VS in the network of replicate2.
#'
#' All statistics assume a joint module assignment (typically derived from the consensus network) but compare topological properties directly between the replicate-wise networks. In this approach, a module is always defined as the same set of genes, but the adjacencies/edge weights among these genes could differ from replicate to replicate; poorly preserved modules are expected to have many, while well-preserved modules are expected to have few such differences.
#'
#' If \code{jackknife} is set to TRUE (the default), the function creates all possible jackknifed versions of each input module by removing each target gene assigned to that module (the regulator is never excluded), then it calculates the preservation statistics for all of these jackknifed module versions in addition to the original module. This way, a confidence interval can be calculated for each module and statistic (see \code{\link{summarizeJackknifeStats}}). Later on in the pipeline, jackknifing can also provide information about which target genes within a conserved/diverged module are particularly responsible for the conservation/divergence (see \code{\link{findConservedDivergedTargets}}). If jackknifing is not desired, please set \code{jackknife} to FALSE. This will also substantially reduce running times.
#'
#' @param pruned_modules Data frame of pruned modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' }
#' @param network_list A named list of \code{\link{igraph}} objects containing the networks of all replicates.
#' @param stats Character or character vector specifying which preservation statistics to calculate (one or more of "cor_adj", "cor_kIM", default: c("cor_adj", "cor_kIM")).
#' @param replicate2species A data frame specifying which species each replicate belongs to, required columns:
#' \describe{
#' \item{replicate}{Character, name of the replicate.}
#' \item{species}{Character, name of the species.}
#' }
#' If NULL (default), the output will contain no species information.
#' @param jackknife Logical specifying whether jackknifing should be performed or not (default: TRUE).
#' @param signed Logical indicating whether the networks in \code{network_list} are signed (default: FALSE, see also \code{\link{normalizeEdgeWeights}}). If set to FALSE and \code{network_list} contains the edge attribute "direction", the edge weights in the network with direction "-" are negated for the calculation of "cor_adj".
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
#' \item{replicate1, replicate2}{Character, the names of the replicates compared.}
#' \item{species1, species2}{Character, the names of the species \code{replicate1} and \code{replicate2} belong to, respectively (only present if \code{replicate2species} is not NULL).}
#' \item{\{\{nameOfStat\}\}}{Numeric, one or more columns containing the preservation statistics specified by the parameter \code{stats}.}
#' }
#' @export
#' @examples pres_stats_jk <- calculatePresStats(pruned_modules, network_list, "cor_kIM", replicate2species)
#' @references
#' Langfelder, P., Luo, R., Oldham, M. C., & Horvath, S. (2011). Is my network module preserved and reproducible? PLoS Computational Biology, 7(1), 1001057.
calculatePresStats <- function(pruned_modules, network_list, stats = c("cor_adj", "cor_kIM"), replicate2species = NULL, jackknife = TRUE, signed = FALSE, n_cores = 1L, corr_method = "spearman") {

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

  if (is.null(stats) || !all(stats %in% c("cor_adj", "cor_kIM")))
    stop("The argument \"stats\" should be one or more of \"cor_adj\", \"cor_kIM\".")

  if (is.null(replicate2species)) {

    warning("The argument \"replicate2species\" is missing, no species information can be added to the preservation statistics. If you would like to proceed to the tree reconstruction step of the pipeline, please specify \"replicate2species\" or add the columns \"species1\" and \"species2\" to the output manually.")

  } else {

    if (!is.data.frame(replicate2species))
      stop("The argument \"replicate2species\" should be a data frame.")

    if (any(!(c("replicate", "species") %in% colnames(replicate2species))))
      stop("The argument \"replicate2species\" should contain the columns \"replicate\" and \"species\".")

    if (any(!(names(network_list) %in% unique(replicate2species$replicate))))
      stop("The names of \"network_list\" should match the replicate names in the column \"replicate\" of \"replicate2species\".")

  }

  if (!inherits(jackknife, "logical") || length(jackknife) != 1)
    stop("The argument \"jackknife\" should be a logical value.")

  if (!inherits(signed, "logical") || length(signed) != 1)
    stop("The argument \"signed\" should be a logical value.")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") && !(inherits(n_cores, "numeric") && n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  if (is.null(corr_method) || !corr_method %in% c("spearman", "pearson", "kendall"))
    stop("The argument \"corr_method\" should be one of \"spearman\", \"pearson\", \"kendall\".")

  # avoid NSE notes in R CMD check
  replicate = module = replicate1 = replicate2 = NULL

  # replicate names
  replicate_names <- names(network_list)

  # if direction is present in the network objects, correlate direction*weight instead of weight
  if ("direction" %in% names(igraph::edge_attr(network_list[[1]])) && !signed && "cor_adj" %in% stats) {

    message("Calculating cor_adj using weight*direction")

    network_list <- lapply(network_list, function(network) {

      directions <- igraph::E(network)$direction

      if (any(!directions %in% c("+", "-") & !is.na(directions)))
        stop("The column \"direction\" of the networks in \"network_list\" should contain only the following values: \"+\", \"-\", NA.")

      igraph::E(network)$weight = igraph::E(network)$weight * ifelse(directions == "+" | is.na(directions), 1, -1)
      network

    })

  }

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

      foreach::foreach(replicate = replicate_names) %:%

      foreach::foreach(module = module_names,
                       .combine = c) %dopar% {

        adjMat <- adjMat_list[[replicate]]
        genes <- module_gene_list[[module]]
        adj <- list(vectorize(adjMat[genes, genes]))
        names(adj) <- module
        adj

                       }

    names(adj_list) <- replicate_names

    doParallel::stopImplicitCluster()

  }

  if ("cor_kIM" %in% stats) {

    doParallel::registerDoParallel(n_cores)

    kIM_list <-

      foreach::foreach(replicate = replicate_names) %:%

      foreach::foreach(module = module_names,
                       .combine = c) %dopar% {

                         adjMat <- abs(adjMat_list[[replicate]])
                         genes <- module_gene_list[[module]]
                         kIM <- list(Matrix::rowSums(adjMat[genes, genes]))
                         names(kIM) <- module
                         kIM

                       }

    names(kIM_list) <- replicate_names

    doParallel::stopImplicitCluster()

  }


  # calculate preservation statistics between all possible pairs of replicates for all modules
  doParallel::registerDoParallel(n_cores)

  pres_stats <-

    foreach::foreach(replicate1 = replicate_names[1:(length(replicate_names) - 1)],
                     .combine = rbind) %:%

    foreach::foreach(replicate2 = replicate_names[(which(replicate_names == replicate1) + 1):length(replicate_names)],
                     .combine = rbind) %dopar%

    {

      # add replicate info
      pres_stats_replicate1_replicate2 <- pruned_modules_sum %>%
          dplyr::mutate(replicate1 = replicate1,
                        replicate2 = replicate2)

      # calculate preservation statistics
      if ("cor_adj" %in% stats) {

        adj_replicate1 <- adj_list[[replicate1]]
        adj_replicate2 <- adj_list[[replicate2]]

        pres_stats_replicate1_replicate2 <- pres_stats_replicate1_replicate2 %>%
          dplyr::mutate(cor_adj = unname(sapply(.data[[id_column]], function(id) {

            stats::cor(adj_replicate1[[id]],
                       adj_replicate2[[id]],
                       method = corr_method)

          })))

      }

      if ("cor_kIM" %in% stats) {

        kIM_replicate1 <- kIM_list[[replicate1]]
        kIM_replicate2 <- kIM_list[[replicate2]]

        pres_stats_replicate1_replicate2 <- pres_stats_replicate1_replicate2 %>%
          dplyr::mutate(cor_kIM = unname(sapply(.data[[id_column]], function(id) {

            stats::cor(kIM_replicate1[[id]],
                       kIM_replicate2[[id]],
                       method = corr_method)

          })))

      }

      pres_stats_replicate1_replicate2

    }

  doParallel::stopImplicitCluster()

  # if replicate2species is provided, add species info
  if(!is.null(replicate2species)) {

    suppressMessages(
      pres_stats <- pres_stats %>%
        dplyr::left_join(replicate2species %>% dplyr::rename(replicate1 = replicate, species1 = species)) %>%
        dplyr::left_join(replicate2species %>% dplyr::rename(replicate2 = replicate, species2 = species)) %>%
        dplyr::relocate(.data[["species1"]], .data[["species2"]], .after = "replicate2")
    )

    # factor species to keep the order fixed in later functions
    species <- unique(c(pres_stats$species1, pres_stats$species2))
    pres_stats$species1 <- factor(pres_stats$species1, species)
    pres_stats$species2 <- factor(pres_stats$species2, species)

  }

  # factor replicates to keep the order fixed in later functions
  pres_stats$replicate1 <- factor(pres_stats$replicate1, replicate_names)
  pres_stats$replicate2 <- factor(pres_stats$replicate2, replicate_names)

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
#' @param mat1 Adjacency matrix of replicate1.
#' @param mat2 Adjacency matrix of replicate2.
#' @param corr_method Character, the method for the calculation of correlation, one of "spearman", "pearson", "kendall".
#'
#' @return Numeric, the correlation coefficient of the intramodular adjacencies between replicate1 and replicate2.
#' @noRd
calculateCorAdj <- function(g, mat1, mat2, corr_method = corr_method) {

  suppressMessages(stats::cor(vectorize(mat1[g, g]),
                              vectorize(mat2[g, g]),
                              method = corr_method))

}


#' Calculate correlation of intramodular connectivities
#'
#' @param g Character vector, the genes of the module (regulator and targets) for which the correlation needs to be calculated.
#' @param mat1 Adjacency matrix of replicate1.
#' @param mat2 Adjacency matrix of replicate2.
#' @param corr_method Character, the method for the calculation of correlation, one of "spearman", "pearson", "kendall".
#'
#' @return Numeric, the correlation coefficient of the intramodular connectivites between replicate1 and replicate2.
#' @noRd
calculateCorKIM <- function(g, mat1, mat2, corr_method = corr_method) {

  suppressMessages(stats::cor(Matrix::rowSums(mat1[g, g]),
                              Matrix::rowSums(mat2[g, g]),
                              method = corr_method))

}
