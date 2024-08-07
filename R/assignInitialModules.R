#' Assign initial modules
#'
#' Assigns initial modules by taking the top targets of each provided transcriptional regulator.
#'
#' In the CroCoNet approach, modules are centered around transcriptional regulators and assigned in 2 main steps: 1) large initial modules are created by selecting a fixed number of target genes per regulator (performed by this function), and 2) the initial modules are pruned to keep only the best targets of each regulator (performed by \code{\link{pruneModules}}).
#'
#' The module assignment is recommended to be done based on the consensus network. The regulators that provide the starting point of the module assignment can be selected based on prior biological knowledge, or the combination of prior biological knowledge and the data (see also \code{\link{getRegulators}}). As the default, all transcription factors with at least 1 annotated motif in the JASPAR 2024 vertebrate core collection are used.
#'
#' The function creates as many modules as there are regulators, each containing the regulator and its \code{N} best target genes. When choosing the best targets, the genes are ranked based on how strongly they are connected to the regulator (regulator-target adjecency/edge weight). The parameter \code{N} should be greater than or equal to the minimum number of targets per regulator across all regulators.
#'
#' @param network \code{\link{igraph}} object, the consensus network across all species and clones.
#' @param regulators Character vector of transcriptional regulators.
#' @param N Integer, the initial module size, i.e. the number of target genes to keep for each regulator.
#'
#' @return A data frame of initial modules with 3 core columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, member gene of the regulator's initial module.}
#' \item{weight}{Numeric, consensus edge weight/adjacency, the weighted average of clonewise adjacencies.}
#' }
#' If the input \code{network} had more edge attributes, those appear as additional columns in the data frame.
#' @export
#'
#' @examples
#' initial_modules <- assignInitialModules(consensus_network, regulators, N = 250)
assignInitialModules <- function(network, regulators = intersect(V(network)$name, jaspar_core_TRs), N = 3000L) {

  if (!inherits(network, "igraph"))
    stop("The argument \"network\" should be of class \"igraph\".")

  if (!inherits(regulators, "character"))
    stop("The argument \"regulators\" should be a character vector.")

  if (length(N) != 1 || (!inherits(N, "integer") & !(inherits(N, "numeric") & N == round(N))) | N < 1)
    stop("The argument \"N\" should be a positive integer.")

  if (all(!(regulators %in% V(network)$name)))
    stop("None of the regulators are present in the network.")

  df <- igraph::as_data_frame(network)

  all_targets <- dplyr::bind_rows(df,
                                  df %>% dplyr::rename(from2 = .data[["from"]], to2 = .data[["to"]]) %>% dplyr::rename(from = .data[["to2"]], to = .data[["from2"]])) %>%
    dplyr::filter(.data[["from"]] %in% regulators) %>%
    dplyr::rename(regulator = .data[["from"]], target = .data[["to"]])

  if (min(table(all_targets$regulator)) < N) {

    warning("One or more regulators have less targets than the specified initial module size, consider decreasing N.")

  }

  all_targets %>%
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::slice_max(order_by = .data[["weight"]], n = N) %>%
    dplyr::arrange(.data[["target"]], .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(regulator = factor(.data[["regulator"]], regulators)) %>%
    dplyr::arrange(.data[["regulator"]])

}
