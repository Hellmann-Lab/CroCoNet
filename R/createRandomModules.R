#' Create random modules
#'
#' Creates a matching random module for each actual module.
#' @param pruned_modules Data frame of the pruned modules with columns 'regulator', 'target' and 'weight' containing the transcriptional regulators, their target genes and the edge weights between each regulator-target pair.
#' @param network_genes Character vector of all genes in the network.
#' Outputs a random module for each module in 'pruned_modules'. The random modules have the same regulators and contain the same number of target genes as the original modules, but these target genes are randomly drawn from 'network_genes'.
#' @param seed Integer, the seed to use for the random sampling.
#' @return Data frame of the random modules.
#' @export
#'
#' @examples random_modules <- createRandomModules(pruned_modules, genes)
createRandomModules <- function(pruned_modules, network_genes, seed = 0) {

  if (!is.data.frame(pruned_modules))
    stop("The argument \"pruned_modules\" should be a data frame.")

  if (any(!(c("regulator", "module_size", "target") %in% colnames(pruned_modules))))
    stop("The argument \"pruned_modules\" should contain the columns \"regulator\", \"module_size\" and \"target\".")

  if (!inherits(network_genes, "character"))
    stop("The argument \"network_genes\" should be a character vector.")

  # set seed and schedule the restoration of the old seed
  restoreOldSeed()
  set.seed(seed)

  # for each regulator, randomly sample as many targets from all network genes as there are actual targets
  pruned_modules %>%
    dplyr::select(.data[["regulator"]], .data[["target"]], .data[["module_size"]]) %>%
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::mutate(target = sample(setdiff(network_genes, .data[["regulator"]]), dplyr::n())) %>%
    dplyr::ungroup()

}
