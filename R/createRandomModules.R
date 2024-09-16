#' Create random modules
#'
#' Creates a matching random module for each actual module.
#'
#' The function outputs a random module for each module in \code{pruned_modules}. The random modules have the same regulators and contain the same number of target genes as the original modules, but these target genes are randomly drawn from \code{network_genes}.
#'
#' In the next steps of the pipeline, the actual modules are compared to these random modules in terms of the statistics calculated to check whether the 2 groups of modules behave in general differently (see \code{\link{plotPresStatDistributions}}, \code{\link{plotPresStats}}, \code{\link{plotTreeStatDistributions}} and \code{\link{plotTreeStats}}) and to remove those individual actual modules that show too similar characteristics to the random modules (see \code{\link{filterModuleTrees}}).
#'
#' @param pruned_modules Data frame of pruned modules, required columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' }
#' @param network_genes Character vector of all genes in the network.
#' @param seed Integer, the seed to use for the random sampling (default: 42).
#' @return Data frame of the random modules with the following columns:
#'\describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{module_size}{Integer, the number of genes assigned to a regulator (only present if the column is also present in the input \code{pruned_modules}).}
#' \item{target}{Character, member gene of the regulator's random module.}
#' }
#' @export
#'
#' @examples random_modules <- createRandomModules(pruned_modules, genes)
createRandomModules <- function(pruned_modules, network_genes, seed = 42) {

  # check input data
  if (!is.data.frame(pruned_modules))
    stop("The argument \"pruned_modules\" should be a data frame.")

  if (any(!(c("regulator", "module_size", "target") %in% colnames(pruned_modules))))
    stop("The argument \"pruned_modules\" should contain the columns \"regulator\", \"module_size\" and \"target\".")

  if (!inherits(network_genes, "character"))
    stop("The argument \"network_genes\" should be a character vector.")

  if (!inherits(seed, "numeric"))
    stop("The argument \"seed\" should be a numeric.")

  # set seed and schedule the restoration of the old seed
  old_seed <- globalenv()$.Random.seed
  on.exit(suspendInterrupts({
    if (is.null(old_seed)) {
      rm(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      assign(".Random.seed", value = old_seed, envir = globalenv(), inherits = FALSE)
    }
  }), add = TRUE)
  set.seed(seed)

  # for each regulator, randomly sample as many targets from all network genes as there are actual targets
  pruned_modules %>%
    dplyr::select(dplyr::any_of(c("regulator", "target", "module_size"))) %>%
    dplyr::group_by(.data[["regulator"]]) %>%
    dplyr::mutate(target = sample(setdiff(network_genes, .data[["regulator"]]), dplyr::n())) %>%
    dplyr::ungroup()

}
