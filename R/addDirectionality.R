#' Add directionality of interactions
#'
#' Determines the directionality of each network connection based on a modified Spearman's correlation between the 2 genes' expression profiles (positive expression correlation - activating interaction, negative expression correlation - repressing interaction).
#'
#' If the networks were inferred using a method that cannot distinguish positive and negative regulatory interactions, it might be useful to add this information for lines of analysis where it makes sense to separate the activated and repressed target genes of a regulator (e.g. for the calculation of eigengenes).
#'
#' The calculation relies on the approximate version of the Spearman's rho, significance testing and blocking implemented by [scran::correlatePairs]. The results are summarized as 3 new edge attributes in the igraph objects: rho (approximate Spearman's correlation coefficient), p.adj (BH-corrected approximate p-value) and direction (+ or -).
#' @param network_list A named list of \code{\link{igraph}} objects containing the networks per clone and potentially the consensus network. The list element containing the consensus network should have the name 'consensus'.
#' @param sce_list A named list of \code{\link{SingleCellExperiment-class}} objectd containing the expression data (raw counts, normaized counts and metadata) for all network genes per clone.
#' @param n_cores Integer, the number of cores (default: 1).
#' @param assay Character, the name of the assay in the SingleCellExperiment objects that should be used for the calculation of gene-gene correlations (default: "logcounts").
#'
#' @return A named list of \code{\link{igraph}} objects after adding the information about the direction of interactions. In addition to the edge attributes in the input \code{network_list}, it contains 3 new attributes in each network:
#' \describe{
#' \item{rho}{Numeric, the approximate Spearman's correlation coefficient of the 2 genes' expression profiles that form the edge.}
#' \item{p.adj}{Numeric, BH-corrected approximate p-value of rho.}
#' \item{direction}{Character, the direction of the interaction between the 2 genes that form the edge ("+" or "-").}
#' }
#'
#' @export
#'
#' @examples
#' network_list_scaled_filt_withCons_withDir <- addDirectionality(network_list_scaled_filt_withCons,
#'                                                                sce_list)
addDirectionality <- function(network_list, sce_list, n_cores = 1L, assay = "logcounts") {

  # check input data
  if (!inherits(network_list, "list"))
    stop("The argument \"network_list\" should be a named list.")

  if (is.null(names(network_list)))
    stop("The argument \"network_list\" should be a named list.")

  if (!inherits(sce_list, "list"))
    stop("The argument \"sce_list\" should be a named list.")

  if (is.null(names(sce_list)))
    stop("The argument \"sce_list\" should be a named list.")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") & !(inherits(n_cores, "numeric") & n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  if ((!inherits(assay, "character") || length(assay) != 1))
    stop("The argument \"assay\" should be a string.")

  if (length(sce_list) != length(network_list[names(network_list) != "consensus"]))
    stop("The arguments \"network_list\" and \"sce_list\" should be of the same length (except when \"network_list\" contains the list element \"consensus\", in this case it should be 1 element longer than \"sce_list\").")

  if (!all(sort(names(sce_list)) == sort(setdiff(names(network_list), "consensus"))))
    stop("The names of \"network_list\" should match the names of \"sce_list\" (except for the list element \"consensus\" that does not have to have an equivalent element in \"sce_list\").")

  if (any(!sapply(sce_list, function(sce) {inherits(sce, "SingleCellExperiment")})))
    stop("All elements of \"sce_list\" should be of class \"SingleCellExperiment\".")

  if (any(!sapply(network_list, function(net) {inherits(net, "igraph")})))
    stop("All elements of \"network_list\" should be of class \"igraph\".")

  if (any(!sapply(sce_list, function(sce) {assay %in% names(SummarizedExperiment::assays(sce))})))
    stop("One or more SingleCellExperiment objects in \"sce_list\" do not contain the specified assay.")

  # set to NULL due to NSE notes in R CMD check
  . = from = to = rho = p.adj = FDR = direction = p.value = clone_name = NULL

  # convert igraphs to data tables
  dt_list <- convertToDT(network_list)

  # get directionality for each clonewise network
  if (all(names(sce_list) %in% names(dt_list))) {

    clone_names <- names(sce_list)

    dt_list_with_dir <- foreach::foreach(clone_name = clone_names) %do% {

       # expression matrix for the clone
       expr <- SummarizedExperiment::assay(sce_list[[clone_name]], assay)

       # connections in the clonewise network
       pairings <- as.matrix(dt_list[[clone_name]][, .(from, to)])

       # set seed and schedule the restoration of the old seed
       restoreOldSeed()
       set.seed(0)

       # calculate correlations
       corr = suppressWarnings(data.table::as.data.table(scran::correlatePairs(expr, pairings = pairings, BPPARAM = BiocParallel::MulticoreParam(workers = n_cores))))

       # add correlations to the network data table
       dt_list[[clone_name]][corr, on = c(from = "gene1", to = "gene2")][
         , `:=` (p.adj = data.table::fifelse(is.na(FDR), 0, FDR), direction = data.table::fifelse(rho > 0, "+", "-"), p.value = NULL, FDR = NULL)]

    }
    names(dt_list_with_dir) <- clone_names

  } else {

    dt_list_with_dir <- list()

  }

  # get directionality for the consensus network
  if ("consensus" %in% names(dt_list)) {

    # expression matric including all clones (correlations will be calculated separately within each clone and then only the statistics will be combined -> per clone normalisation is ok, no neeed to normalise again everything together)
    expr_combined <- rlist::list.cbind(lapply(sce_list, function(sce) { SummarizedExperiment::assay(sce,assay) }))

    # connections in the consensus network
    pairings_combined <- as.matrix(dt_list[["consensus"]][, .(from, to)])

    # block info
    block <- dplyr::bind_rows(lapply(sce_list, function(sce) {as.data.frame(SummarizedExperiment::colData(sce))})) %>%
      dplyr::select(.data[["clone"]]) %>%
      tibble::rownames_to_column("cell") %>%
      tibble::deframe()

    # set seed and schedule the restoration of the old seed
    old_seed <- globalenv()$.Random.seed
    on.exit(assign(".Random.seed", value = old_seed, envir = globalenv(), inherits = FALSE), add = TRUE)
    set.seed(0)

    # calculate correlations
    corr <- suppressWarnings(data.table::as.data.table(scran::correlatePairs(expr_combined, pairings = pairings_combined, block = block, BPPARAM = BiocParallel::MulticoreParam(workers = n_cores))))

    # add correlations to the network data table
    dt_list_with_dir[["consensus"]] <- dt_list[["consensus"]][corr, on = c(from = "gene1", to = "gene2")][
      , `:=` (p.adj = data.table::fifelse(is.na(FDR), 0, FDR), direction = data.table::fifelse(rho > 0, "+", "-"), p.value = NULL, FDR = NULL)]

  }

  # convert data tables back to igraphs
  convertToGraph(dt_list_with_dir, V(network_list[[1]])$name, n_cores)

}
