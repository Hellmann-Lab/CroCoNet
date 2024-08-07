#' Load the GRNBoost2 output files as a list of igraph objects
#'
#' Loads the TSV files that the network inference algorithm GRNBoost2 outputs and summarizes them as a list of igraph objects.
#' If GRNBoost2 was run several times on the same data to circumvent the stochastic nature of the algorithm, the edge weights are averaged across the different runs, and a single combined igraph object is returned per clone. The TSV files corresponding to the same clone but different runs are expected to have the same base name: nameOfClone_indexOfRun.tsv.
#' In addition, edges inferred between the same gene pair but in opposite directions are also averaged. Rarely occurring edges can be removed altogether by specifying 'min_occurrence'.
#'
#' @param path Path where the GRNBoost2 output files can be found.
#' @param clone_names The names of the clones/cell lines that were input for network reconstruction. Expected to be the base names of the GRNBoost2 TSV files.
#' @param rep The number of GRNBoost2 output files per dataset, i.e. the number of independent runs. The GRNBoost2 TSV files are expected to be indexed from 1 to 'rep' for each dataset.
#' @param min_occurrence The minimum number of occurrences an edges has to have across runs and directions to be kept (the highest possible number of occurrences is 2×(# of runs)).
#' @param n_cores Number of cores.
#'
#' @return A named list that contains the networks of each clone in an igraph format, with edge weights averaged across different runs and directions, and rarely occurring edges removed. The list names are taken from the argument \code{clone_names}.
#' @export
#'
#' @examples network_list_raw <- loadGRNBoost2output(system.file("extdata", package = "CroCoNet"))
loadGRNBoost2output <- function(path, clone_names = NULL, rep = NULL, min_occurrence = 2L, n_cores = 1L) {

  if (is.null(clone_names))
    clone_names <- unique(sapply(list.files(path), function(file) {strsplit(file, "\\_")[[1]][1]}))

  if (is.null(rep))
    rep <- max(as.integer(sapply(list.files(path), function(file) {strsplit(file, "_|.tsv")[[1]][2]})))

  if (!inherits(path, "character") || length(path) != 1)
    stop("The argument \"path\" should be a string specifying the path to the directory where the GRNBoost2 output files can be found.")

  if (!file.exists(path))
    stop("The specified path does not exist.")

  if (length(list.files(path)) == 0)
    stop("The specified directory is empty.")

  if (!inherits(clone_names, "character"))
    stop("The argument \"clone_names\" should be a character vector.")

  if (length(rep) != 1 || ((!inherits(rep, "integer") & !(inherits(rep, "numeric") & rep == round(rep))) || rep < 1))
    stop("The argument \"rep\" should be a positive integer.")

  if (length(min_occurrence) != 1 || ((!inherits(min_occurrence, "integer") & !(inherits(min_occurrence, "numeric") & min_occurrence == round(min_occurrence))) || min_occurrence < 1))
    stop("The argument \"min_occurrence\" should be a positive integer.")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") & !(inherits(n_cores, "numeric") & n_cores == round(n_cores))) | n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  . = clone_name = i = TF = target = importance = from = to = weight = n_supporting_edges = NULL

  doParallel::registerDoParallel(n_cores)

  dt_list_raw <- foreach::foreach(clone_name = clone_names) %:%
                 foreach::foreach(i = 1:rep) %dopar% {

                   if (file.exists(paste0(path, "/", clone_name, "_", i, ".tsv")))
                     data.table::fread(paste0(path, "/", clone_name, "_", i, ".tsv"))

            }

  doParallel::stopImplicitCluster()

  names(dt_list_raw) <- clone_names

  doParallel::registerDoParallel(n_cores)

  dt_list <- foreach::foreach(clone_name = clone_names) %dopar% {

                    data.table::rbindlist(dt_list_raw[[clone_name]])[
                      , c("from", "to", "weight") := .(TF, target, importance)][
                        , c("from", "to") := .(pmin(from, to), pmax(from, to))][
                          , .(weight = sum(weight) / (2*rep), n_supporting_edges = .N), by = .(from, to)][
                            n_supporting_edges >= min_occurrence]

                     }

  doParallel::stopImplicitCluster()

  names(dt_list) <- clone_names

  dt_combined <- data.table::rbindlist(dt_list)
  network_genes <- sort(unique(c(dt_combined$from, dt_combined$to)))

  convertToGraph(dt_list, network_genes, n_cores)

}
