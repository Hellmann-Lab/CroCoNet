#' Load networks as a list of igraph objects
#'
#' Loads the TSV files containing the gene-gene edge weights and summarizes them as a list of \code{\link{igraph}} objects per clone (i.e. cell line or biological replicate within a species).
#'
#' The functions reads and formats the TSV files in the directory specified by \code{path}. The base names of the TSV files of interest can be specified using the parameter \code{clone_names} (in this case the files are expected to be named using the format: nameOfClone.tsv or nameOfClone_index.tsv). If  \code{clone_names} is set to NULL (default), all TSV files in \code{path} are read and the clone names are deduced from the file names by stripping anything starting with "_" and ".tsv".
#'
#' As the first 3 columns, each input TSV file is expected to contain the first gene that forms the edge, the second gene that forms the edge and the edge weight. There can be additional columns as well, these will be preserved as edge attributes in the \code{\link{igraph}} objects.
#'
#' If the edge weights were calculated on different subsamplings of cells per clone, or a network inference algorithm involving stochastic steps (e.g. GRNBoost2) was run several times on each clone, the edge weights are averaged across these subsamplings/runs, and a single combined \code{\link{igraph}} object is returned per clone. The number of subsamplings/runs has to be specified by the parameter \code{rep}. The TSV files corresponding to the same clone but different subsamplings/runs are expected to have the same base name with different indices: nameOfClone_index.tsv.
#'
#' If the network inference method produces an output with directed edges, i.e. geneA-geneB and geneB-geneA can both be present, the parameter \code{directed} should be set to TRUE (this is the case e.g. for GRNBoost2). In this case the edge weights inferred between the same gene pair but in opposite directions are averaged. If the edge in one of the directions is missing, it is regarded as a 0. For correlation-based methods (e.g. \code{\link{correlatePairs}}), \code{directed} has to most likely be set to FALSE. In both cases, \code{loadNetworks} returns an undirected network.
#'
#' Rarely occurring edges can be removed altogether by specifying \code{min_occurrence}. This is not relevant if the network inference was done only once per clone (\code{rep} = 1), therefore in this case the value of \code{min_occurrence} is ignored. If the network inference was done several times for each clone (\code{rep} > 1), the highest possible number of occurrences for each edge is 2×\code{rep} in case of a directed network inference method and \code{rep} in case of an undirected network inference method. If an edge occurs less often than the specified value of \code{min_occurrence}, the edge is removed. This can be helpful to 1) denoise the networks and 2) decrease the computational power needed for the next steps.
#'
#' @param path Character, the path where the files containing the gene-gene edge weights can be found.
#' @param clone_names Character vector, the names of the clones that were used for network reconstruction. These are expected to be the base names of the TSV files (format: nameOfClone.tsv or nameOfClone_index.tsv). If set to NULL (default), the names are deduced from the file names by stripping anything starting with "_" and ".tsv".
#' @param rep Integer, the number of output files per clone, e.g. the number of different subsamplings or the number of independent runs in case of a stochastic network inference algorithm (default: 1). If \code{rep} > 1, the TSV files are expected to be indexed from 1 to \code{rep} for each clone (format: nameOfClone_index.tsv).
#' @param directed Logical indicating whether the network inference method produces a directed output, i.e. whether geneA-geneB and geneB-geneA can both be present among the edges (default: TRUE). For GRNBoost2 this has to be left as TRUE, while for correlation-based methods this has to most likely be set to FALSE.
#' @param min_occurrence Integer, the minimum number of occurrences an edges has to have across runs/subsamplings to be kept (default: 2). Disregarded if \code{rep} = 1.
#' @param n_cores Integer, the number of cores (default: 1).
#'
#' @return A named list of \code{\link{igraph}} objects containing the networks per clone. The list names are taken from \code{clone_names}. Each \code{\link{igraph}} object contains the following edge attributes:
#' \describe{
#' \item{weight}{Numeric, the edge weight (taken from the 3rd column of the input TSV files).}
#' \item{n_supporting_edges}{Integer, the number of times an edge appears across the different subsamplings/runs and the 2 possible directions (only present if the network inference was run several times per clone). The highest possible number of occurrence for each edge is 2×\code{rep} in case of a directed network inference method and \code{rep} in case of an undirected network inference method.}
#' }
#' Starting from the 4th column of the input TSV files, all columns are preserved as edge attributes with unchanged names.
#' @export
#'
#' @examples
#' network_list_raw <- loadNetworks(path = system.file("extdata", package = "CroCoNet"),
#'                                  rep = 10)
#' @references
#' Moerman, T., Santos, S. A., González-Blas, C. B., Simm, J., Moreau, Y., Aerts, J., & Aerts, S. (2019). GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks. Bioinformatics , 35(12), 2159–2161.
#'
#' Lun, A. T. L., McCarthy, D. J., & Marioni, J. C. (2016). A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor. F1000Research, 5, 2122.
loadNetworks <- function(path, clone_names = NULL, rep = 1L, directed = TRUE, min_occurrence = 2L, n_cores = 1L) {

  # check input data
  if (is.null(clone_names))
    clone_names <- unique(sapply(list.files(path), function(file) {strsplit(file, "_|.tsv")[[1]][1]}))

  if (!inherits(path, "character") || length(path) != 1)
    stop("The argument \"path\" should be a string specifying the path to the directory where the output files can be found.")

  if (!file.exists(path))
    stop("The specified path does not exist.")

  if (length(list.files(path)) == 0)
    stop("The specified directory is empty.")

  if (any(!grepl(".tsv", list.files(path))))
    stop("The directory should contain TSV files only.")

  if (!inherits(clone_names, "character"))
    stop("The argument \"clone_names\" should be a character vector.")

  if (length(rep) != 1 || ((!inherits(rep, "integer") && !(inherits(rep, "numeric") && rep == round(rep))) || rep < 1))
    stop("The argument \"rep\" should be a positive integer.")

  if (!inherits(directed, "logical") || length(directed) != 1)
    stop("The argument \"directed\" should be a logical value.")

  if (length(min_occurrence) != 1 || ((!inherits(min_occurrence, "integer") && !(inherits(min_occurrence, "numeric") && min_occurrence == round(min_occurrence))) || min_occurrence < 1))
    stop("The argument \"min_occurrence\" should be a positive integer.")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") && !(inherits(n_cores, "numeric") && n_cores == round(n_cores))) | n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  # avoid NSE notes in R CMD check
  . = clone_name = i = gene1 = gene2 = rho = p.value = FDR = from = to = weight = n_supporting_edges = NULL

  # if the network inference was run only once per clone
  if (rep == 1) {

    doParallel::registerDoParallel(n_cores)

    dt_list <- foreach::foreach(clone_name = clone_names) %dopar% {

      # load TSV files as data tables and rename columns
      dt <- data.table::fread(paste0(path, "/", clone_name, ".tsv")) %>%
        data.table::setnames(., 1:3, c("from", "to", "weight"))

      # if the networks are directed, average across the 2 directions
      if (directed) {

        dt <- dt[
          , c("from", "to") := .(pmin(from, to), pmax(from, to))][
            , .(weight = sum(weight) / 2), by = .(from, to)]

      }

        dt

      }

    doParallel::stopImplicitCluster()

    names(dt_list) <- clone_names

  # if the network inference was run several times per clone
  } else {

    doParallel::registerDoParallel(n_cores)

    dt_list <- foreach::foreach(clone_name = clone_names) %dopar% {

      # load TSV files as data tables and rename columns
      dt_all_reps <- lapply(1:rep, function(i) {

        if (file.exists(paste0(path, "/", clone_name, "_", i, ".tsv")))
          data.table::fread(paste0(path, "/", clone_name, "_", i, ".tsv")) %>%
          data.table::setnames(., 1:3, c("from", "to", "weight"))

      }) %>% data.table::rbindlist()

      # if the networks are directed, average across the runs and the 2 directions
      if (directed) {

        dt_all_reps[
          , c("from", "to") := .(pmin(from, to), pmax(from, to))][
            , .(weight = sum(weight) / (2*rep), n_supporting_edges = .N), by = .(from, to)][
              n_supporting_edges >= min_occurrence][
                order(from, to)]

        # if the networks are undirected, average across the runs
      } else {

        dt_all_reps[
          , .(weight = sum(weight) / rep, n_supporting_edges = .N), by = .(from, to)][
            n_supporting_edges >= min_occurrence][
              order(from, to)]

      }

    }

    doParallel::stopImplicitCluster()

    names(dt_list) <- clone_names

  }

  # get all network genes
  dt_combined <- data.table::rbindlist(dt_list)
  network_genes <- sort(unique(c(dt_combined$from, dt_combined$to)))

  # convert data tables to igraphs
  convertToGraph(dt_list, network_genes, n_cores)

}
