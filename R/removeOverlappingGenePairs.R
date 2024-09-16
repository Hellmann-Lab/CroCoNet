#' Remove gene pairs with overlapping annotations
#'
#' Removes the gene pairs (edges) that have overlapping annotations in the genomes of at least 1 species from all networks.
#'
#' Mapping and counting is problematic for overlapping genomic features: it is difficult to tell apart which read belongs to which gene. Often parts of the reads from one gene are assigned to the other gene, leading to correlated expression profiles simply due to genomic position. This only has a marginal effect on the results of a DE analysis, but can cause false positive edges with very high edge weights in case of a network analysis.
#'
#' This functions circumvents such potential artefacts by removing all edges between overlapping genes. As the first step, the positions of the network genes are determined in the genome of each species. Depending on the value of \code{gene_col}, the names of the network nodes are matched to the entries in either the "gene_name" or the "gene_id" column of the provided GTF files. As the second step, gene pairs with overlapping annotations are identified in each genome. Finally, these gene pairs are removed from all networks (not just from the networks of the species where they were found to overlap!).
#'
#' For all remaining edges of the networks, the function adds the genomic distance between the 2 genes that form the edge as the edge attribute "genomic_dist". If the genes are annotated on different chromosomes, the distance is set to \code{Inf}. This information can be used for further sanity checking, e.g. to check the relationship between edge weight and genomic proximity.
#'
#' @param network_list A named list of \code{\link{igraph}} objects containing the networks of all clones.
#' @param gtf_list A named list of GRanges objects containing the genome annotations of all species.
#' @param clone2species A data frame specifying which species each clone belongs to, required columns:
#' \describe{
#' \item{clone}{Character, name of the clone.}
#' \item{species}{Character, name of the species.}
#' }
#' @param gene_col Character specifying the type of identifier the network nodes have, one of "gene_name", "gene_id". The function looks for the names of the network nodes in the corresponding column of the GTF files.
#' @param n_cores Integer, the number of cores (default: 1).
#'
#' @return A named list of \code{\link{igraph}} objects containing the networks of all clones, after the removal of gene pairs with overlapping annotations. A new edge attribute is added to all \code{\link{igraph}} objects:
#' \describe{
#' \item{genomic_dist}{Numeric, the genomic distance of the 2 genes that form the edge (\code{Inf} if the 2 genes are annotated on different chromosomes/contigs).}
#' }
#' @export
#'
#' @examples
#' network_list_filt <- removeOverlappingGenePairs(network_list_raw,
#'                                                 gtf_list,
#'                                                 clone2species,
#'                                                 "gene_name")
removeOverlappingGenePairs <- function(network_list, gtf_list, clone2species, gene_col = c("gene_name", "gene_id"), n_cores = 1L) {

  # check input data
  if (!inherits(network_list, "list"))
    stop("The argument \"network_list\" should be a named list.")

  if (is.null(names(network_list)))
    stop("The argument \"network_list\" should be a named list.")

  if (any(!sapply(network_list, function(net) {inherits(net, "igraph")})))
    stop("All elements of \"network_list\" should be of class \"igraph\".")

  if (!inherits(gtf_list, "list"))
    stop("The argument \"gtf_list\" should be a named list.")

  if (is.null(names(gtf_list)))
    stop("The argument \"gtf_list\" should be a named list.")

  if (any(!sapply(gtf_list, function(gtf) {inherits(gtf, "GRanges")})))
    stop("All elements of \"gtf_list\" should be of class \"GRanges\".")

  if (!is.data.frame(clone2species))
    stop("The argument \"clone2species\" should be a data frame.")

  if (any(!(c("clone", "species") %in% colnames(clone2species))))
    stop("The argument \"clone2species\" should contain the columns \"clone\" and \"species\".")

  if (any(!(names(network_list) %in% unique(clone2species$clone))) || any(!(unique(clone2species$clone) %in% names(network_list))))
    stop("The names of \"network_list\" should match the clone names in the column \"clone\" of \"clone2species\".")

  if (any(!(names(gtf_list) %in% unique(clone2species$species))) || any(!(unique(clone2species$species) %in% names(gtf_list))))
    stop("The names of \"gtf_list\" should match the species names in the column \"species\" of \"clone2species\".")

  if (is.null(gene_col) || !gene_col %in% c("gene_name", "gene_id"))
    stop("The argument \"gene_col\" should be one of \"gene_name\", \"gene_id\".")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") && !(inherits(n_cores, "numeric") && n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  # avoid NSE notes in R CMD check
  species_name = seqnames = start = end = gene = type = gene.x = gene.y = from = to = seqnames.from = seqnames.to = start.from = start.to = end.from = end.to = weight = n_supporting_edges = clone = genomic_dist = . = NULL

  # convert igraphs to data tables
  dt_list <- convertToDT(network_list)
  column_names <- colnames(dt_list[[1]])

  # get all network genes (assumes that all networks have the same set of vertices, if the networks were originally loaded by loadNetworks this will be fulfilled)
  network_genes <- V(network_list[[1]])$name

  # get all gene pairs that overlap in any of the genomes
  doParallel::registerDoParallel(n_cores)

  dt_with_dist <- foreach::foreach(species_name = names(gtf_list)) %dopar% {

    # annotations
    gtf_filt <- gtf_list[[species_name]] %>%
      plyranges::select(gene = .data[[gene_col]], .data[["type"]]) %>%
      plyranges::filter(gene %in% network_genes & type == "gene")

    # convert to data table
    gtf_filt_dt <- data.table::as.data.table(gtf_filt)[, .(seqnames, start, end, gene)]

    # find overlapping gene pairs
    ovlp <- data.table::as.data.table(
      plyranges::join_overlap_intersect(gtf_filt,
                                        gtf_filt)
      )[
        !(gene.x == gene.y)][
          , .(gene.x, gene.y)]

    # combine the interactions detected in the different clones into 1 data table
    dist <- data.table::rbindlist(dt_list[clone2species[clone2species$species == species_name, "clone"]], idcol = "clone")

    # add genomic position of "from"
    dist <- data.table::setnames(gtf_filt_dt, c("seqnames.from", "start.from", "end.from", "from"))[dist, on = "from"]

    # add genomic position of "to"
    dist <- data.table::setnames(gtf_filt_dt, c("seqnames.to", "start.to", "end.to", "to"))[dist, on = "to"]

    # calculate genomic distance of from and to (if they overlap, distance = 0, if they are different chromosomes, distance = Inf)
    dist <- dist[, genomic_dist := as.double(pmin(abs(start.from - end.to), abs(start.to - end.from)))][
      ovlp, genomic_dist := 0, on = .(from = gene.x, to = gene.y)][
        seqnames.from != seqnames.to, genomic_dist := Inf][
          , c("clone", column_names, "genomic_dist"), with = FALSE]

    # if a gene name had 2 gene IDs, summarise the genomic positions into a single line
    if (gene_col == "gene_name") {

      dist <- dist[, .(genomic_dist = min(genomic_dist)), by = c("clone", column_names)]

    }

    dist

                                        } %>% data.table::rbindlist()

  doParallel::stopImplicitCluster()

  # remove overlapping gene pairs
  dt_with_dist_filt <- dt_with_dist %>%
    dplyr::group_by(.data[["from"]], .data[["to"]]) %>%
    dplyr::filter(sum(.data[["genomic_dist"]] == 0) == 0) %>%
    dplyr::ungroup() %>%
    data.table::as.data.table()

  # print how many gene pairs will be removed
  message("Number of gene pairs that overlap in any of the genomes:")
  message(nrow(dplyr::distinct(dt_with_dist, .data[["from"]], .data[["to"]])) - nrow(dplyr::distinct(dt_with_dist_filt, .data[["from"]], .data[["to"]])))

  # convert to list format again
  split.data.table <- utils::getFromNamespace("split.data.table", "data.table")
  dt_list_with_dist_filt <- split.data.table(dt_with_dist_filt, by = "clone", keep.by = FALSE)

  # reorder
  dt_list_with_dist_filt <- dt_list_with_dist_filt[names(network_list)]

  # convert data tables back to igraphs
  convertToGraph(dt_list_with_dist_filt, network_genes, n_cores)

}
