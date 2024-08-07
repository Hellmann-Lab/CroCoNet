#' Remove gene pairs with overlapping annotations
#'
#' Removes gene pairs that have overlapping annotations in any of the species' genomes from all networks.
#' @param network_list A named list that contains the networks of each clone in an igraph format.
#' @param gtf_list A named list that contains the genome annotation of each species in a GenomicRanges format.
#' @param clone2species A data frame with columns 'clone' and 'species' that specifies which species each clone belongs to. The names of clones and species should match the names of 'network_list' and gtf_list'.
#' @param gene_col One of 'gene_name' or 'gene_id', specifies what type of identifier the network nodes have.
#' @param n_cores Number of cores.
#'
#' @return A named list of networks after the removal of gene pairs with overlapping annotations.
#' @export
#'
#' @examples
#' network_list_scaled_filt <- removeOverlappingGenePairs(network_list_scaled,
#'                                                        gtf_list,
#'                                                        clone2species,
#'                                                        "gene_name")
removeOverlappingGenePairs <- function(network_list, gtf_list, clone2species, gene_col = c("gene_name", "gene_id"), n_cores = 1L) {

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

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") & !(inherits(n_cores, "numeric") & n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")


  species_name = seqnames = start = end = gene = type = gene.x = gene.y = from = to = seqnames.from = seqnames.to = start.from = start.to = end.from = end.to = weight = n_supporting_edges = clone = distance = . = NULL

  dt_list <- convertToDT(network_list)

  network_genes <- V(network_list[[1]])$name # assuming that all networks have the same set of vertices

  doParallel::registerDoParallel(n_cores)

  dt_with_dist <- foreach::foreach(species_name = names(gtf_list)) %dopar% {

    # annotations
    gtf_filt <- gtf_list[[species_name]] %>%
      plyranges::select(gene = .data[[gene_col]], .data[["type"]]) %>%
      plyranges::filter(gene %in% network_genes & type == "gene")
    gtf_filt_dt <- data.table::as.data.table(gtf_filt)[, .(seqnames, start, end, gene)]

    # find overlapping genes
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
    dist <- dist[, .(from, to, seqnames.from, seqnames.to, weight, n_supporting_edges, clone, distance = as.double(pmin(abs(start.from - end.to), abs(start.to - end.from))))][
      ovlp, distance := 0, on = .(from = gene.x, to = gene.y)][
        seqnames.from != seqnames.to, distance := Inf][
          , .(clone, from, to, weight, n_supporting_edges, distance)]

    dist

                                        } %>% data.table::rbindlist()

  doParallel::stopImplicitCluster()

  # remove gene pairs that are overlapping in any of the genomes
  dt_with_dist_filt <- dt_with_dist %>%
    dplyr::group_by(.data[["from"]], .data[["to"]]) %>%
    dplyr::filter(sum(.data[["distance"]] == 0) == 0) %>%
    dplyr::ungroup() %>%
    data.table::as.data.table()

  # print how many gene pairs/entries will be removed
  message("Number of gene pairs that overlap in any of the species:")
  message(nrow(dplyr::distinct(dt_with_dist, .data[["from"]], .data[["to"]])) - nrow(dplyr::distinct(dt_with_dist_filt, .data[["from"]], .data[["to"]])))
  message("Number of entries removed (same gene pair but different clone = different entries):")
  message(nrow(dt_with_dist) - nrow(dt_with_dist_filt))

  # convert to list format again
  split.data.table <- utils::getFromNamespace("split.data.table", "data.table")
  dt_list_with_dist_filt <- split.data.table(dt_with_dist_filt, by = "clone", keep.by = FALSE)

  # reorder
  dt_list_with_dist_filt <- dt_list_with_dist_filt[names(network_list)]

  convertToGraph(dt_list_with_dist_filt, network_genes, n_cores)

}
