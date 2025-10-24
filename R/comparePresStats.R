#' Compare preservation statistics
#'
#' Compares the two topology-based preservation statistics - correlation of adjacencies (cor_adj) and correlation of intramodular connectivities (cor_kIM) - in terms of their ability to distinguish actual modules from random ones and to capture the expected decrease in module preservation with increasing phylogenetic distance.
#'
#' As part of the CroCoNet approach, pairwise module preservation scores are calculated between replicates, both within and across species (see \code{\link{calculatePresStats}}) to gain information about the cross-species differences but also about the within-species diversity of the modules. These correlation-based preservation statistics quantify how well the module topology is preserved between the networks of two replicates. While cor_adj compares fine-grained topology at the level of adjacencies per edge, cor_kIM compares higher-level topology based on intramodular connectivities, gene-level summaries of the individual adjacencies. The statistics are calculated not just for the actual, biologically meaningful modules, but also for random modules with matching sizes.
#'
#' The function plots two distributions for both cor_adj and cor_kIM: 1) the difference in preservation between each actual and corresponding random module, and 2) the inverse correlation between preservation and phylogenetic distance for each actual module. The higher these values are, the better the preservation statistic perfoms, since the actual modules are expected to be more preserved than the random modules, and all modules, but especially the actual ones, are expected to be more preserved between closely related species than between phylogenetically distant species. By comparing the distributions between cor_adj and cor_kIM, the user can select the better preservation statistic for the downstream steps of the workflow (tree reconstruction and quantification of module conservation).
#'
#' @param pres_stats Data frame of the preservation statistics for the actual (pruned) modules. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{replicate1, replicate2}{Character, the names of the replicates compared.}
#' \item{species1, species2}{Character, the names of the species \code{replicate1} and \code{replicate2} belongs to, respectively.}
#' \item{cor_adj}{Numeric, correlation of adjacencies per module and replicate pair.}
#' \item{cor_kIM}{Numeric, correlation of intramodular connectivities per module and replicate pair.}
#' }
#' @param random_pres_stats Data frame of the preservation statistics for the random modules. Required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{replicate1, replicate2}{The names of the replicates compared.}
#' \item{species1, species2}{The names of the species \code{replicate1} and \code{replicate2} belongs to, respectively.}
#' \item{cor_adj}{Numeric, correlation of adjacencies per random module and replicate pair.}
#' \item{cor_kIM}{Numeric, correlation of intramodular connectivities per random module and replicate pair.}
#' }
#' @param tree Object of class \code{\link{phylo}}, the phylogenetic tree of the species.
#' @param colors Character vector of length 2, the colors for \code{cor_adj} and \code{cor_kIM}.
#' @param font_size Numeric, font size (default: 14).
#'
#' @return A boxplot as a \code{\link{ggplot}} object comparing how well each preservation statistic can distinguish actual and random modules and capture phylogenetic information.
#' @export
#'
#' @examples comparePresStats(pres_stats, random_pres_stats, tree)
#' @family functions to plot preservation statistics
comparePresStats <- function(pres_stats, random_pres_stats, tree, colors = NULL, font_size = 14) {

  # check input data
  if (!is.data.frame(pres_stats))
    stop("The argument \"pres_stats\" should be a data frame.")

  if (any(!c("regulator", "replicate1", "replicate2", "species1", "species2", "cor_adj", "cor_kIM") %in% colnames(pres_stats)))
    stop("The argument \"pres_stats\" should contain the columns \"regulator\", \"replicate1\", \"replicate2\", \"species1\", \"species2\", \"cor_adj\" and \"cor_kIM\".")

  if (!is.data.frame(random_pres_stats))
    stop("The argument \"random_pres_stats\" should be a data frame.")

  if (any(!c("regulator", "replicate1", "replicate2", "species1", "species2", "cor_adj", "cor_kIM") %in% colnames(random_pres_stats)))
    stop("The argument \"random_pres_stats\" should contain the columns \"regulator\", \"replicate1\", \"replicate2\", \"species1\", \"species2\", \"cor_adj\" and \"cor_kIM\".")

  if (!inherits(tree, "phylo"))
    stop("The argument \"tree\" should be a phylo object.")

  if (any(!unique(c(pres_stats$species1, pres_stats$species2)) %in% tree$tip.label))
    stop("One or more species in \"pres_stats\" are not present in the phylogenetic tree. Please make sure the the species names in the columns \"species1\" and \"species2\" of \"pres_stats\" match the tip labels of \"tree\".")

  if (any(!unique(c(random_pres_stats$species1, random_pres_stats$species2)) %in% tree$tip.label))
    stop("One or more species in \"random_pres_stats\" are not present in the phylogenetic tree. Please make sure the the species names in the columns \"species1\" and \"species2\" of \"random_pres_stats\" match the tip labels of \"tree\".")

  if (!is.null(colors) && (!inherits(colors, "character") || any(!areColors(colors)) || length(colors) != 2))
    stop("The argument \"colors\" should be a character vector of valid color representations, with length of 2.")

  if (!inherits(font_size, "numeric") || length(font_size) != 1 || font_size <= 0)
    stop("The argument \"font_size\" should be a positive numeric value.")

  # avoid NSE notes in R CMD check
  . = NULL

  # if no colors are provided, take the default
  if (is.null(colors))
    colors <- c(cor_adj = "palevioletred", cor_kIM = "deeppink4")

  # differences between actual and rando modules
  random_vs_actual_pres <- dplyr::bind_rows("actual_module" = pres_stats,
                                            "random_module" = random_pres_stats,
                                            .id = "module_set") %>%
    tidyr::pivot_longer(cols = c("cor_kIM", "cor_adj"), names_to = "statistic", values_to = "value") %>%
    dplyr::select(dplyr::all_of(c("module_set", "regulator", "replicate1", "replicate2", "species1", "species2", "statistic", "value"))) %>%
    tidyr::pivot_wider(names_from = "module_set", values_from = "value") %>%
    dplyr::mutate(actual_random_diff = .data[["actual_module"]] -  .data[["random_module"]])

  # plot differences
  p1 <- random_vs_actual_pres %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[["statistic"]], y =  .data[["actual_random_diff"]], fill =  .data[["statistic"]])) +
    ggplot2::geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(color = "black", size = 14)) +
    ggplot2::ylab(expression(Delta * ~preservation~score[~actual - random])) +
    ggpubr::geom_signif(comparisons = list(c("cor_kIM", "cor_adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.08)))

  # phylogenetic distances
  species_names <- levels(pres_stats$species1)
  phylo_dists <- ape::cophenetic.phylo(tree) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("species1") %>%
    tidyr::pivot_longer(cols = 2:ncol(.), names_to = "species2", values_to = "distance") %>%
    dplyr::mutate(species1 = factor(.data[["species1"]], species_names),
                  species2 = factor(.data[["species2"]], species_names)) %>%
    dplyr::filter(as.integer(.data[["species1"]]) < as.integer( .data[["species2"]])) %>%
    dplyr::transmute(species_pair = paste0(.data[["species1"]], "_", .data[["species2"]]),  .data[["distance"]]) %>%
    tibble::deframe()

  # calculate correlations between preservation and phylogenetic distance
  pres_stats_to_plot <- pres_stats %>%
    tidyr::pivot_longer(cols = c("cor_adj", "cor_kIM"), names_to = "statistic", values_to = "value") %>%
    dplyr::mutate(species_pair = ifelse(.data[["species1"]] ==  .data[["species2"]], "within-\nspecies", paste0( .data[["species1"]], '_',  .data[["species2"]])),
                  phylo_dist = phylo_dists[.data[["species_pair"]]],
                  phylo_dist = tidyr::replace_na(.data[["phylo_dist"]], 0)) %>%
    dplyr::group_by(dplyr::across(c("regulator", "statistic"))) %>%
    dplyr::summarize(corr_pres_phyl = stats::cor(.data[["value"]], .data[["phylo_dist"]], method = "pearson"))

  # plot inverse correlations
  p2 <- (pres_stats_to_plot %>%
    ggplot2::ggplot(ggplot2::aes(x =  .data[["statistic"]], y = -.data[["corr_pres_phyl"]] , fill = .data[["statistic"]])) +
    ggplot2::geom_boxplot(color = "grey20", linewidth = 0.1, outlier.alpha = 0.5, outlier.size = 0.2) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(color = "black", size = 14)) +
    ggplot2::ylab(expression(-italic(r)[preservation~score * ", " * phylogenetic~distance])))  +
    ggpubr::geom_signif(comparisons = list(c("cor_kIM", "cor_adj")), map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05, "n.s." =1), textsize = 4, tip_length = 0.01, size = 0.3, vjust = 0.3)  +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.08)))

  p1 + p2

}
