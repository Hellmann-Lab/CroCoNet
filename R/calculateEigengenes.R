#' Calculate module eigengenes
#'
#' Calculates the module eigengene for each input module.
#'
#' A concept adapted from WGCNA, the eigengene summarizes the expression profile of an entire module, and it is calculated as the first principal component of the module expression data. In principle, it is a weighted average of the individual genes' expression profiles.
#'
#' As the first step, the logcounts are subsetted to keep only module member genes and the resulting count matrix is scaled and centered per gene. Next, singular value decomposition is performed on the scaled and centered count matrix using \code{svd}, with \code{nu} = 1 and \code{nv} = 1 (only 1 left and 1 right singular vector computed). The right singular vector is taken as the eigengene. Finally, this eigengene is aligned along the average expression of the module: if the correlation of the two vectors is negative, the eigengene is negated, if the correlation is positive, the eigengene is kept as it is.
#'
#' If a module contains both activated and repressed targets of the transcriptional regulator, calculating the eigengene (or any other summary expression profiles) across both directions of regulation does not make biological sense and leads to the dilution of signal. It makes more sense to calculate the eigengene either for the activated targets only (\code{direction_of_regulation} = "+_only") or for the activated and repressed targets separately (\code{direction_of_regulation} = "+-_separately"). In these cases, the \code{module} column of the output will specify not only the name of the transcriptional regulator but also the direction of regulation (format: nameOfRegulator(+) or nameOfRegulator(-)). If an eigengene across all targets is desired irrespective of the direction of regulation, \code{direction_of_regulation} should be set to "all_together".
#'
#' If the aim is to compare the eigengenes across species, it is recommended to calculate the eigengenes per species by setting \code{per_species} to TRUE. In this case, the scaling, centering and SVD is performed for each species separately.
#'
#' If the user plans to plot the eigengene along pseudotime and/or cell types, the corresponding columns of the \code{sce} object should be specified in the arguments \code{pseudotime_column} and \code{cell_type_column}, and then the pseudotime and cell type information of each cell will be added to the output.
#'
#' @param module_names Character vector, the names of the modules for which the eigengenes should be calculated.
#' @param pruned_modules  Data frame of the pruned modules, required columns:
#' \describe{
#' \item{regulator}{Character, transcriptional regulator.}
#' \item{target}{Character, target gene of the transcriptional regulator (member of the regulator's pruned module).}
#' \item{direction}{Character specifying the direction of interaction between the regulator and the target, either "+" or "-" (only required if \code{direction_of_regulation} is set to "+_only" or "+-_separately").}
#' }
#' @param sce \code{\link{SingleCellExperiment-class}} object containing the expression data (raw counts, logcounts and metadata) for all network genes. Required metadata columns:
#'\describe{
#' \item{species}{Character, the name of the species.}
#' \item{\{\{pseudotime_column\}\}}{Numeric, inferred pseudotime (optional).}
#' \item{\{\{cell_type_column\}\}}{Character, cell type annotation (optional).}
#' }
#' @param direction_of_regulation Character specifying how positively and negatively regulated targets of the same transcriptional regulator should be treated, one of "+_only", "+-_separately" and "all_together". If "+_only", the eigengene is calculated only for the positively regulated targets, the negatively regulated targets are simply removed. If "+-_separately", the eigengene is calculated separately for the positively and negatively regulated targets. If "all_together", the eigengene is calculated for all targets, irrespective of the direction of regulation.
#' @param per_species Logical, if FALSE (default), the SVD is performed across all cells, if TRUE, the SVD is performed per species.
#' @param pseudotime_column Character, the name of the pseudotime column in the metadata of \code{sce} (default: "pseudotime", if there is no pseudotime column, it should be set to NULL).
#' @param cell_type_column Character, the name of the cell type annotation column in the metadata of \code{sce} (default: "cell_type", if there is no cell type column, it should be set to NULL).
#' @param n_cores Integer, the number of cores (default: 1).
#'
#' @return A data frame of eigengenes with the following columns:
#'\describe{
#' \item{cell}{Character, the cell barcode.}
#' \item{species}{Character, the name of the species.}
#' \item{\{\{pseudotime_column\}\}}{Numeric, inferred pseudotime (only present if the argument \code{pseudotime_column} is not NULL).}
#' \item{\{\{cell_type_column\}\}}{Character, cell type annotation (only present if the argument \code{cell_type_column} is not NULL).}
#' \item{module}{Character, transcriptional regulator and in case the eigengene was calculated for the positively or negatively regulated targets only, the direction of interaction (format: nameOfRegulator(+) or nameOfRegulator(-)).}
#' \item{eigengene}{Numeric, the eigengene (i.e. the first principal component of the scaled and centered logcounts) of the module. In case \code{per_species} is TRUE, it is calculated for each species separately.}
#' \item{mean_expr}{Numeric, the mean of the scaled and centered logcounts across all genes in the module.}
#' \item{regulator_expr}{Numeric, the scaled and centered logcounts of the regulator.}
#' }
#' @export
#'
#' @examples
#' eigengenes <- calculateEigengenes(regulators, pruned_modules, sce)
#' eigengenes_per_dir <- calculateEigengenes(regulators, pruned_modules, sce, "+-_separately")
#' eigengenes_per_species <- calculateEigengenes(regulators, pruned_modules, sce, per_species = TRUE)
calculateEigengenes <- function(module_names, pruned_modules, sce, direction_of_regulation = "+_only", per_species = FALSE, pseudotime_column = "pseudotime", cell_type_column = "cell_type", n_cores = 1L) {

  if (!inherits(module_names, "character"))
    stop("The argument \"module_names\" should be a character vector.")

  if (any(!module_names %in% unique(pruned_modules$regulator)))
    stop("One or more modules cannot be found in \"pruned_modules\". Please make sure that the all elements of the argument \"module_names\" are present in the column \"regulator\" of \"pruned_modules\".")

  if (!is.data.frame(pruned_modules))
    stop("The argument \"pruned_modules\" should be a data frame.")

  if (any(!(c("regulator", "target") %in% colnames(pruned_modules))))
    stop("The argument \"pruned_modules\" should contain the columns \"regulator\" and \"target\".")

  if (!inherits(sce, "SingleCellExperiment"))
    stop("The argument \"sce\" should be of class \"SingleCellExperiment\".")

  if (!("logcounts" %in% names(SummarizedExperiment::assays(sce))))
    stop("The argument \"sce\" should contain the assay \"logcounts\".")

  if (is.null(sce$species))
    stop("The argument \"sce\" should contain the metadata column \"species\".")

  if (is.null(direction_of_regulation) || !direction_of_regulation %in% c("+_only", "+-_separately", "all_together"))
    stop("The argument \"direction_of_regulation\" should be one of \"+_only\", \"+-_separately\", \"all_together\".")

  if(direction_of_regulation %in% c("+_only", "+-_separately")) {

    if (!("direction" %in% colnames(pruned_modules)))
      stop("If \"direction_of_regulation\" is set to \"+_only\" or \"+-_separately\", \"pruned_modules\" should contain the column \"direction\".")

  }

  if (!inherits(per_species, "logical") || length(per_species) != 1)
    stop("The argument \"per_species\" should be a logical value.")

  if (!is.null(pseudotime_column)) {

    if (!inherits(pseudotime_column, "character") | length(pseudotime_column) > 1)
      stop("The argument \"pseudotime_column\" should be a string.")
    if (is.null(sce[[pseudotime_column]]))
      stop("The argument \"pseudotime_column\" should be the name of a column in the metadata of \"sce\".")

  }

  if (!is.null(cell_type_column)) {

    if (!inherits(cell_type_column, "character") | length(cell_type_column) > 1)
      stop("The argument \"cell_type_column\" should be a string.")
    if (is.null(sce[[cell_type_column]]))
      stop("The argument \"cell_type_column\" should be the name of a column in the metadata of \"sce\".")

  }

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") & !(inherits(n_cores, "numeric") & n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  module_name = dir = species_name = NULL

  # get metadata
  metadata <- SingleCellExperiment::colData(sce) %>%
    as.data.frame() %>%
    dplyr::select(.data[["species"]]) %>%
    tibble::rownames_to_column("cell")

  if (!is.null(pseudotime_column)) metadata[, pseudotime_column] <- sce[[pseudotime_column]]
  if (!is.null(cell_type_column)) metadata[, cell_type_column] <- sce[[cell_type_column]]

  # get logcounts
  logcnts <- SingleCellExperiment::logcounts(sce)

  # direction of interactions
  if (direction_of_regulation == "+_only") {

    pruned_modules <- pruned_modules %>%
      dplyr::filter(.data[["direction"]] == "+")

  } else if (direction_of_regulation == "all_together") {

    pruned_modules <- pruned_modules %>%
      dplyr::mutate(direction = "")

  }
  directions <- unique(pruned_modules$direction)

  # calculate eigengenes per module across all species
  if (!per_species) {

    doParallel::registerDoParallel(n_cores)

    # loop through all input modules
    eigengenes <- foreach::foreach(module_name = module_names,
                                   .combine = dplyr::bind_rows,
                                   .multicombine = TRUE) %:%

                  foreach::foreach(dir = directions,
                                   .combine = dplyr::bind_rows,
                                   .multicombine = TRUE) %dopar% {

                                     # get module member genes (if activated_only = T: regulator + activated targets, if activated_only = F: regulator + all targets)
                                       module_genes <- pruned_modules %>%
                                         dplyr::filter(.data[["regulator"]] == module_name & .data[["direction"]] == dir) %>%
                                         dplyr::pull(.data[["target"]]) %>%
                                         c(module_name)

                                       # extract module expression data (in the species of interest only)
                                       expr_module <- as.matrix(logcnts[module_genes,])

                                       # if any gene has 0s only, add some very small random noise
                                       for (i in which(rowSums(expr_module) == 0)) {

                                         expr_module[i, ] <- expr_module[i, ] + stats::rnorm(ncol(expr_module), 0, 1e-16)

                                       }

                                       # scale and center
                                       expr_module <- t(scale(t(expr_module)))

                                       # average scaled and center module expression
                                       avg_expr_module <- colMeans(expr_module)
                                       regulator_expr <- expr_module[module_name, ]

                                       # first PC of the scaled and center module expression data
                                       pc1 = svd(expr_module, nu = 1, nv = 1)$v[, 1]

                                       # align first PC to the average expression
                                       cor_avg_pc1 = stats::cor(avg_expr_module, pc1)
                                       # cor_reg_pc1 = stats::cor(regulator_expr, pc1)
                                       # if ((abs(cor_avg_pc1) > abs(cor_reg_pc1) && cor_avg_pc1 < 0) || (abs(cor_avg_pc1) < abs(cor_reg_pc1) && cor_reg_pc1 < 0 && dir == "+") || (abs(cor_avg_pc1) < abs(cor_reg_pc1) && cor_reg_pc1 > 0 && dir == "-")) pc1 = -pc1
                                       if (cor_avg_pc1 < 0) pc1 = -pc1

                                       # output metadata + eigengene + average expression + regulator expression
                                       metadata %>%
                                         dplyr::mutate(module = ifelse(direction_of_regulation == "all_together", module_name, paste0(module_name, "(", dir, ")")),
                                                       eigengene = pc1,
                                                       mean_expr = unname(avg_expr_module),
                                                       regulator_expr = unname(regulator_expr))

                                   }

    doParallel::stopImplicitCluster()

  # calculate eigengenes per module per species
  } else {

    species_names <- unique(sce$species)

    doParallel::registerDoParallel(n_cores)

    # loop through all input modules and all species
    eigengenes <- foreach::foreach(module_name = module_names,
                                   .combine = dplyr::bind_rows,
                                   .multicombine = TRUE) %:%

                  foreach::foreach(species_name = species_names,
                                   .combine = dplyr::bind_rows,
                                   .multicombine = TRUE) %:%

                  foreach::foreach(dir = directions,
                                   .combine = dplyr::bind_rows,
                                   .multicombine = TRUE) %dopar% {

                                     # get module member genes (if activated_only = T: regulator + activated targets, if activated_only = F: regulator + all targets)
                                     module_genes <- pruned_modules %>%
                                       dplyr::filter(.data[["regulator"]] == module_name & .data[["direction"]] == dir) %>%
                                       dplyr::pull(.data[["target"]]) %>%
                                       c(module_name)

                                     # extract module expression data (in the species of interest only)
                                     expr_module <- as.matrix(logcnts[module_genes, sce$species == species_name])

                                     # if any gene has 0s only, add some very small random noise
                                     for (i in which(rowSums(expr_module) == 0)) {

                                       expr_module[i, ] <- expr_module[i, ] + stats::rnorm(ncol(expr_module), 0, 1e-16)

                                     }

                                     # scale and center
                                     expr_module <- t(scale(t(expr_module)))

                                     # average scaled and center module expression
                                     avg_expr_module <- colMeans(expr_module)
                                     regulator_expr <- expr_module[module_name, ]

                                     # first PC of the scaled and center module expression data
                                     pc1 = svd(expr_module, nu = 1, nv = 1)$v[, 1]

                                     # align first PC to the average expression
                                     cor_avg_pc1 = stats::cor(avg_expr_module, pc1)
                                     # cor_reg_pc1 = stats::cor(regulator_expr, pc1)
                                     # if ((abs(cor_avg_pc1) > abs(cor_reg_pc1) && cor_avg_pc1 < 0) || (abs(cor_avg_pc1) < abs(cor_reg_pc1) && cor_reg_pc1 < 0 && dir == "+") || (abs(cor_avg_pc1) < abs(cor_reg_pc1) && cor_reg_pc1 > 0 && dir == "-")) pc1 = -pc1
                                     if (cor_avg_pc1 < 0) pc1 = -pc1

                                     # output metadata + eigengene + average expression + regulator expression
                                     metadata %>%
                                       dplyr::filter(.data[["species"]] == species_name) %>%
                                       dplyr::mutate(module = ifelse(direction_of_regulation == "all_together", module_name, paste0(module_name, "(", dir, ")")),
                                                     eigengene = pc1,
                                                     mean_expr = unname(avg_expr_module),
                                                     regulator_expr = unname(regulator_expr))

                                   }

    doParallel::stopImplicitCluster()

  }

  eigengenes %>%
    dplyr::mutate(module = factor(.data[["module"]], unique(.data[["module"]])))

}
