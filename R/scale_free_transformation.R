#' Test soft-thresholding power values
#'
#' Applies a range of soft-thresholding powers to the networks provided and determines the goodness of fit for each soft power and network.
#'
#' @param network_list A named list that contains the networks of each clone and the consensus network in an \code{\link{igraph}} format.
#' @param powers A numeric or integer vector specifying the soft-thresholding powers to test.
#' @param n_cores Integer, the number of cores (default: 1).
#'
#' @return A data frame containing the \eqn{R^2} values of the scale-free model fit for each network and soft-thresholding power.
#' @export
#'
#' @examples scale_free_fit <- testPowers(network_list, c(1, seq(2, 10, by = 2)))
testPowers <- function(network_list, powers = 1:20, n_cores = 1L) {

  # check input data
  if (!inherits(network_list, "list"))
    stop("The argument \"network_list\" should be a named list.")

  if (is.null(names(network_list)))
    stop("The argument \"network_list\" should be a named list.")

  if (any(!sapply(network_list, function(net) {inherits(net, "igraph")})))
    stop("All elements of \"network_list\" should be of class \"igraph\".")

  if ((!inherits(powers, "numeric") & !inherits(powers, "integer")) || length(powers) < 1 || any(powers < 1))
    stop("The argument \"powers\" should be a numeric vector with values greater than or equal to 1.")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") && !(inherits(n_cores, "numeric") && n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  power = clone_name = NULL

  clone_names <- names(network_list)

  doParallel::registerDoParallel(n_cores)

  scale_free_fit <- foreach::foreach(power = powers,
                                     .combine = rbind) %:%
    foreach::foreach(clone_name = clone_names,
                     .combine = c) %dopar% {

                       # get connectivities
                       network <- network_list[[clone_name]]
                       igraph::E(network)$weight <- (igraph::E(network)$weight)^power
                       k <- igraph::strength(network)

                       if (sum(k != 0) == 0) {

                         return(NA)

                       } else {

                         return(scaleFreeFitRSquared(k))

                       }

                     }
  doParallel:: stopImplicitCluster()

  rownames(scale_free_fit) <- as.character(powers)
  colnames(scale_free_fit) <- clone_names

  scale_free_fit %>%
    as.data.frame()

}


#' Select power for scale-free transfrmation
#'
#' @param scale_free_fit A data frame containing the \eqn{R^2} values of the scale-free model fit for each network and soft-thresholding power.
#' @param R2_cutoff Numeric, the \eqn{R^2} cutoff each network needs to pass when selecting the power (default: 0.85).
#' @param rescue Logical indicating whether to still select a power if none of the powers achieve \eqn{R^2} values greater than the specified \code{R2_cut} for all networks. If TRUE (the default), the power with the highest minimum \eqn{R^2} across the networks is output in case the cutoff is not met, if FALSE, the function gives an error in case the cutoff is not met.
#'
#' @return Numeric, the best power to use for the scale-free transformation.
#' @export
#'
selectPower <- function(scale_free_fit, R2_cutoff = 0.85, rescue = TRUE) {

  # check input data
  if (!is.data.frame(scale_free_fit))
    stop("The argument \"scale_free_fit\" should be a data frame.")

  if (!inherits(R2_cutoff, "numeric") || length(R2_cutoff) != 1 || R2_cutoff < 0 || R2_cutoff >= 1)
    stop("The argument \"R2_cutoff\" should be a numeric value between 0 and 1.")

  if (!inherits(rescue, "logical") || length(rescue) != 1)
    stop("The argument \"rescue\" should be a logical value.")

  scale_free_fit <- tidyr::drop_na(scale_free_fit)

  goodPowerFilter <- rowSums(scale_free_fit < R2_cutoff) == 0

  if(sum(goodPowerFilter) == 0) {

    if (rescue) {

      goodPowerFilter <- which.max(apply(scale_free_fit, 1, min))
      warning("None of the provided powers achieves a high enough R2 for all networks. Choosing the power with the highest minimum R2 across networks.")

    } else {

      stop("None of the provided powers achieves a high enough R2 for all networks. Please consider including other power values, decreasing \"R2_cutoff\" or setting \"rescue\" to FALSE.")

    }

  }

  as.numeric(rownames(scale_free_fit)[goodPowerFilter][1])

}


#' Apply soft-thresholding
#'
#' @param network_list A named list of \code{\link{igraph}} objects containing the networks of all clones.
#' @param power Numeric, the power to raise the edge weights to.
#' @param n_cores Integer, the number of cores (default: 1).
#'
#' @return A named list of \code{\link{igraph}} objects containing the networks of all clones after the scale-free transformation on the edge weights
#' @export
#'
applyPower <- function(network_list, power, n_cores = 1L) {

  # check input data
  if (!inherits(network_list, "list"))
    stop("The argument \"network_list\" should be a named list.")

  if (is.null(names(network_list)))
    stop("The argument \"network_list\" should be a named list.")

  if (any(!sapply(network_list, function(net) {inherits(net, "igraph")})))
    stop("All elements of \"network_list\" should be of class \"igraph\".")

  if (!inherits(power, "numeric") || power < 1)
    stop("The argument \"power\" should be a numeric value greater than or equal to 1.")

  if (length(n_cores) != 1 || (!inherits(n_cores, "integer") && !(inherits(n_cores, "numeric") && n_cores == round(n_cores))) || n_cores < 1)
    stop("The argument \"n_cores\" should be a positive integer.")

  doParallel::registerDoParallel(n_cores)

  network_list_scale_free <- foreach::foreach(network = network_list) %dopar% {

                                     igraph::E(network)$weight <- (igraph::E(network)$weight)^power
                                     network

  }
  doParallel::stopImplicitCluster()

  names(network_list_scale_free) <- names(network_list)

  network_list_scale_free

}


#' Calculate \eqn{R^2} value of the scale-free model fit for the provided connectivities
#'
#' @param k Numeric vector, the connectivities of a network.
#'
#' @return Numeric, the \eqn{R^2} value of the scale-free model fit.
scaleFreeFitRSquared <- function(k) {

  discretized.k = cut(k, 10)
  dk = tapply(k, discretized.k, mean)
  p.dk = as.vector(tapply(k, discretized.k, length)/length(k))
  breaks1 = seq(from = min(k), to = max(k),
                length = 10 + 1)
  hist1 = graphics::hist(k, breaks = breaks1, plot = FALSE, right = TRUE)
  dk2 = hist1$mids
  dk = ifelse(is.na(dk), dk2, dk)
  dk = ifelse(dk == 0, dk2, dk)
  p.dk = ifelse(is.na(p.dk), 0, p.dk)
  log.dk = as.vector(log10(dk))
  log.p.dk= as.numeric(log10(p.dk + 1e-09))
  lm1 = try(stats::lm(log.p.dk ~ log.dk));
  if (inherits(lm1, "try-error")) browser();

  round(summary(lm1)$r.squared, 2)

}
