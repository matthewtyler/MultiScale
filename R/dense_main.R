#' Multidimensional scaling WITHOUT large amounts of missing data.
#'
#' The user should never call this function directly! Use the
#' MultiScale function instead. The dense (as opposed to sparse)
#' algorithm includes the latent variables for missing vote
#' observations when updating the global parameters alpha, beta,
#' gamma. It is less robust in the presence of many missing votes
#'
#' @param prior A list of priors for the parameters ab, gamma:
#' \itemize{
#' \item \code{sigma.inv.ab} a \eqn{(D+1)\times (D+1)} positive definite matrix corresponding to the inverse of the prior covariance matrix for \eqn{(\alpha, \beta)}.
#' \item \code{mu.ab} a \eqn{(D+1)} vector corresponding to the prior mean of \eqn{(\alpha, \beta)}.
#' \item \code{sigma.inv.gamma} a \eqn{D \times D} definite matrix corresponding to the inverse of the prior covariance matrix for \eqn{\gamma}.
#' \item \code{mu.gamma} a \eqn{D} vector corresponding to the prior mean of \eqn{\gamma}.
#' }
#' @param data A list of data values:
#' \itemize{
#' \item \code{Y} an \eqn{N \times J} matrix of \eqn{\pm 1} or \code{NA}.
#' \item \code{N}, \code{J} the dimensions of \code{Y}.
#' \item \code{D} the dimensions of political conflict being modeled.
#' }
#' @param init A list of initialization values for \code{alpha}, \code{beta}, \code{gamma}.
#' @param max.iter The maximum number of iterations.
#' @param tol The algorithm stops after the current iteration and the previous iteration parameters have correlation 1 - \code{tol}.
#' @param verbose print useful warnings and updates on algorithm progress.
#' @return A list of estimated parameters \code{alpha}, \code{beta}, \code{gamma}

multiscale_dense <- function(prior, data, init, max.iter = 100, tol = 1e-04, verbose = TRUE) {

    if (verbose) {
        message("WARNING:")
        message("The parameters of this model are NOT locally identified. Importantly, estimated beta and gamma will not be necessarily be the same when this algorithm is re-used. For more information, see Rivers, Douglas. 2003. Identification of Multidimensional Spatial Voting Models.")

        stopifnot(menu(c("Yes", "No"), title = "This function is sometimes faster than MultiScaleSparse, but works performs abysmally in the presence of much missing data. Proceed?") ==
            1)
    }

    alpha <- init$alpha
    beta <- init$beta
    gamma <- init$gamma
    old.params <- c(alpha, beta, gamma)

    for (iter in 1:max.iter) {
        Z <- update.Z(alpha, beta, gamma, prior, data)
        gamma <- update.gamma(Z, alpha, beta, prior, data)
        ab <- update.ab(Z, gamma, prior, data)
        alpha <- c(ab[, 1])
        beta <- as.matrix(ab[, -1])

        new.params <- c(alpha, beta, gamma)
        param.cor <- cor(new.params, old.params)

        if (param.cor > 1 - tol) {
            message(paste("Converged after", iter, "iterations!"))
            break
        } else if (iter == max.iter) {
            message(paste("Maximum number of iterations reached at", param.cor, "correlation"))
        }

        old.params <- new.params
        if (verbose)
            message(paste("Iteration:", iter))
    }


    return(list(iter = iter, alpha = alpha, beta = beta, gamma = gamma))
}
