#' Make Prior Means and Precisions
#'
#' @export
#' @inheritParams multiscale_sparse
#' @return A list of estimated prior variables
#' #' \itemize{
#' \item \code{sigma.inv.ab} a \eqn{(D+1)\times (D+1)} positive definite matrix corresponding to the inverse of the prior covariance matrix for \eqn{(\alpha, \beta)}.
#' \item \code{mu.ab} a \eqn{(D+1)} vector corresponding to the prior mean of \eqn{(\alpha, \beta)}.
#' \item \code{sigma.inv.gamma} a \eqn{D \times D} definite matrix corresponding to the inverse of the prior covariance matrix for \eqn{\gamma}.
#' \item \code{mu.gamma} a \eqn{D} vector corresponding to the prior mean of \eqn{\gamma}.
#' }

make_prior <- function(data) {

    prior <- list()

    if (data$D > 0) {
        prior$sigma.inv.gamma <- diag(data$D)
        prior$mu.gamma <- rep(0, data$D)
        prior$sigma.inv.ab <- (1/25.0) * diag(data$D + 1)
        prior$mu.ab <- rep(0, data$D + 1)
    }

    return(prior)
}

#' Make Parameter Initialization Values
#'
#' @export
#' @inheritParams multiscale_sparse
#' @return A list of initialized parameter values for \code{alpha},
#'     \code{beta}, and \code{gamma}

make_starts <- function(data) {
    alpha <- rnorm(data$J)
    beta <- matrix(rnorm(data$J * data$D), nrow = data$J)
    gamma <- matrix(rnorm(data$N * data$D), nrow = data$N)
    return(list(alpha = alpha, beta = beta, gamma = gamma))
}


#' Fit an intercept-only binary choice model
#'
#' Do not call directly! Used by multiscale
#'
#' @inheritParams multiscale_sparse
#'
#' @return A list only containing intercepts \code{alpha}.


fit_intercepts <- function(data) {
    ## intercept from probit, use as start value if fullbayes
    probalph = rep(NA, data$J)
    for (j in 1:data$J) {
        DV <- (data$Y[, j] + 1)/2
        DV <- DV[complete.cases(DV)]  # get errors without this line for some reason
        mod = glm(DV ~ 1, family = binomial(link = "probit"))
        probalph[j] = coef(mod)[1]
    }

    return(list(alpha = probalph))
}


#' Fit a zero- or multi-dimensional spatial voting model.
#'
#' The sparse (as opposed to dense) algorithm ignores the latent
#' variables for missing vote observations when updating the global
#' parameters alpha, beta, gamma. It is robust in the presence of many
#' missing votes.
#'
#' @param method One of \code{c("sparse", "intercepts")},
#'     with \code{"sparse"} as the default. If \code{"intercepts"},
#'     estimate an intercept-only model. Otherwise, estimate a model
#'     with positive dimension.
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
#' @export
#' @examples \dontrun{
#' data(s109) # from pscl package
#' data(s109)
#' Y <- s109$votes  ## WARNING: data must be +/- 1, or NA
#' Y[Y %in% 1:3] <- 1
#' Y[Y %in% 4:6] <- -1
#' Y[Y %in% c(0, 7:9)] <- NA
#' data <- list(Y = Y, N = dim(Y)[1], J = dim(Y)[2], D = 2)
#' prior <- make_prior(data)
#' init <- make_starts(data)
#' lout <- multiscale(method = "sparse", prior = prior, data = data, init = init)
#' cor(lout$gamma)
#' }

multiscale <- function(method == "sparse", prior, data, init, max.iter = 250, tol = 1e-04, verbose = TRUE) {

    if (method == "intercepts") {
        fit_intercepts(data)
    } else if (method != "sparse") {
        stop("Not a valid MultiScale method! Try 'sparse' or 'intercepts'.")
    }

    if (verbose) {
        message("")
        message(">>>WARNING:<<<")
        message("")
        message("The parameters of this model are NOT locally identified. Importantly, estimated beta and gamma will not be necessarily be the same when this algorithm is re-used. For more information, see Rivers, Douglas. 2003. Identification of Multidimensional Spatial Voting Models.")
        message("")
    }

    data$rows.obs <- list()
    data$cols.obs <- list()

    for (i in 1:data$N) {
        data$cols.obs[[i]] <- which(!is.na(data$Y[i, ]))
    }
    for (j in 1:data$J) {
        data$rows.obs[[j]] <- which(!is.na(data$Y[, j]))
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
