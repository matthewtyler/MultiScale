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

    if (ndim > 0) {
        prior$sigma.inv.gamma <- diag(data$D)
        prior$mu.gamma <- rep(0, data$d)
        prior$sigma.inv.ab <- (1/25.0) * diag(data$D)
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

#' Loss function for Bayes MAP estimation
#'
#' Do not call directly! Used by fit_intercepts
#'
#' @inheritParams multiscale_sparse


probit_log_prob <- function(params, data, priors) {

    alpha <- params

    lp <- 0
    for (j in 1:data$J) {
        lp <- lp + dnorm(alpha[j],
                         mean = prior$mu.ab,
                         sd = sqrt(1/ prior$ sigma.inv.ab),
                         log = TRUE)

        ll.i <- ifelse(data$Y[, j] == 1, pnorm(alpha[j], mean = 0, sd = 1, log = TRUE), ifelse(data$Y[, j] == -1, log(1 - pnorm(alpha[j])), NA))

        lp <- lp + sum(ll.i, na.rm = TRUE)
    }

    return(-1 * lp)
}

#' Fit an intercept-only binary choice model
#'
#' Do not call directly! Used by multiscale
#'
#' @inheritParams multiscale_sparse
#' @param bayes Logical, set to true if the MAP estimate is desired
#'     instead of the MLE.
#'
#' @return A list only containing intercepts \code{alpha}.


fit_intercepts <- function(data, prior = NULL, bayes = FALSE, ...) {
    ## intercept from probit, use as start value if fullbayes
    probalph = rep(NA, data$J)
    for (j in 1:data$J) {
        DV <- (data$Y[, j] + 1)/2
        DV <- DV[complete.cases(DV)]  # get errors without this line for some reason
        mod = glm(DV ~ 1, family = binomial(link = "probit"))
        probalph[j] = coef(mod)[1]
    }

    if (bayes) {
        message("Finding MAP estimate for 0-d model\n\n")

        out <- optim(par = probalph, fn = probit.log.prob, data = data, priors = priors,
            method = "BFGS")

        if (out$convergence != 0) {
            warning("Warning! MAP estimation didn't converge after", out$counts,
                "iterations. Convergence code:", out$convergence)
        }
        probalpha <- out$par

    }

    return(list(alpha = probalph))
}


#' Fit a zero- or multi-dimensional spatial voting model.
#'
#' Wrapper function for \code{multiscale_sparse} or
#' \code{multiscalle_dense}.
#'
#' @param method One of \code{c("sparse", "dense", "intercepts")},
#'     with \code{"sparse"} as the default. If \code{"intercepts"},
#'     estimate an intercept-only model. Otherwise, estimate a model
#'     with positive dimension, choosing \code{"dense"} only if there
#'     is not a lot of missing data.
#' @param ... Parameters to pass to estimation methods.
#' @inheritParams multiscale_sparse
#' @return A list of estimated parameter values for \code{alpha},
#'     \code{beta}, and \code{gamma} --- the latter two not included
#'     if \code{"intercepts"} is chosen
#' @export
#' @examples \dontrun{
#' data(s109) # from pscl package
#' data(s109)
#' Y <- s109$votes  ## WARNING: data must be +/- 1, or NA
#' Y[Y %in% 1:3] <- 1
#' Y[Y %in% 4:6] <- -1
#' Y[Y %in% c(0, 7:9)] <- NA
#' data <- list(Y = Y, N = dim(Y)[1], J = dim(Y)[2], D = 2)
#' prior <- make_priors(data)
#' init <- make_starts(data)
#' lout <- multiscale(method = "sparse", prior = prior, data = data, init = init)
#' var(lout$gamma)
#' }

multiscale <- function(method = "sparse", ...) {
    if (method == "sparse") {
        multiscale_sparse(...)
    } else if (method == "dense") {
        multiscale_dense(...)
    } else if (method == "intercepts") {
        fit_intercepts(...)
    } else {
        stop("Not a valid MultiScale method! Try 'sparse', 'dense', or 'intercepts'.")
    }
}
