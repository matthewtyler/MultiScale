#' Update Gamma Parameters (Dense)
#'
#' The user should never call this function directly!
#'
#' @param Z An \eqn{N \times J} matrix of latent variables.
#' @param alpha A \eqn{J} vector of intercepts.
#' @param beta An \eqn{J \times D} matrix of question slope
#'     parameters.
#' @inheritParams multiscale_dense
#'
#' @return An updated \eqn{N \times D} matrix of actor ideal points.

update_gamma <- function(Z, alpha, beta, prior, data) {

    post.prec <- solve(prior$sigma.inv.gamma + t(beta) %*% beta)

    prior.mean <- prior$sigma.inv.gamma %*% prior$mu.gamma
    debias.Z <- sweep(Z, 2, alpha)
    data.mean <- debias.Z %*% beta
    post.loc <- sweep(data.mean, 2, -prior.mean)

    new.gamma <- post.loc %*% post.prec

    return(new.gamma)
}

#' Update Alpha, Beta Parameters (Dense)
#'
#' The user should never call this function directly!
#'
#' @param Z An \eqn{N \times J} matrix of latent variables.
#' @param gamma An \eqn{N \times D} matrix of actor ideal points.
#' @inheritParams multiscale_dense
#'
#' @return An updated \eqn{J \times (D+1)} matrix of question
#'     intercepts and slopes

update_ab <- function(Z, gamma, prior, data) {

    gamma.design <- cbind(rep(1, data$N), gamma)  # add intercept
    post.prec <- solve(prior$sigma.inv.ab + crossprod(gamma.design))
    prior.mean <- prior$sigma.inv.ab %*% prior$mu.ab
    data.mean <- t(Z) %*% gamma.design
    post.loc <- sweep(data.mean, 2, -prior.mean)

    new.ab <- post.loc %*% post.prec

    return(new.ab)
}

#' Update Latent Variables (Dense)
#'
#' The user should never call this function directly!
#'
#' @param alpha A \eqn{J} vector of intercepts.
#' @param beta An \eqn{J \times D} matrix of question slope
#'     parameters.
#' @param gamma An \eqn{N \times D} matrix of actor ideal points.
#' @inheritParams multiscale_dense
#'
#' @return An updated \eqn{N \times J} matrix of latent variables

update_Z <- function(alpha, beta, gamma, prior, data) {

    M <- sweep(gamma %*% t(beta), 2, -alpha)
    M[M > 15] <- 15  # keep quantiles in range
    M[M < -15] <- -15  #
    add.to.M <- array(0, dim = c(data$N, data$J))

    for (i in 1:data$N) {
        for (j in 1:data$J) {
            if (is.na(data$Y[i, j])) {
            } else if (data$Y[i, j] == 1) {
                add.to.M[i, j] <- exp(dnorm(M[i, j], log = TRUE) - pnorm(M[i, j],
                  log.p = TRUE))
            } else if (data$Y[i, j] == -1) {
                add.to.M[i, j] <- -exp(dnorm(M[i, j], log = TRUE) - pnorm(-M[i, j],
                  log.p = TRUE))
            } else {
                message(paste0("ERROR: cell ", i, ",", j, " has non-(1,-1,NA) value ",
                  data$Y[i, j]))
            }
        }
    }

    return(M + add.to.M)
}
