#' Update Gamma Parameters (Sparse)
#'
#' The user should never call this function directly! The sparse
#' version ignores latent variables for unobserved votes
#'
#' @param Z An \eqn{N \times J} matrix of latent variables.
#' @param alpha A \eqn{J} vector of intercepts.
#' @param beta An \eqn{J \times D} matrix of question slope
#'     parameters.
#' @inheritParams multiscale_sparse
#'
#' @return An updated \eqn{N \times D} matrix of actor ideal points.
update.gamma <- function(Z, alpha, beta, prior, data) {

    prior.mean <- prior$sigma.inv.gamma %*% prior$mu.gamma

    new.gamma <- array(NA, dim = c(data$N, data$D))
    for (i in 1:data$N) {
        use.cols <- data$cols.obs[[i]]
        post.prec <- solve(prior$sigma.inv.gamma + crossprod(beta[use.cols, ]))
        post.loc <- prior.mean + t(beta[use.cols, ]) %*% (Z[i, use.cols] - alpha[use.cols])

        new.gamma[i, ] <- post.prec %*% post.loc
    }

    return(new.gamma)
}

#' Update Alpha, Beta Parameters (Sparse)
#'
#' The user should never call this function directly! The sparse
#' version ignores latent variables for unobserved votes
#'
#' @param Z An \eqn{N \times J} matrix of latent variables.
#' @param gamma An \eqn{N \times D} matrix of actor ideal points.
#' @inheritParams multiscale_sparse
#'
#' @return An updated \eqn{J \times (D+1)} matrix of question
#'     intercepts and slopes


update.ab <- function(Z, gamma, prior, data) {

    gamma.design <- cbind(rep(1, data$N), gamma)  # add intercept
    prior.mean <- prior$sigma.inv.ab %*% prior$mu.ab

    new.ab <- array(NA, dim = c(data$J, data$D + 1))
    for (j in 1:data$J) {
        use.rows <- data$rows.obs[[j]]
        post.prec <- solve(prior$sigma.inv.ab + crossprod(gamma.design[use.rows,]))
        post.loc <- prior.mean + t(gamma.design[use.rows, ]) %*% Z[use.rows, j]
        new.ab[j, ] <- post.prec %*% post.loc
    }

    return(new.ab)
}

#' Update Latent Variables (Sparse)
#'
#' The user should never call this function directly!
#'
#' @param alpha A \eqn{J} vector of intercepts.
#' @param beta An \eqn{J \times D} matrix of question slope
#'     parameters.
#' @param gamma An \eqn{N \times D} matrix of actor ideal points.
#' @inheritParams multiscale_sparse
#'
#' @return An updated \eqn{N \times J} matrix of latent variables


update.Z <- function(alpha, beta, gamma, prior, data) {

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
