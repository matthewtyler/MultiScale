# MultiScale

Extends the EM algorithm of Imai, Lo, and Olmsted (2016, *American
Political Science Review*) for the quadratic utility, Gaussian
error-differences spatial voting model to multiple dimensions and
potentially sparse data.

Example code:
```R
devtools::install_github("matthewtyler/MultiScale")
library(MultiScale)
library(pscl)

set.seed(20347)

data(s109) # from pscl package
Y <- s109$votes  ## WARNING: data must be +/- 1, or NA
Y[Y %in% 1:3] <- 1
Y[Y %in% 4:6] <- -1
Y[Y %in% c(0, 7:9)] <- NA

# Two Dimensions
data <- list(Y = Y, N = dim(Y)[1], J = dim(Y)[2], D = 2)
prior <- make_prior(data)
init <- make_starts(data)
lout <- multiscale(method = "sparse", prior = prior, data = data, init = init)
cor(lout$gamma)

# Sparse vs. Dense Data Comparison
Y.sparse <- Y
Y.sparse[as.logical(rbinom(length(Y), 1, 0.8))] <- NA
data.sparse <- list(Y = Y.sparse, N = dim(Y)[1], J = dim(Y)[2], D = 1)
data$D  <- 1

prior <- make_prior(data)
init <- make_starts(data)

lout.dense <- multiscale(method = "sparse", prior = prior, data = data, init = init)
lout.sparse <- multiscale(method = "sparse", prior = prior, data = data.sparse, init = init)
cor(lout.dense$gamma, lout.sparse$gamma) ## 0.995
```
