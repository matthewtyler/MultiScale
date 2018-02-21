# MultiScale

Extends the EM algorithm of Imai, Lo, and Olmsted (2016, *Annual
Political Science Review*) for the quadratic utility, Gaussian
error-differences spatial voting model to multiple dimensions and
potentially sparse data.

Example code:
```R
devtools::install_github("matthewtyler/MultiScale")
library(MultiScale)
library(pscl)

data(s109) # from pscl package
Y <- s109$votes  ## WARNING: data must be +/- 1, or NA
Y[Y %in% 1:3] <- 1
Y[Y %in% 4:6] <- -1
Y[Y %in% c(0, 7:9)] <- NA
data <- list(Y = Y, N = dim(Y)[1], J = dim(Y)[2], D = 2)
prior <- make_prior(data)
init <- make_starts(data)
lout <- multiscale(method = "dense", prior = prior, data = data, init = init)
cor(lout$gamma)
```
