% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse_main.R
\name{fit_intercepts}
\alias{fit_intercepts}
\title{Fit an intercept-only binary choice model}
\usage{
fit_intercepts(data)
}
\arguments{
\item{data}{A list of data values:
\itemize{
\item \code{Y} an \eqn{N \times J} matrix of \eqn{\pm 1} or \code{NA}.
\item \code{N}, \code{J} the dimensions of \code{Y}.
\item \code{D} the dimensions of political conflict being modeled.
}}
}
\value{
A list only containing intercepts \code{alpha}.
}
\description{
Do not call directly! Used by multiscale
}
