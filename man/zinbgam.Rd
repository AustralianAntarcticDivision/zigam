% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM_gam.R
\name{zinbgam}
\alias{zinbgam}
\title{Zero-inflated Negative Binomial GAM}
\usage{
zinbgam(
  mu.formula,
  pi.formula,
  data,
  knots = NULL,
  mu = NULL,
  pi = NULL,
  theta = 1,
  gamma.pi = 1,
  gamma.mu = 1,
  select = FALSE,
  method = "GCV.Cp",
  min.em = 5,
  max.em = 50,
  tol = 0.01
)
}
\arguments{
\item{mu.formula}{formula for the count model}

\item{pi.formula}{formula for the binary model}

\item{data}{a data frame or list containing the model response
variable and covariates required by the formula.}

\item{knots}{an optional list of knot values to be used for basis
construction.}

\item{mu}{initial mu vector}

\item{pi}{intial pi vector}

\item{theta}{initial theta value}

\item{gamma.pi}{binary model gamma}

\item{gamma.mu}{count model gamma}

\item{select}{include model selection penalty}

\item{method}{method for selecting smoothing parameters}

\item{min.em}{minimum number of EM iterations}

\item{max.em}{maximum number of EM iterations}

\item{tol}{tolerance (default=1.0E-2)}
}
\description{
Fit a zero-inflated Negative Binomial Generalized Additive Model
using the EM Algorithm
}
