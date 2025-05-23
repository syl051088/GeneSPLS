% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{spls_cpp}
\alias{spls_cpp}
\title{Sparse Partial Least Squares Regression}
\usage{
spls_cpp(
  x,
  y,
  K,
  eta,
  kappa,
  select,
  fit,
  scale_x,
  scale_y,
  eps,
  maxstep,
  trace
)
}
\arguments{
\item{x}{Predictor matrix (n x p)}

\item{y}{Response matrix (n x q)}

\item{K}{Number of latent components (default: 1)}

\item{eta}{Sparsity parameter, 0 < eta < 1 (default: 0.99)}

\item{kappa}{If y is multivariate, 0 < kappa <= 0.5 (default: 0.5)}

\item{select}{Method for variable selection: "pls2" (update Y) or "simpls" (update X) (default: "pls2")}

\item{fit}{PLS estimation algorithm (only "widekernelpls" implemented) (default: "widekernelpls")}

\item{scale_x}{Whether to scale the predictor matrix (default: TRUE)}

\item{scale_y}{Whether to scale the response matrix (default: FALSE)}

\item{eps}{Convergence criterion for direction vector calculation (default: 1e-4)}

\item{maxstep}{Maximum number of iterations for direction vector calculation (default: 100)}

\item{trace}{Whether to print progress (default: FALSE)}
}
\value{
A list containing model information including:
\item{betahat}{Regression coefficients}
\item{A}{Active set of predictors}
\item{projection}{Projection matrix}
}
\description{
C++ implementation of the Sparse Partial Least Squares algorithm
}
