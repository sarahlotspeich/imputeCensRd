% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stride_surv.R
\name{stride_surv}
\alias{stride_surv}
\title{Kernel-smoothed Nelson-Aalen estimator of conditional survival}
\usage{
stride_surv(t, W, Delta, Z)
}
\arguments{
\item{t}{Value at which the survival function is to be estimated.}

\item{W}{Numeric vector of observed covariate values (including censored opens).}

\item{Delta}{Numeric vector of censoring indicators to accompany \code{W}. Note that \code{Delta = 0} is interpreted as a censored observation.}

\item{Z}{Numeric vector of additional fully observed (continuous) covariate to accompany \code{W} and \code{Delta}.}
}
\value{
Dataframe
}
\description{
Uses the STRIDE implementation to estimate survival of \code{W} given \code{Z} nonparametrically using the kernel-smoothed Nelson-Aalen estimator.
}
