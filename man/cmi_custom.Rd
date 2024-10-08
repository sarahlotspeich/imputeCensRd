% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmi_custom.R
\name{cmi_custom}
\alias{cmi_custom}
\title{Custom conditional mean imputation (CMI) for a censored predictor with user-specified survival function}
\usage{
cmi_custom(W, Delta, Z, data, useSURV, trapezoidal_rule = FALSE)
}
\arguments{
\item{W}{Column name of observed predictor values (including censored opens).}

\item{Delta}{Column name of censoring indicators. Note that \code{Delta = 0} is interpreted as a censored observation.}

\item{Z}{Column name of additional fully observed covariates.}

\item{data}{Dataframe or named matrix containing columns \code{W}, \code{Delta}, and \code{Z}.}

\item{useSURV}{Assumed survival function for \code{W} given \code{Z}. The only arguments to \code{useSURV} should be \code{W} and \code{Z}, in that order.}

\item{trapezoidal_rule}{A logical input for whether the trapezoidal rule should be used to approximate the integral in the imputed values. Default is \code{FALSE}.}
}
\value{
\item{imputed_data}{A copy of \code{data} with added column \code{imp} containing the imputed values.}
\item{code}{Indicator of algorithm status (\code{TRUE} or \code{FALSE}).}
}
\description{
Custom conditional mean imputation (CMI) for a censored predictor using (externally calculated) user-specified conditional survival.
}
