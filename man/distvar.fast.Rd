% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distcov_fast.R
\name{distvar.fast}
\alias{distvar.fast}
\title{Calculate an unbiased estimate of the squared distance variance}
\usage{
distvar.fast(X)
}
\arguments{
\item{X}{numeric vector containing samples}
}
\value{
double giving the square of the distance covariance
}
\description{
Calculate an unbiased estimate of the squared distance variance
}
\details{
The calculation of this estimate is based on the algorithm in Huo 2014, i.e., it gives an unbiased estimate of the square of the distance covariance. The calculation uses U-statistics and requires only O(n log n) calculations.
}
