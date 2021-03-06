% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distcov_test.R
\name{distcov.test}
\alias{distcov.test}
\title{Performs a distance covariance test}
\usage{
distcov.test(X, Y, test = "permutation", b = 499L, affine = FALSE,
  bias.corr = TRUE, type.X = "sample", type.Y = "sample",
  metr.X = "euclidean", metr.Y = "euclidean", use = "all")
}
\arguments{
\item{X}{contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".}

\item{Y}{see X.}

\item{test}{specifies the type of test that is performed, "permutation" performs a Monte Carlo Permutation test. "gamma" performs a test based on a gamma approximation of the test statistic under the null.}

\item{b}{specifies the number of random permutations used for the permutation test. Ignored when test="gamma"}

\item{affine}{logical; indicates if the affinely transformed distance covariance should be calculated or not.}

\item{bias.corr}{logical; indicates if the bias corrected version of the sample distance covariance should be calculated, currently ignored when test="gamma"}

\item{type.X}{either "sample" or "distance"; specifies the type of input for X.}

\item{type.Y}{see type.X.}

\item{metr.X}{specifies the metric which should be used for X to analyse the distance covariance. TO DO: Provide details for this.}

\item{metr.Y}{see metr.X.}

\item{use}{: "all" uses all observations, "complete.obs" excludes NA's}
}
\value{
list with two elements, dcov gives the distance covariance between X and Y, pval gives the p-value of the corresponding test
}
\description{
Performs a distance covariance test
}
