#' Calculate the distance covariance
#'
#' @param X contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".
#' @param Y see X.
#' @param affine logical; indicates if the affinely transformed distance covariance should be calculated or not.
#' @param bias_corr logical; indicates if the bias corrected version of the sample distance covariance should be calculated.
#' @param type.X either "sample" or "distance"; specifies the type of input for X.
#' @param type.Y see type.X.
#' @param metr.X specifies the metric which should be used for X to analyse the distance covariance. TO DO: Provide details for this.
#' @param metr.Y see metr.X.
#' @param bandwidth currently not implemented.
#' @param use : "all" uses all observations, "complete.obs" excludes NA's
#' @return numeric giving the distance covariance between samples X and Y.
#' @export
distcov <- function(X, Y, affine = FALSE, bias_corr = TRUE, type.X = "sample",
                    type.Y = "sample", metr.X = "euclidean", metr.Y = "euclidean",
                    bandwidth = 1, use="all") {





    #extract dimensions and sample sizes
    ss.dimX <- extract_np(X,type.X)
    ss.dimY <- extract_np(Y,type.Y)

    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension

    m <- ss.dimY$Sample.Size
    q <- ss.dimY$Dimension

    if (n != m) {
      stop("Samples X and Y must have the same sizes!")
    }

    if  (use == "complete.obs") {
        ccX <- ccY <- cc <- 1:n
        if (type.X == "sample") {
            ccX <- which(complete.cases(X))
        }
        if (type.Y == "sample") {
            ccY <- which(complete.cases(Y))
        }
        cc <- intersect(ccX, ccY)
        if (type.X == "sample" && p == 1) {
            X <- X[cc]
        } else if (type.X == "sample" && p > 1) {
            X <- X[cc, ]
        }
        if (type.Y == "sample" && p == 1) {
            Y <- Y[cc]
        } else if (type.X == "sample" && p > 1) {
            Y <- Y[cc, ]
        }
        n <- m <- length(cc)
    }


    if (bias_corr == TRUE && type.X == "sample" && type.Y == "sample" &&
        metr.X == "euclidean" && metr.Y == "euclidean" && n > 1000 && p == 1L && q == 1L) {
        dcov2 <- distcov_fast(X, Y)
        dcov <- sqrt(abs(dcov2)) * sign(dcov2)
        return(dcov)
    }


    ## normalize samples if calculation of affinely invariant distance covariance is desired
    if (affine == TRUE) {
      if (p > n | q > n) {
        stop("Affinely invariant distance covariance cannot be calculated for p>n")
      }
      if (type.X == "distance" | type.Y == "distance") {
        stop("Affinely invariant distance covariance cannot be calculated for type distance")
      }
      if (p > 1) {
        X <- X %*% solve(mroot(var(X)))
      } else {
        X <- X / sd(X)
      }
      if (q > 1) {
        Y <- Y %*% solve(mroot(var(Y)))
      } else {
        Y <- Y / sd(Y)
      }
    }


    ## if distance matrix is given
    if (type.X == "distance") {
      distX <- X
    } else {

        distX <- distmat(X,metr.X,alpha,n,p)
    }




    ## if distance matrix is given
    if (type.Y == "distance") {
      distY <- Y
    }else {
            distY <- distmat(Y,metr.Y,alpha,m,q)
    }

    ##calculate rowmeans
    cmX <- colmeans(distX)
    cmY <- colmeans(distY)

    ##calculate means of total matrix
    mX <- .Internal(mean(cmX))
    mY <- .Internal(mean(cmY))

    if (bias_corr == TRUE) {
      term1 <- matrix_prod_sum(distX, distY) / n / (n - 3)
      term2 <- n ^ 4 * mX * mY / n / (n - 1) / (n - 2) / (n - 3)
      term3 <- n ^ 2 * vector_prod_sum(cmX,  cmY) / n / (n - 2) / (n - 3)
      dcov2 <- term1 + term2 - 2 * term3
    }


    if (bias_corr == FALSE) {
      term1 <- matrix_prod_sum(distX, distY) / n ^ 2
      term2 <- mX * mY
      term3 <- vector_prod_sum(cmX, cmY) / n
      dcov2 <- term1 + term2 - 2 * term3
    }


    ## distance covariance (alternative construction if dcov2 is negative due to bias correction)
    dcov <- sqrt(abs(dcov2)) * sign(dcov2)
    return(dcov)
  }


#' Calculates the distance correlation
#'
#' @param X contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".
#' @param Y see X.
#' @param affine logical; indicates if the affinely transformed distance correlation should be calculated or not.
#' @param bias_corr logical; indicates if the bias corrected version of the sample distance correlation should be calculated.
#' @param type.X either "sample" or "distance"; specifies the type of input for X.
#' @param type.Y see type.X.
#' @param metr.X specifies the metric which should be used for X to analyse the distance correlation. TO DO: Provide details for this.
#' @param metr.Y see metr.X.
#' @param bandwidth currently not implemented.
#' @param use : "all" uses all observations, "complete.obs" excludes NA's
#' @return numeric giving the distance correlation between samples X and Y.
#' @export

distcorr <- function(X, Y, affine = FALSE, bias_corr = TRUE, type.X = "sample",
                     type.Y = "sample", metr.X = "euclidean", metr.Y = "euclidean",
                     bandwidth = 1, use = "all") {

    #extract dimensions and sample sizes
    ss.dimX <- extract_np(X,type.X)
    ss.dimY <- extract_np(Y,type.Y)

    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension

    m <- ss.dimY$Sample.Size
    q <- ss.dimY$Dimension

    if (use == "complete.obs") {
        ccX <- ccY <- cc <- 1:n
        if (type.X == "sample") {
            ccX <- which(complete.cases(X))
        }
        if (type.Y == "sample") {
            ccY <- which(complete.cases(Y))
        }
        cc <- intersect(ccX, ccY)
        if (type.X == "sample" && p == 1) {
            X <- X[cc]
        } else if (type.X == "sample" && p > 1) {
            X <- X[cc, ]
        }
        if (type.Y == "sample" && p == 1) {
            Y <- Y[cc]
        } else if (type.X == "sample" && p > 1) {
            Y <- Y[cc, ]
        }
        n <- m <- length(cc)
    }


    if (n != m) {
        stop("Samples X and Y must have the same sizes!")
    }

    if (bias_corr == TRUE && type.X == "sample" && type.Y == "sample" &&
        metr.X == "euclidean" && metr.Y == "euclidean" && n > 175 && p == 1L && q == 1L) {
            dcov2 <- distcov_fast(X, Y)
            dvarX2 <- distcov_fast(X, X)
            dvarY2 <- distcov_fast(Y, Y)
            dcorr2 <- dcov2 / sqrt(dvarX2 * dvarY2)
            dcorr <- sqrt(abs(dcorr2)) * sign(dcorr2)
            return(dcorr)
    }


    ## normalize samples if calculation of affinely invariant distance covariance is desired
    if (affine == TRUE) {
        if (p > n | q > n) {
            stop("Affinely invariant distance covariance cannot be calculated for p>n")
        }
        if (type.X == "distance" | type.Y == "distance") {
            stop("Affinely invariant distance covariance cannot be calculated for type distance")
        }
        if (p > 1) {
            X <- X %*% solve(mroot(var(X)))
        } else {
            X <- X / sd(X)
        }
        if (q > 1) {
            Y <- Y %*% solve(mroot(var(Y)))
        } else {
            Y <- Y / sd(Y)
        }
    }



    ## if distance matrix is given
    if (type.X == "distance") {
        distX <- X
    } else {

        distX <- distmat(X,metr.X,n,p)
    }



    ## if distance matrix is given
    if (type.Y == "distance") {
        distY <- Y
    } else {
        distY <- distmat(Y,metr.Y,m,q)
    }

    ##calculate rowmeans
    cmX <- colmeans(distX)
    cmY <- colmeans(distY)

    ##calculate means of total matrix
    mX <- .Internal(mean(cmX))
    mY <- .Internal(mean(cmY))

    if (bias_corr == TRUE) {
        term1 <- matrix.sum(hadamard_product(distX, distY)) / n / (n - 3)
        term2 <- n ^ 4 * mX * mY / n / (n - 1) / (n - 2) / (n - 3)
        term3 <- n ^ 2 * sum(vector_product(cmX,  cmY)) / n / (n - 2) / (n - 3)
        dcov2 <- term1 + term2 - 2 * term3
        dvarX <- distvar(X, affine, bias_corr, type.X, metr.X, bandwidth, use = "all")
        dvarY <- distvar(Y, affine, bias_corr, type.Y, metr.Y, bandwidth, use = "all")
        dcorr2 <- dcov2 / dvarX / dvarY
    }


    if (bias_corr == FALSE) {
        term1 <- matrix.sum(hadamard_product(distX, distY)) / n ^ 2
        term2 <- mX * mY
        term3 <- sum(vector_product(cmX, cmY)) / n
        dcov2 <- term1 + term2 - 2 * term3
        dvarX <- distvar(X, affine, bias_corr, type.X, metr.X, bandwidth, use = "all")
        dvarY <- distvar(Y, affine, bias_corr, type.Y, metr.Y, bandwidth, use = "all")
        dcorr2 <- dcov2 / dvarX / dvarY
    }


    ## distance covariance (alternative construction if dcov2 is negative due to bias correction)
    dcorr <- sqrt(abs(dcorr2)) * sign(dcorr2)
    return(dcorr)
}


#' Calculates the distance variance
#'
#' @param X contains either the first sample or its corresponding distance matrix.
#'
#' In the first case, this input can be either a vector of positive length,
#'
#' a matrix with one column or a data.frame with one column.
#'
#' In this case, type.X must be specified as "sample".
#'
#' In the second case, the input must be a distance matrix corresponding to the sample of interest.
#'
#' In this second case, type.X must be "distance".
#' @param affine logical; indicates if the affinely transformed distance variance should be calculated or not.
#' @param bias_corr logical; indicates if the bias corrected version of the sample distance variance should be calculated.
#' @param type.X either "sample" or "distance"; specifies the type of input for X.
#' @param metr.X specifies the metric which should be used for X to analyse the distance variance TO DO: Provide details for this.
#' @param bandwidth currently not implemented.
#' @param use : "all" uses all observations, "complete.obs" excludes NA's
#' @return numeric giving the distance variance of the sample X..
#' @export
distvar <- function(X, affine = FALSE, bias_corr = TRUE, type.X = "sample",
                             metr.X = "euclidean", bandwidth = 1, use = "all") {

    #extract dimensions and sample sizes
    ss.dimX <- extract_np(X,type.X)

    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension

    if (use == "complete.obs") {
        ccX <-  1:n
        if (type.X == "sample") {
            ccX <- which(complete.cases(X))
        }
        if (type.X == "sample" && p == 1) {
            X <- X[ccX]
        } else if (type.X == "sample" && p > 1) {
            X <- X[ccX, ]
        }
        n <- length(ccX)
    }


    if (bias_corr == TRUE &&
        type.X == "sample" &&
        metr.X == "euclidean" && n > 175 && p == 1L) {
        dvar2 <- distvar_fast(X)
        dvar <- sqrt(abs(dvar2)) * sign(dvar2)
        return(dvar)
    }

      ## normalize samples if calculation of affinely invariant distance covariance is desired
    if (affine == TRUE) {
        if (p > n) {
            stop("Affinely invariant distance variance cannot be calculated for p>n")
        }
        if (type.X == "distance") {
            stop("Affinely invariant distance variance cannot be calculated for type distance")
        }
        if (p > 1) {
            X <- X %*% solve(mroot(var(X)))
        } else {
            X <- X / sd(X)
        }
    }


    ## if distance matrix is given
    if (type.X == "distance") {
        distX <- X
    } else {
        distX <- distmat(X,metr.X,n,p)
    }

    ##calculate rowmeans
    cmX <- colmeans(distX)

    ##calculate means of total matrix
    mX <- .Internal(mean(cmX))

    if (bias_corr == TRUE) {
        term1 <- matrix_prod_sum(distX, distX) / n / (n - 3)
        term2 <- n ^ 4 * mX ^ 2 / n / (n - 1) / (n - 2) / (n - 3)
        term3 <- n ^ 2 * vector_prod_sum(cmX,  cmX) / n / (n - 2) / (n - 3)
        dvar2 <- term1 + term2 - 2 * term3
    }


    if (bias_corr == FALSE) {
        term1 <- matrix_prod_sum(distX, distX) / n ^ 2
        term2 <- mX * mX
        term3 <- vector_prod_sum(cmX, cmX) / n
        dvar2 <- term1 + term2 - 2 * term3
    }
    ## distance covariance (alternative construction if dcov2 is negative due to bias correction)
    dvar <- sqrt(abs(dvar2)) * sign(dvar2)
    return(dvar)
}

## function to calculate matrix roots
mroot <- function(A) {
    e <- eigen(A)
    V <- e$vectors
    V %*% diag(e$values) %*% t(V)


    B <- V %*% diag(sqrt(e$values)) %*% t(V)
    return(B)
}



distmat <- function(X,
                    metr.X = "euclidean",
                    alpha = 1,
                    n,
                    p)
{
    if (metr.X == "euclidean") {
        distX <- Dist(X)
    } else if (metr.X == "gaussian") {
        distX <- 1 - gausskernel(X, sigma = bandwidth)
    } else if (metr.X == "discrete") {
        distX <- 1 * (Dist(X) > 0)
    } else {
        if (p == 1) {
            distX <-
                outer(1:n, 1:n,  function(i, j)
                    Vectorize(match.fun(metr.X))(X[i], X[j]))
        }
        else {
            distX <- matrix(ncol = n, nrow = n)
            for (i in 1:n) {
                for (j in i:n) {
                    distX[i, j] <- distX[j, i] <- match.fun(metr)(X[i, ], X[j, ])
                }
            }
        }
    }
    return(distX)
}


centmat <-
    function(X,
             metr.X = "euclidean",
             type.X = "sample",
             bias_corr = TRUE,
             n,
             p) {
        ## if distance matrix is given
        if (type.X == "distance") {
            distX <- X
        } else {
            distX <- distmat(X, metr.X, alpha, n, p)
        }
        cmX <- colmeans(distX)
        mX <- .Internal(mean(cmX))

        if (bias_corr == TRUE) {
            cmX <- n * cmX / (n - 2)
            mX  <- n ^ 2 * mX / (n - 1) / (n - 2)
        }

        res <- normalize_matrix(distX, cmX, mX)

        if (bias_corr == TRUE) {
            diag(res) <- rep(0, n)
            res <- sqrt(n / (n - 3)) * res
        }

        return(res)
    }


extract_np <- function(X,type.X) {
    if (type.X=="sample") {
        if (is.vector(X))
        {
         n <- length(X)
         p <- 1L
         } else if (is.matrix(X)) {
        n <- nrow(X)
        p <- ncol(X)
         } else {return("X must be a vector or matrix for type sample")}
  } else if (type.X=="distance") {if (is.matrix(X)) {
        n <- nrow(X)
        p <- 1L
     } else {return("X must be a matrix for type distance")}
  }
    return(list("Sample.Size"=n,"Dimension"=p))
}

