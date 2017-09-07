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
#' @return numeric giving the distance covariance between samples X and Y.
#' @export

distcov <-
  function(X,
           Y,
           affine = FALSE,
           bias_corr = TRUE,
           type.X = "sample",
           type.Y = "sample",
           metr.X = "euclidean",
           metr.Y = "euclidean",
           bandwidth = 1)
  {
    ## extract dimensions and sample size
    if (type.X == "sample" && type.Y == "sample") {
      n <- length(X)
      m <- length(Y)
    } else {
      n <- nrow(as.matrix(X))
      m <- nrow(as.matrix(Y))
    }
    
    
    if (n != m) {
      stop("Samples X and Y must have the same sizes!")
    }
    if (bias_corr == TRUE &&
        type.X == "sample" &&
        type.Y == "sample" &&
        metr.X == "euclidean" && metr.Y == "euclidean" && n > 175) {
      dcov2 <- distcov_fast(X, Y)
      dcov <- sqrt(abs(dcov2)) * sign(dcov2)
      return(dcov)
    }
    p <- ncol(as.matrix(X))
    q <- ncol(as.matrix(Y))
    
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
    }
    
    ## if sample is given
    if (type.X == "sample") {
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
              distX[i, j] <- distX[j, i] <- match.fun(metr.X)(X[i,], X[j,])
            }
          }
        }
      }
    }
    
    ## if distance matrix is given
    if (type.Y == "distance") {
      distY <- Y
    }
    
    ## if sample is given
    if (type.Y == "sample") {
      if (metr.Y == "euclidean") {
        distY <- Dist(Y)
      } else if (metr.Y == "gaussian") {
        distY <- 1 - gausskernel(Y, sigma = bandwidth)
      } else if (metr.X == "discrete") {
        distY <- 1 * (Dist(Y) > 0)
      } else {
        if (q == 1) {
          distY <-
            outer(1:n, 1:n,  function(i, j)
              Vectorize(match.fun(metr.Y))(Y[i], Y[j]))
        } else {
          distY <- matrix(ncol = n, nrow = n)
          for (i in 1:n) {
            for (j in i:n) {
              distY[i, j] <- distY[j, i] <- match.fun(metr.Y)(Y[i,], Y[j,])
            }
          }
        }
      }
    }
    
    if (bias_corr == FALSE) {
      
      ##calculate rowmeans
      cmX <- colmeans(distX)
      cmY <- colmeans(distY)
      
      ##calculate means of total matrix
      mX <- .Internal(mean(cmX))
      mY <- .Internal(mean(cmY))
      
      ## Centered distance matrices
      distX_copy <- data.table::copy(distX)
      distY_copy <- data.table::copy(distY)
      A <- normalize_matrix(distX, cmX, mX)
      B <- normalize_matrix(distY, cmY, mY)
      dcov2 <- matrix.sum(hadamard_product(A, B)) / n ^ 2
      
    } else if (bias_corr == TRUE) {
      ## Bias correction if desired, now for bias correction from paper by Huo
      if (n <= 3) {
        stop("Sample size must at least be 4 for the bias corrected sample distance covariance.")
      }
      csX <- colsums(distX)
      csY <- colsums(distY)
      cmX <- csX / (n - 2)
      cmY <- csY / (n - 2)
      mX <- sum(csX) / (n - 2) / (n - 1)
      mY <- sum(csY) / (n - 2) / (n - 1)
      A <- normalize_matrix(distX, cmX, mX)
      B <- normalize_matrix(distY, cmY, mY)
      diag(A) <- diag(B) <- 0
      dcov2 <- matrix.sum(hadamard_product(A, B)) / n / (n - 3)
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
#' @return numeric giving the distance correlation between samples X and Y.
#' @export

distcorr <-
  function(X,
           Y,
           affine = FALSE,
           bias_corr = TRUE,
           type.X = "sample",
           type.Y = "sample",
           metr.X = "euclidean",
           metr.Y = "euclidean",
           bandwidth = 1) {
    
    ## extract dimensions and sample size
    if (type.X == "sample" && type.Y == "sample") {
      n <- length(X)
      m <- length(Y)
    } else {
      n <- nrow(as.matrix(X))
      m <- nrow(as.matrix(Y))
    }
    
    if (n != m) {
      stop("Samples X and Y must have the same sizes!")
    }
    
    if (bias_corr == TRUE &&
        type.X == "sample" &&
        type.Y == "sample" &&
        metr.X == "euclidean" && metr.Y == "euclidean" && n > 175) {
      dcov2 <- distcov_fast(X, Y)
      dvarX2 <- distcov_fast(X, X)
      dvarY2 <- distcov_fast(Y, Y)
      dcorr2 <- dcov2 / sqrt(dvarX2 * dvarY2)
      dcorr <- sqrt(abs(dcorr2)) * sign(dcorr2)
      return(dcorr)
    }
    
    p <- ncol(as.matrix(X))
    q <- ncol(as.matrix(Y))
    
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
    if (type.X == "distance")
    {
      distX <- X
    }
    
    ## if sample is given
    if (type.X == "sample") {
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
              distX[i, j] <- distX[j, i] <- match.fun(metr.X)(X[i, ], X[j, ])
            }
          }
        }
      }
    }
    
    ## if distance matrix is given
    if (type.Y == "distance") {
      distY <- Y
    }
    
    ## if sample is given
    if (type.Y == "sample") {
      if (metr.Y == "euclidean") {
        distY <- Dist(Y)
      } else if (metr.Y == "gaussian") {
        distY <- 1 - gausskernel(Y, sigma = bandwidth)
      } else if (metr.X == "discrete") {
        distY <- 1 * (Dist(Y) > 0)
      } else {
        if (q == 1) {
          distY <-
            outer(1:n, 1:n,  function(i, j)
              Vectorize(match.fun(metr.Y))(Y[i], Y[j]))
        } else {
          distY <- matrix(ncol = n, nrow = n)
          for (i in 1:n) {
            for (j in i:n) {
              distY[i, j] <- distY[j, i] <- match.fun(metr.Y)(Y[i, ], Y[j, ])
            }
          }
        }
      }
    }
    
    if (bias_corr == FALSE) {
      
      ##calculate rowmeans
      cmX <- colmeans(distX)
      cmY <- colmeans(distY)
      
      ##calculate means of total matrix
      mX <- .Internal(mean(cmX))
      mY <- .Internal(mean(cmY))
      
      ## Centered distance matrices
      distX_copy <- data.table::copy(distX)
      distY_copy <- data.table::copy(distY)
      A <- normalize_matrix(distX, cmX, mX)
      B <- normalize_matrix(distY, cmY, mY)
      A_copy <- data.table::copy(A)
      B_copy <- data.table::copy(B)
      dcov2 <- matrix.sum(hadamard_product(A, B))
      dvarX2 <- matrix.sum(hadamard_product(A_copy, A_copy))
      dvarY2 <- matrix.sum(hadamard_product(B_copy, B_copy))
      dcorr2 <- dcov2 / sqrt(dvarX2 * dvarY2)
      
    } else if (bias_corr == TRUE) {
      
      ## Bias correction if desired, now for bias correction from paper by Huo
      csX <- colsums(distX)
      csY <- colsums(distY)
      cmX <- csX / (n - 2)
      cmY <- csY / (n - 2)
      mX <- sum(csX) / (n - 2) / (n - 1)
      mY <- sum(csY) / (n - 2) / (n - 1)
      A <- normalize_matrix(distX, cmX, mX)
      B <- normalize_matrix(distY, cmY, mY)
      A_copy <- data.table::copy(A)
      B_copy <- data.table::copy(B)
      diag(A) <- diag(B) <- 0
      dcov2 <- matrix.sum(hadamard_product(A, B))
      dvarX2 <- matrix.sum(hadamard_product(A_copy, A_copy))
      dvarY2 <- matrix.sum(hadamard_product(B_copy, B_copy))
      dcorr2 <- dcov2 / sqrt(dvarX2 * dvarY2)
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
#' @return numeric giving the distance variance of the sample X..
#' @export

distvar <-
  function(X,
           affine = FALSE,
           bias_corr = TRUE,
           type.X = "sample",
           metr.X = "euclidean",
           bandwidth = 1) {
    
    ## extract dimensions and sample size
    if (type.X == "sample") {
      n <- length(X)
    } else {
      n <- nrow(as.matrix(X))
    }
    
    if (bias_corr == TRUE &&
        type.X == "sample" &&
        metr.X == "euclidean" && n > 175) {
      dvar2 <- distcov_fast(X, X)
      dcorr <- sqrt(abs(dvar2)) * sign(dvar2)
      return(dcorr)
    }
    
    p <- ncol(as.matrix(X))
        
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
    if (type.X == "distance")
    {
      distX <- X
    }
    
    ## if sample is given
    if (type.X == "sample") {
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
              distX[i, j] <- distX[j, i] <- match.fun(metr.X)(X[i, ], X[j, ])
            }
          }
        }
      }
    }
    
    
    
    if (bias_corr == FALSE) {
      
      ##calculate rowmeans
      cmX <- colmeans(distX)
          
      ##calculate means of total matrix
      mX <- .Internal(mean(cmX))
          
      ## Centered distance matrices
      distX_copy <- data.table::copy(distX)
      A <- normalize_matrix(distX, cmX, mX)
      dvar2 <- matrix.sum(hadamard_product(A, A))
            
    } else if (bias_corr == TRUE) {
      
      ## Bias correction if desired, now for bias correction from paper by Huo
      csX <- colsums(distX)
      cmX <- csX / (n - 2)
      mX <- sum(csX) / (n - 2) / (n - 1)
      A <- normalize_matrix(distX, cmX, mX)
      diag(A) <- 0
      dvar2 <- matrix.sum(hadamard_product(A, A))
    }
    
    ## distance covariance (alternative construction if dcov2 is negative due to bias correction)
    dvar <- sqrt(abs(dvar2)) * sign(dvar2)
    return(dcorr)
  }



## function to calculate matrix roots
mroot <- function(A) {
  e <- eigen(A)
  V <- e$vectors
  V %*% diag(e$values) %*% t(V)
  
  
  B <- V %*% diag(sqrt(e$values)) %*% t(V)
  return(B)
}
