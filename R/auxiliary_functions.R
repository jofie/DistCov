distvar.meanoutput <-
  function(X,
           affine = FALSE,
           bias.corr = TRUE,
           type.X = "sample",
           metr.X = "euclidean",
           use = "all") {

    #extract dimensions and sample sizes
    ss.dimX <- extract_np(X, type.X)

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


    if (bias.corr == TRUE &&
        type.X == "sample" &&
        metr.X == "euclidean" && n > 175 && p == 1L) {

        dvar2mean <- distvar.fast.meanoutput(X)
        dvar2 <- dvar2mean$dvar2
        dvar <- sqrt(abs(dvar2)) * sign(dvar2)
        return(list("dvar" = dvar, "mean" = dvar2mean$mean))
  }

    ## normalize samples if calculation of affinely invariant distance covariance is desired
    if (affine == TRUE) {
      if (p > n) {
        stop("Affinely invariant distance variance cannot be calculated for p>n")
      }
      if (type.X == "distance") {
        stop("Affinely invariant distance variance cannot be calculated for type 'distance'")
      }
      if (p > 1) {
        X <- X %*% Rfast::spdinv(mroot(var(X)))
      } else {
        X <- X / sd(X)
      }
    }


    ## if distance matrix is given
    if (type.X == "distance") {
      distX <- X
    } else {
      distX <- distmat(X, metr.X, n, p)
}

    ##calculate rowmeans
    cmX <- Rfast::colmeans(distX)

    ##calculate means of total matrix
    mX <- .Internal(mean(cmX))

    if (bias.corr == TRUE) {
      term1 <- matrix_prod_sum(distX, distX) / n / (n - 3)
      term2 <- n ^ 4 * mX ^ 2 / n / (n - 1) / (n - 2) / (n - 3)
      term3 <-
        n ^ 2 * vector_prod_sum(cmX,  cmX) / n / (n - 2) / (n - 3)
      dvar2 <- term1 + term2 - 2 * term3
    }


    if (bias.corr == FALSE) {
      term1 <- matrix_prod_sum(distX, distX) / n ^ 2
      term2 <- mX * mX
      term3 <- vector_prod_sum(cmX, cmX) / n
      dvar2 <- term1 + term2 - 2 * term3
    }
    ## distance covariance (alternative construction if dcov2 is negative due to bias correction)
    dvar <- sqrt(abs(dvar2)) * sign(dvar2)
    return(list("dvar" = dvar, "mean" = mX))
  }


distvar.fast.meanoutput <- function(X) {
  n <- length(X)
  temp <- IX <- 1:n

  IX0 <- Rfast::Order(X)
  vX <- X[IX0]
  IX[IX0] <- temp


  sX <- cumsum(vX)
  alphaX <- IX - 1
  betaX <- sX[IX] - vX[IX]
  Xdot <- sum(X)
  aidot <- Xdot + (2 * alphaX - n) * X - 2 * betaX
  Saa <- vector_prod_sum(aidot, aidot)
  
  adotdot <- 2 * vector_prod_sum(alphaX, X) - 2 * sum(betaX)
  
  gamma.1  <- PartialSum2D(X, X, rep(1, n))
  gamma.X  <- PartialSum2D(X, X, X)
  gamma.XX <- PartialSum2D(X, X, X * X)
  
  aijaij <- specific_vector_prod_sum(X, X, gamma.1, gamma.X, gamma.X, gamma.XX)
  dvar <-
    aijaij / n / (n - 3) - 2 * Saa / n / (n - 2) / (n - 3) + adotdot * adotdot / n / (n - 1) / (n - 2) / (n - 3)
  return(list("dvar2" = dvar, "mean" = adotdot / n ^ 2))
}

#' Calculate the square root of a semi-positive definite matrix.
#'
#' @param A square matrix.
#'
#' @return Square matrix that is a square root of A.
mroot <- function(A) {
  e <- eigen(A)
  V <- e$vectors
  V %*% diag(e$values) %*% t(V)
  B <- V %*% diag(sqrt(e$values)) %*% t(V)
  return(B)
}

#' Calculate the Gaussian kernel matrix for a given vector.
#'
#' @param X a numeric vector.
#' @param sigma the bandwidth of the kernel.
#'
#' @return The distance matrix corresponding to X.
gausskernel <- function(X, bandwidth) {
    return(exp(-1 * Rfast::Dist(X) ^ 2 / bandwidth))
}


#' Calculate the distance matrix of a given vector and a given metric.
#'
#' @param X a numeric vector or a numeric matrix.
#' @param metr.X metric that should be used to compute the distance matrix.
#' @param n number of samples, i.e. the number of rows of X..
#' @param p number of repetitions, i.e. the number of columns of X.
#' @param ... additional parameters that are used for other metrics (e.g., the bandwidth for Gaussian kernels)
#' @details For metr.X the following metrices are built in: euclidean, gaussian and discrete. However,
#' it is possible to use a function taking two numerical arguments as metr.X.
#'
#' @return The distance matrix corresponding to X.
distmat <- function(X,
                    metr.X = "euclidean",
                    n,
                    p,
                    ...)
{
    args <- list(...)
    bandwidth = args$bandwidth
   if (metr.X == "euclidean" && p == 1) {
    distX <- Rfast::Dist(X)
  } else if (metr.X == "gaussian" && p == 1) {
    distX <- 1 - gausskernel(X, bandwidth)
  } else if (metr.X == "discrete" && p == 1) {
    distX <- 1 * (Rfast::Dist(X) > 0)
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
  return(distX)
}

#' Calculate a centralized version of the distance matrix.
#'
#' @param X a numeric vector or a numeric matrix.
#' @param metr.X metric that should be used to compute the distance matrix.
#' @param bias.corr logical; indicates if the corresponding estimator of the covariance matrix is biased or unbiased.
#' @param ... additional parameters that are used for other metrics (e.g., the bandwidth for Gaussian kernels)
#' @details For metr.X the following metrices are built in: euclidean, gaussian and discrete. However,
#' it is possible to use a function taking two numerical arguments as metr.X.
#'
#' @return The centralized distance matrix corresponding to X.
centmat <- function(X,
                    metr.X = "euclidean",
                    type.X = "sample",
                    bias.corr = TRUE,
                    n,
                    p,
                    ...) {
   ## if distance matrix is given
  if (type.X == "distance") {
    distX <- X
  } else {
    distX <- distmat(X, metr.X, n, p, ...)
  }
  cmX <- Rfast::colmeans(distX)
  mX <- .Internal(mean(cmX))

  if (bias.corr == TRUE) {
    cmX <- n * cmX / (n - 2)
    mX  <- n ^ 2 * mX / (n - 1) / (n - 2)
  }

  res <- normalize_matrix(distX, cmX, mX)

  if (bias.corr == TRUE) {
    diag(res) <- rep(0, n)
    res <- sqrt(n / (n - 3)) * res
  }

  return(res)
}

#' Extract the dimensions of X.
#'
#' @param X a numeric vector or a numeric matrix.
#' @param type.X either "sample" or "distance". If type.X = "sample", X must be
#' a numeric vector or numeric matrix with the corresponding observations. If metr.X = "distance",
#' X must be a distance matrix.
#'
#' @return The centralized distance matrix corresponding to X.
extract_np <- function(X, type.X) {
  if (type.X == "sample") {
    if (is.vector(X))
    {
      n <- length(X)
      p <- 1L
    } else if (is.matrix(X)) {
      n <- nrow(X)
      p <- ncol(X)
    } else {
      stop("X must be a vector or matrix for type 'sample'!")
    }
  } else if (type.X == "distance") {
    if (is.matrix(X)) {
      n <- nrow(X)
      p <- 1L
    } else {
      stop("X must be a matrix for type 'distance'!")
    }
  } else {
    stop("type.X must be either 'sample' or 'distance'.")
  }
  return(list("Sample.Size" = n, "Dimension" = p))
}
