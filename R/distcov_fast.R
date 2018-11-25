#' Calculate an unbiased estimate of the squared distance covariance
#'
#' @param X numeric vector containing samples
#' @param Y numeric vector containing samples
#' @return double giving the square of the distance covariance
#' @details The calculation of this estimate is based on the algorithm in Huo 2014, i.e., it gives an unbiased estimate of the square of the distance covariance. The calculation uses U-statistics and requires only O(n log n) calculations.
#' @export
distcov.fast <- function(X, Y) {
  if (!is.numeric(X) || !is.numeric(Y) || !is.vector(X) || !is.vector(Y)) {
    stop("X and Y must be numeric vectors for fast distance covariance algorithms!")
  }
    n <- length(X)
    if (n == 0 || n != length(Y)) {
      stop("X and Y are not allowed to be empty and must have the same number of entries!")
    }
    temp <- IX <- IY <- 1:n

    IX0 <- Rfast::Order(X)
    vX <- X[IX0]
    IX[IX0] <- temp

    IY0 <- Rfast::Order(Y)
    vY <- Y[IY0]
    IY[IY0] <- temp

    sX <- cumsum(vX)
    sY <- cumsum(vY)
    
    alphaX <- IX - 1
    alphaY <- IY - 1
    betaX <- sX[IX] - vX[IX]
    betaY <- sY[IY] - vY[IY]
    
    Xdot <- sum(X)
    Ydot <- sum(Y)

    aidot <- Xdot + (2 * alphaX - n) * X - 2 * betaX
    bidot <- Ydot + (2 * alphaY - n) * Y - 2 * betaY
    Sab <- vector_prod_sum(aidot, bidot) # sum(aidot * bidot)

    adotdot <- 2 * vector_prod_sum(alphaX, X) - 2 * sum(betaX)
    bdotdot <- 2 * vector_prod_sum(alphaY, Y) - 2 * sum(betaY)

    gamma.1  <- PartialSum2D(X, Y, rep(1, n))
    gamma.X  <- PartialSum2D(X, Y, X)
    gamma.Y  <- PartialSum2D(X, Y, Y)
    gamma.XY <- PartialSum2D(X, Y, X * Y)

    aijbij <- specific_vector_prod_sum(X, Y, gamma.1, gamma.X, gamma.Y, gamma.XY)
    dCov <- aijbij / n / (n - 3) - 2 * Sab / n / (n - 2) / (n - 3) + adotdot * bdotdot / n / (n - 1) / (n - 2) / (n - 3)
    return(dCov)
}


#' Calculate an unbiased estimate of the squared distance variance
#'
#' @param X numeric vector containing samples
#' @return double giving the square of the distance covariance
#' @details The calculation of this estimate is based on the algorithm in Huo 2014, i.e., it gives an unbiased estimate of the square of the distance covariance. The calculation uses U-statistics and requires only O(n log n) calculations.
distvar.fast <- function(X) {
  if (!is.numeric(X) || !is.vector(X)) {
    stop("X must be a numeric vector for fast distance covariance algorithms!")
  }
    n <- length(X)
    
    if (n == 0) {
      stop("X is empty!")
    }
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
    dVar <- aijaij / n / (n - 3) - 2 * Saa / n / (n - 2) / (n - 3) + adotdot * adotdot / n / (n - 1) / (n - 2) / (n - 3)
    return(dVar)
}



#' Calculate partial sums
#'
#' @param X numeric vector containing samples
#' @param Y numeric vector containing samples
#' @param Z numeric vector
#' @return numeric vector
#' @details This is a help function used in calculating the distance covariance using the algorithm by Huo 2014.
PartialSum2D <- function(X, Y, Z) {
    n <- length(X)
    temp <- IX <- IY <- 1:n
    IX0 <- Rfast::Order(X)
    IX[IX0] <- temp

    Y.sort <- Y[IX0]
    Z.sort <- Z[IX0]
    IY0 <- Rfast::Order(Y.sort)
    IY[IY0] <- temp
    Y.sort.2 <- IY

    sY <- cumsum(Z.sort[IY0]) - Z.sort[IY0]
    sX <- cumsum(Z.sort) - Z.sort
    Zdot <- sum(Z.sort)

    gamma1 <- DyadUpdate(Y.sort.2, Z.sort)

    gamma <- Zdot - Z.sort - 2 * sY[IY] - 2 * sX + 4 * gamma1
    gamma <- gamma[IX]
    return(gamma)
}

