#' Calculate an unbiased estimate of the squared distance covariance
#'
#' @param X numeric vector containing samples
#' @param Y numeric vector containing samples
#' @return double giving the square of the distance covariance
#' @details The calculation of this estimate is based on the algorithm in Huo 2014, i.e., it gives an unbiased estimate of the square of the distance covariance. The calculation uses U-statistics and requires only O(n log n) calculations.
distcov_fast <- function(X, Y) {
    n <- length(X)
    temp <- IX <- IY <- 1:n

    vX <- Sort(X)
    IX0 <- sort_index(X) + 1
    IX[IX0] <- temp

    vY <- Sort(Y)
    IY0 <- sort_index(Y) + 1
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
    Sab <- sum(aidot * bidot)

    adotdot <- 2 * sum(alphaX * X) - 2 * sum(betaX)
    bdotdot <- 2 * sum(alphaY * Y) - 2 * sum(betaY)

    gamma_1  <- PartialSum2D(X, Y, rep(1, n))
    gamma_X  <- PartialSum2D(X, Y, X)
    gamma_Y  <- PartialSum2D(X, Y, Y)
    gamma_XY <- PartialSum2D(X, Y, X * Y)

    aijbij <- sum(X * Y * gamma_1 + gamma_XY - X * gamma_Y - Y * gamma_X)
    dCov <- aijbij / n / (n - 3) - 2 * Sab / n / (n - 2) / (n - 3) + adotdot * bdotdot / n / (n - 1) / (n - 2) / (n - 3)
    return (dCov)
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
    vX <- Sort(X)
    IX0 <- sort_index(X) + 1
    IX[IX0] <- temp

    X <- X[IX0]
    Y <- Y[IX0]
    Z <- Z[IX0]
    vY <- Sort(Y)
    IY0 <- sort_index(Y) + 1
    IY[IY0] <- temp
    Y <- IY

    sY <- cumsum(Z[IY0]) - Z[IY0]
    sX <- cumsum(Z) - Z
    Zdot <- sum(Z)

    gamma1 <- DyadUpdate(Y, Z)

    gamma <- Zdot - Z -2 * sY[IY] - 2 * sX + 4 * gamma1
    gamma <- gamma[IX]
}
