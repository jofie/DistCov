#' Calculate an unbiased estimate of the squared distance covariance
#'
#' @param X numeric vector containing samples
#' @param Y numeric vector containing samples
#' @return double giving the square of the distance covariance
#' @details The calculation of this estimate is based on the algorithm in Huo 2014, i.e., it gives an unbiased estimate of the square of the distance covariance. The calculation uses U-statistics and requires only O(n log n) calculations.
distcov_fast <- function(X, Y) {
    n <- length(X)
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

    gamma_1  <- PartialSum2D(X, Y, rep(1, n))
    gamma_X  <- PartialSum2D(X, Y, X)
    gamma_Y  <- PartialSum2D(X, Y, Y)
    gamma_XY <- PartialSum2D(X, Y, X * Y)

    aijbij <- specific_vector_prod_sum(X, Y, gamma_1, gamma_X, gamma_Y, gamma_XY)
    dCov <- aijbij / n / (n - 3) - 2 * Sab / n / (n - 2) / (n - 3) + adotdot * bdotdot / n / (n - 1) / (n - 2) / (n - 3)
    return(dCov)
}


#' Calculate an unbiased estimate of the squared distance variance
#'
#' @param X numeric vector containing samples
#' @return double giving the square of the distance covariance
#' @details The calculation of this estimate is based on the algorithm in Huo 2014, i.e., it gives an unbiased estimate of the square of the distance covariance. The calculation uses U-statistics and requires only O(n log n) calculations.
distvar_fast <- function(X) {
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

    gamma_1  <- PartialSum2D(X, X, rep(1, n))
    gamma_X  <- PartialSum2D(X, X, X)
    gamma_XX <- PartialSum2D(X, X, X * X)

    aijaij <- specific_vector_prod_sum(X, X, gamma_1, gamma_X, gamma_X, gamma_XX)
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
    IX0 <- Order(X)
    IX[IX0] <- temp

    #X_sort <- X[IX0]
    Y_sort <- Y[IX0]
    Z_sort <- Z[IX0]
    IY0 <- Rfast::Order(Y_sort)
    IY[IY0] <- temp
    Y_sort_2 <- IY

    sY <- cumsum(Z_sort[IY0]) - Z_sort[IY0]
    sX <- cumsum(Z_sort) - Z_sort
    Zdot <- sum(Z_sort)

    gamma1 <- DyadUpdate(Y_sort_2, Z_sort)

    gamma <- Zdot - Z_sort - 2 * sY[IY] - 2 * sX + 4 * gamma1
    gamma <- gamma[IX]
    return(gamma)
}

