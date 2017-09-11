distvar.meanoutput <- function(X, affine = FALSE, bias_corr = TRUE, type.X = "sample",
                    metr.X = "euclidean", bandwidth = 1) {

    #extract dimensions and sample sizes
    ss.dimX <- extract_np(X,type.X)

    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension

    if (bias_corr == TRUE &&
        type.X == "sample" &&
        metr.X == "euclidean" && n > 175 && p==1L) {
        dvar2 <- distvar_fast.meanoutput(X)
        dvar <- sqrt(abs(dvar2$dvar2)) * sign(dvar2$dvar2)
        return(list("dvar2"=dvar,"mean"=dvar2$mean))
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
        distX <- distmat(X,metr.X,n)
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
    return(list("dvar"=dvar,"mean"=mX))
}


distvar_fast.meanoutput <- function(X) {
    n <- length(X)
    temp <- IX <- 1:n

    vX <- Sort(X)
    IX0 <- sort_index(X) + 1
    IX[IX0] <- temp


    sX <- cumsum(vX)
    alphaX <- IX - 1
    betaX <- sX[IX] - vX[IX]
    Xdot <- sum(X)

    aidot <- Xdot + (2 * alphaX - n) * X - 2 * betaX
    Saa <- sum(aidot^2)

    adotdot <- 2 * sum(alphaX * X) - 2 * sum(betaX)

    gamma_1  <- PartialSum2D(X, X, rep(1, n))
    gamma_X  <- PartialSum2D(X, X, X)
    gamma_XX <- PartialSum2D(X, X, X * X)

    aijaij <- sum(X^2 * gamma_1 + gamma_XX - X * gamma_X - X * gamma_X)
    dvar <- aijaij / n / (n - 3) - 2 * Saa / n / (n - 2) / (n - 3) + adotdot * adotdot / n / (n - 1) / (n - 2) / (n - 3)
    return (list("dvar2"=dvar,"mean"=adotdot / n^2))
}


