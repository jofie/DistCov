distvar.meanoutput <- function(X, affine = FALSE, bias_corr = TRUE, type.X = "sample",
                               metr.X = "euclidean", alpha = 1, use = "all") {

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
        dvar2mean <- distvar_fast.meanoutput(X)
        dvar2 <- dvar2mean$dvar2
        dvar <- sqrt(abs(dvar2)) * sign(dvar2)
        return(list("dvar"=dvar,"mean"=dvar2mean$mean))
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
        distX <- distmat(X, metr.X, alpha, n, p)
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
    return (list("dvar2"=dVar,"mean"=adotdot / n^2))
}


