#' Performs a distance covariance test
#'
#' @param X contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".
#' @param Y see X.
#' @param test specifies the type of test that is performed, "permutation" performs a Monte Carlo Permutation test. "gamma" performs a test based on a gamma approximation of the test statistic under the null.
#' @param b specifies the number of random permutations used for the permutation test. Ignored when test="gamma"
#' @param affine logical; indicates if the affinely transformed distance covariance should be calculated or not.
#' @param bias.corr logical; indicates if the bias corrected version of the sample distance covariance should be calculated, currently ignored when test="gamma"
#' @param type.X either "sample" or "distance"; specifies the type of input for X.
#' @param type.Y see type.X.
#' @param metr.X specifies the metric which should be used for X to analyse the distance covariance. TO DO: Provide details for this.
#' @param metr.Y see metr.X.
#' @param use : "all" uses all observations, "complete.obs" excludes NA's
#' @return list with two elements, dcov gives the distance covariance between X and Y, pval gives the p-value of the corresponding test
#' @export

distcov.test <- function(X,
                         Y,
                         test = "permutation",
                         b = 499L,
                         affine = FALSE,
                         bias.corr = TRUE,
                         type.X = "sample",
                         type.Y = "sample",
                         metr.X = "euclidean",
                         metr.Y = "euclidean",
                         use = "all") {

    #extract dimensions and sample sizes
    ss.dimX <- extract_np(X, type.X)
    ss.dimY <- extract_np(Y, type.Y)

    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension

    m <- ss.dimY$Sample.Size
    q <- ss.dimY$Dimension

    if (n != m) {

        stop("Samples X and Y must have the same sizes!")
    }

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

        if (type.X == "distance") {
            X <- X[cc,cc]
        }
        if (type.Y == "distance") {
            Y <- Y[cc,cc]
        }
    }


    ## normalize samples if calculation of affinely invariant distance covariance is desired
    if (affine == TRUE) {
        if (p > n | q > n) {
            stop("Affinely invariant distance covariance cannot be calculated for p>n")
        }
        if (type.X == "distance" | type.Y == "distance") {
            stop("Affinely invariant distance covariance cannot be calculated for type 'distance'")
        }
        if (p > 1) {
            X <- X %*% Rfast::spdinv(mroot(var(X)))
        } else {
            X <- X / sd(X)
        }
        if (q > 1) {
            Y <- Y %*% Rfast::spdinv(mroot(var(Y)))
        } else {
            Y <- Y / sd(Y)
        }
    }





    if (test == "permutation") {
        if (bias.corr == TRUE && type.X == "sample" && type.Y == "sample" &&
            metr.X == "euclidean" &&
            metr.Y == "euclidean" && n > 1e4 && p == 1L && q == 1L) {
            temp <- IX <- IY  <- 1:n

            IX0 <- Rfast::Order(X) + 1
            vX <- X[IX0]
            IX[IX0] <- temp

            IY0 <- Rfast::Order(Y) + 1
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
            Sab <- sum(aidot * bidot)

            adotdot <- 2 * sum(alphaX * X) - 2 * sum(betaX)
            bdotdot <- 2 * sum(alphaY * Y) - 2 * sum(betaY)

            gamma.1  <- PartialSum2D(X, Y, rep(1, n))
            gamma.X  <- PartialSum2D(X, Y, X)
            gamma.Y  <- PartialSum2D(X, Y, Y)
            gamma.XY <- PartialSum2D(X, Y, X * Y)

            aijbij <-
                sum(X * Y * gamma.1 + gamma.XY - X * gamma.Y - Y * gamma.X)
            dcov2 <-
                aijbij / n / (n - 3) - 2 * Sab / n / (n - 2) / (n - 3) + adotdot * bdotdot / n / (n - 1) / (n - 2) / (n - 3)
            dcov <- sign(dcov2) * sqrt(abs(dcov2))

            samples <-
                lapply(1:b, function(t) {
                    .Internal(sample(n, n, FALSE, NULL))
                    }
                )

            reps <- sapply(1:b, function(t)
            {
                vY.samp <- temp

                Y.samp <- Y[samples[[t]]]
                IY.samp <- IY[samples[[t]]]

                vY.samp[IY.samp] <- vY[temp]

                sY.samp <- cumsum(vY.samp)

                alphaY <- IY.samp - 1
                betaY <- sY[IY.samp] - vY[IY.samp]

                bidot <-
                    Ydot + (2 * alphaY - n) * Y.samp - 2 * betaY
                Sab <- sum(aidot * bidot)
                gamma.1  <-
                    PartialSum2D(X, Y.samp, rep(1, n))
                gamma.X  <- PartialSum2D(X, Y.samp, X)
                gamma.Y  <- PartialSum2D(X, Y.samp, Y.samp)
                gamma.XY <-
                    PartialSum2D(X, Y.samp, X * Y.samp)
                aijbij <-
                    sum(X * Y.samp * gamma.1 + gamma.XY - X * gamma.Y - Y.samp * gamma.X)
                res2 <-
                    aijbij / n / (n - 3) - 2 * Sab / n / (n - 2) / (n - 3) + adotdot * bdotdot / n / (n - 1) / (n - 2) / (n - 3)
                res <- sign(res2) * sqrt(abs(res2))
            })


            pval <- (1 + length(which(reps > dcov))) / (1 + b)

        }  else {
            A <- centmat(X = X, 
                         metr.X = metr.X, 
                         type.X = type.X, 
                         bias.corr = bias.corr, 
                         n = n, 
                         p = p)
            B <- centmat(X = Y, 
                         metr.X = metr.Y, 
                         type.X = type.Y, 
                         bias.corr = bias.corr, 
                         n = n, 
                         p = q)
            if (bias.corr == TRUE) {
                dcov2 <- matrix_prod_sum(A , B) / n ^ 2
            } else {
                dcov2 <-
                    (matrix_prod_sum(A , B) + vector_prod_sum(diag(A), diag(B))) / n ^ 2
            }


            dcov <- sign(dcov2) * sqrt(abs(dcov2))

            #samples <- lapply(1:b, function(t) {.Internal(sample(n, n, FALSE, NULL))})

            reps <-
                sapply(1:b, function(t) {
                    sample <- .Internal(sample(n, n, FALSE, NULL))
                    if (bias.corr == TRUE) {
                        res2 <- matrix_prod_sum_sample(A , B, sample) / n ^ 2
                    } else {
                        res2 <-
                            (
                                matrix_prod_sum_sample(A , B, sample) + vector_prod_sum_sample(diag(A), diag(B), sample)
                            ) / n ^ 2
                    }
                    return(sign(res2) * sqrt(abs(res2)))
                })


            pval <- (1 + length(which(reps > dcov))) / (1 + b)
        }
    }



    if (test == "gamma") {
        distvarX <-
            distvar.meanoutput(X = X, 
                               affine = affine, 
                               bias.corr = TRUE, 
                               type.X = type.X, 
                               metr.X = metr.X)
        distvarY <-
            distvar.meanoutput(X = Y, 
                               affine = affine, 
                               bias.corr = TRUE, 
                               type.X = type.Y, 
                               metr.X = metr.Y)

        U1 <- distvarX$dvar ^ 2 * distvarY$dvar ^ 2
        U2 <- distvarX$mean * (n / (n - 1))
        U3 <- distvarY$mean * (n / (n - 1))

        alph <- 1 / 2 * (U2 ^ 2 * U3 ^ 2) / U1
        beta <- 1 / 2 * (U2 * U3) / U1
        dcov <-
            distcov(
                X = X,
                Y = Y,
                affine = affine,
                bias.corr = TRUE,
                type.X = type.X,
                type.Y = type.Y,
                metr.X = metr.X,
                metr.Y = metr.Y,
                use = "all"
            )
        stat <- n * sign(dcov) * dcov ^ 2 + U2 * U3
        pval <- pgamma(stat, alph, beta, lower.tail = FALSE)
    }




    return(list("dcov" = dcov, "pval" = pval))




}




