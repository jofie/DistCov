#' Performs a distance covariance test
#'
#' @param X contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".
#' @param Y see X.
#' @param test specifies the type of test that is performed, "permutation" performs a Monte Carlo Permutation test. "gamma" performs a test based on a gamma approximation of the test statistic under the null.
#' @param b specifies the number of random permutations used for the permutation test. Ignored when test="gamma"
#' @param affine logical; indicates if the affinely transformed distance covariance should be calculated or not.
#' @param bias_corr logical; indicates if the bias corrected version of the sample distance covariance should be calculated, currently ignored when test="gamma"
#' @param type.X either "sample" or "distance"; specifies the type of input for X.
#' @param type.Y see type.X.
#' @param metr.X specifies the metric which should be used for X to analyse the distance covariance. TO DO: Provide details for this.
#' @param metr.Y see metr.X.
#' @param bandwidth currently not implemented.
#' @param use : "all" uses all observations, "complete.obs" excludes NA's
#' @return list with two elements, dcov gives the distance covariance between X and Y, pval gives the p-value of the corresponding test
#' @export

distcov.test <- function(X,
                         Y,
                         test = "permutation",
                         b = 499L,
                         affine = FALSE,
                         bias_corr = TRUE,
                         type.X = "sample",
                         type.Y = "sample",
                         metr.X = "euclidean",
                         metr.Y = "euclidean",
                         bandwidth = 1,
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
        if (bias_corr == TRUE && type.X == "sample" && type.Y == "sample" &&
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

            gamma_1  <- PartialSum2D(X, Y, rep(1, n))
            gamma_X  <- PartialSum2D(X, Y, X)
            gamma_Y  <- PartialSum2D(X, Y, Y)
            gamma_XY <- PartialSum2D(X, Y, X * Y)

            aijbij <-
                sum(X * Y * gamma_1 + gamma_XY - X * gamma_Y - Y * gamma_X)
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
                vY_samp <- temp

                Y_samp <- Y[samples[[t]]]
                IY_samp <- IY[samples[[t]]]

                vY_samp[IY_samp] <- vY[temp]

                sY_samp <- cumsum(vY_samp)

                alphaY <- IY_samp - 1
                betaY <- sY[IY_samp] - vY[IY_samp]

                bidot <-
                    Ydot + (2 * alphaY - n) * Y_samp - 2 * betaY
                Sab <- sum(aidot * bidot)
                gamma_1  <-
                    PartialSum2D(X, Y_samp, rep(1, n))
                gamma_X  <- PartialSum2D(X, Y_samp, X)
                gamma_Y  <- PartialSum2D(X, Y_samp, Y_samp)
                gamma_XY <-
                    PartialSum2D(X, Y_samp, X * Y_samp)
                aijbij <-
                    sum(X * Y_samp * gamma_1 + gamma_XY - X * gamma_Y - Y_samp * gamma_X)
                res2 <-
                    aijbij / n / (n - 3) - 2 * Sab / n / (n - 2) / (n - 3) + adotdot * bdotdot / n / (n - 1) / (n - 2) / (n - 3)
                res <- sign(res2) * sqrt(abs(res2))
            })


            pval <- (1 + length(which(reps > dcov))) / (1 + b)

        }  else {
            A <- centmat(X, metr.X, type.X, bias_corr, n, p)
            B <- centmat(Y, metr.Y, type.Y, bias_corr, n, q)
            if (bias_corr == TRUE) {
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
                    if (bias_corr == TRUE) {
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
            distvar.meanoutput(X, affine, bias_corr = TRUE, type.X, metr.X, bandwidth)
        distvarY <-
            distvar.meanoutput(Y, affine, bias_corr = TRUE, type.Y, metr.Y, bandwidth)

        U1 <- distvarX$dvar ^ 2 * distvarY$dvar ^ 2
        U2 <- distvarX$mean * (n / (n - 1))
        U3 <- distvarY$mean * (n / (n - 1))

        alpha <- 1 / 2 * (U2 ^ 2 * U3 ^ 2) / U1
        beta <- 1 / 2 * (U2 * U3) / U1
        dcov <-
            distcov(
                X,
                Y,
                affine = affine,
                bias_corr = TRUE,
                type.X = type.X,
                type.Y = type.Y,
                metr.X = metr.X,
                metr.Y = metr.Y,
                bandwidth = bandwidth,
                use = "all"
            )
        stat <- n * sign(dcov) * dcov ^ 2 + U2 * U3
        pval <- 1 - pgamma(stat, alpha, beta)
    }




    return(list("dcov" = dcov, "pval" = pval))




}




