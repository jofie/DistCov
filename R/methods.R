#' Calculate the distance covariance
#'
#' @param X contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".
#' @param bias_corr logical; indicates if the bias corrected version of the sample distance covariance should be calculated.
#' @param metr.X specifies the metric which should be used for X to analyse the distance covariance. TO DO: Provide details for this.
#' @param alpha.
#' @param use : "all" uses all observations, "complete.obs" only uses observations that are complete for all variables, "pairwise.complete.obs" uses pairwise complete observations for computing distance covariances
#' @return numeric giving the distance covariance between samples X and Y.
#' @export
distcovmatrix <- function(X, bias_corr = TRUE, metr.X = "euclidean", dcor=FALSE, test="none", perms=499L, alpha=NULL, use = "all") {

  #extract dimensions and sample sizes
  ss.dimX <- extract_np(X, type.X = "sample")

  n <- ss.dimX$Sample.Size
  p <- ss.dimX$Dimension

  if (is.null(alpha)) {
      alpha <- rep(1,p)
  }



  if (use == "complete.obs") {
      ccX <-  1:n
      if (type.X == "sample") {
          ccX <- which(complete.cases(X))
      }
      if (p == 1) {
          X <- X[ccX]
      } else if (type.X == "sample" && p > 1) {
          X <- X[ccX, ]
      }
      n <- length(ccX)
  } else if (use== "pairwise.complete.obs") {
    cc <- lapply(1:p, function(k) which(complete.cases(X[,k])))
  }




  if (length(metr.X == 1))
    metr.X <- rep(metr.X, p)


  distX <- cmX <- mX <- as.list(rep(NA,p))

  for (i in 1:p) {
    distX[[i]] <- distmat(X[,i], metr.X[i], alpha[i], n, p=1)
    cmX[[i]] <- colmeans(distX[[i]])
    mX[[i]] <- .Internal(mean(cmX[[i]]))
  }

  dcov2 <- dcor2 <- pval <- matrix(ncol=p,nrow=p)

  if (use== "pairwise.complete.obs") {
    if (test == "permutation" | test == "gamma") {
        for (i in 1:p) {
            for (j in i:p) {
                ccisec <- intersect(cc[[i]],cc[[j]])
                statpval <- distcov.test(distX[[i]][ccisec,ccisec], distX[[j]][ccisec,ccisec], test = test, b = perms,
                                           bias_corr = bias_corr,
                                           type.X = "distance", type.Y="distance")
                pval[i,j] <- pval[j,i] <- statpval$pval
                dcov2[i,j] <- dcov2[j,i] <- statpval$dcov^2
            }
        }
    }  else {
        for (i in 1:p) {
            for (j in i:p) {
               ccisec <- intersect(cc[[i]],cc[[j]])
               dcov2[i,j] <- dcov2[j,i] <- distcov(distX[[i]][ccisec,ccisec], distX[[j]][ccisec,ccisec],
                                         bias_corr = bias_corr,
                                         type.X = "distance", type.Y="distance")
            }
        }
    }
  } else {
    if (test == "permutation" | test == "gamma") {
        for (i in 1:p) {
            for (j in i:p) {
                statpval <- distcov.test(distX[[i]], distX[[j]], test = test, b = perms, affine = FALSE,
                                          bias_corr = bias_corr,
                                          type.X = "distance", type.Y="distance")
                pval[i,j] <- pval[j,i] <- statpval$pval
                dcov2[i,j] <- dcov2[j,i] <- statpval$dcov^2
            }
        }
   }   else {
       for (i in 1:p) {
           for (j in i:p) {
                if (bias_corr == TRUE) {
                term1 <- matrix_prod_sum(distX[[i]], distX[[j]]) / n / (n - 3)
                term2 <- n ^ 4 * mX[[i]] * mX[[j]] / n / (n - 1) / (n - 2) / (n - 3)
                term3 <- n ^ 2 * vector_prod_sum(cmX[[i]],  cmX[[j]]) / n / (n - 2) / (n - 3)
                dcov2[i,j] <- dcov2[j,i] <- term1 + term2 - 2 * term3
                }
               if (bias_corr == FALSE) {
                   term1 <- matrix_prod_sum(distX[[i]], distX[[j]]) / n ^ 2
                   term2 <- mX[[i]] * mX[[j]]
                   term3 <- vector_prod_sum(cmX[[i]], cmY[[j]]) / n
                   dcov2[i,j] <- dcov2[j,i] <- term1 + term2 - 2 * term3
               }
           }
       }

      }

    }


  if (dcor==TRUE) {
    for (i in 1:p) {
      for (j in 1:p) {
        dcor2[i,j] <- dcov2[i,j] / sqrt((dcov2[i,i] * dcov2[j,j]))
      }
    }
  }

  diag(pval) <- NA

  return(list("dcovmatrix"=dcov2, "dcormatrix"=dcor2, "p.values"=pval))
}

#' Calculate the distance covariance
#'
#' @param X contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".
#' @param bias_corr logical; indicates if the bias corrected version of the sample distance covariance should be calculated.
#' @param metr.X specifies the metric which should be used for X to analyse the distance covariance. TO DO: Provide details for this.
#' @param alpha.
#' @param use : "all" uses all observations, "complete.obs" only uses observations that are complete for all variables, "pairwise.complete.obs" uses pairwise complete observations for computing distance covariances
#' @return numeric giving the distance covariance between samples X and Y.
#' @export
dcsis <- function(response, X, bias_corr = TRUE, metr.response = "euclidean", metr.X, alpha=c(1,1), use="all", top = NULL, threshold = NULL, return.all = FALSE) {

    ss.dimX <- extract_np(X, type.X = "sample")

    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension

    distresp <- distmat(response, metr.response, alpha[1], n, p=1)

    dcorvec <- sapply(1:p, function(k) distcorr(X[,k], distresp, affine = FALSE, bias_corr =bias_corr, type.X = "sample",
                                                type.Y = "distance", metr.X = "euclidean", metr.Y = "response",
                                                alpha = alpha, use = use))

    orddcor <- Rfast::Order(dcorvec, descending = TRUE)

    if (!is.null(top) | is.null(threshold))
        top <- ceiling(n/log(n))
        else if (is.null(top))
        top <- length(which(dcorvec > threshold))

    select <- orddcor[1:top]

    if (return.all == FALSE)
    return(select)
    else
    return(list("order.all" = orddcor, "dcor.all" = dcorvec, "selected" = select, "dcor.selected"=dcorvec[select]))

}

