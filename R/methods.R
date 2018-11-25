#' Calculate the distance covariance
#'
#' @param X contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".
#' @param bias.corr logical; indicates if the bias corrected version of the sample distance covariance should be calculated.
#' @param type.X either "sample" or "distance"; specifies the type of input for X.
#' @param metr.X specifies the metric which should be used for X to analyse the distance covariance. TO DO: Provide details for this.
#' @param use : "all" uses all observations, "complete.obs" only uses observations that are complete for all variables, "pairwise.complete.obs" uses pairwise complete observations for computing distance covariances
#' @return numeric giving the distance covariance between samples X and Y.
#' @export
distcovmatrix <-
  function(X,
           bias.corr = TRUE,
           type.X = "sample",
           metr.X = "euclidean",
           dcor = FALSE,
           test = "none",
           perms = 499L,
           use = "all") {
    #extract dimensions and sample sizes
    if (is.data.frame(X))
      X <- as.matrix(X)
    
    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension
    
   
  if (use == "complete.obs") {

      ccX <-  1:n
      if (type.X == "sample") {
        ccX <- which(complete.cases(X))
      }
      if (p == 1) {
        X <- X[ccX]
      } else if (type.X == "sample" && p > 1) {
        X <- X[ccX,]
      }
      n <- length(ccX)
    } else if (use == "pairwise.complete.obs") {
      cc <- lapply(1:p, function(k)
        which(complete.cases(X[, k])))
    }
    
    
    
    
    if (length(metr.X == 1))
      metr.X <- rep(metr.X, p)
    
    
    distX <- cmX <- mX <- as.list(rep(NA, p))
    
    for (i in 1:p) {
      distX[[i]] <- distmat(X = X[, i], 
                            metr.X = metr.X[i], 
                            n = n, 
                            p = 1)
      cmX[[i]] <- colmeans(distX[[i]])
      mX[[i]] <- .Internal(mean(cmX[[i]]))
    }
    
    dcov2 <- dcor2 <- pval <- matrix(ncol = p, nrow = p)
    
    if (use == "pairwise.complete.obs") {
      if (test == "permutation" | test == "gamma") {
        for (i in 1:p) {
          for (j in i:p) {
            ccisec <- intersect(cc[[i]], cc[[j]])
            statpval <-
              distcov.test(
                X = distX[[i]][ccisec, ccisec],
                Y = distX[[j]][ccisec, ccisec],
                test = test,
                b = perms,
                affine = FALSE,
                bias.corr = bias.corr,
                type.X = "distance",
                type.Y = "distance"
              )
            pval[i, j] <- pval[j, i] <- statpval$pval
            dcov2[i, j] <- dcov2[j, i] <- statpval$dcov ^ 2
          }
        }
      }  else {
        for (i in 1:p) {
          for (j in i:p) {
            ccisec <- intersect(cc[[i]], cc[[j]])
            dcov2[i, j] <-
              dcov2[j, i] <-
              distcov(
                X = distX[[i]][ccisec, ccisec],
                Y = distX[[j]][ccisec, ccisec],
                bias.corr = bias.corr,
                type.X = "distance",
                type.Y = "distance"
              )
          }
        }
      }
    } else {
      if (test == "permutation" | test == "gamma") {
        for (i in 1:p) {
          for (j in i:p) {
            statpval <-
              distcov.test(
                X = distX[[i]],
                Y = distX[[j]],
                test = test,
                b = perms,
                affine = FALSE,
                bias.corr = bias.corr,
                type.X = "distance",
                type.Y = "distance"
              )
            pval[i, j] <- pval[j, i] <- statpval$pval
            dcov2[i, j] <- dcov2[j, i] <- statpval$dcov ^ 2
          }
        }
      }   else {
        for (i in 1:p) {
          for (j in i:p) {
            if (bias.corr == TRUE) {
              term1 <- matrix_prod_sum(distX[[i]], distX[[j]]) / n / (n - 3)
              term2 <-
                n ^ 4 * mX[[i]] * mX[[j]] / n / (n - 1) / (n - 2) / (n - 3)
              term3 <-
                n ^ 2 * vector_prod_sum(cmX[[i]],  cmX[[j]]) / n / (n - 2) / (n - 3)
              dcov2[i, j] <-
                dcov2[j, i] <- term1 + term2 - 2 * term3
            }
            if (bias.corr == FALSE) {
              term1 <- matrix_prod_sum(distX[[i]], distX[[j]]) / n ^ 2
              term2 <- mX[[i]] * mX[[j]]
              term3 <- vector_prod_sum(cmX[[i]], cmY[[j]]) / n
              dcov2[i, j] <-
                dcov2[j, i] <- term1 + term2 - 2 * term3
            }
          }
        }
        
      }
      
    }
    
    
    if (dcor == TRUE) {
      for (i in 1:p) {
        for (j in 1:p) {
          dcor2[i, j] <- dcov2[i, j] / sqrt((dcov2[i, i] * dcov2[j, j]))
        }
      }
    }
    
    diag(pval) <- NA
    
    return(list(
      "dcovmatrix" = dcov2,
      "dcormatrix" = dcor2,
      "p.values" = pval
    ))
  }

#' Calculate the distance covariance
#'
#' @param X contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".
#' @param bias.corr logical; indicates if the bias corrected version of the sample distance covariance should be calculated.
#' @param metr.X specifies the metric which should be used for X to analyse the distance covariance. TO DO: Provide details for this.
#' @param use : "all" uses all observations, "complete.obs" only uses observations that are complete for all variables, "pairwise.complete.obs" uses pairwise complete observations for computing distance covariances
#' @return numeric giving the distance covariance between samples X and Y.
#' @export
dcsis <-
  function(response,
           X,
           bias.corr = TRUE,
           metr.response = "euclidean",
           metr.X,
           use = "all",
           top = NULL,
           threshold = NULL,
           return.all = FALSE) {
    ss.dimX <- extract_np(X = X, type.X = "sample")
    
    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension
    
    distresp <- distmat(X = response, 
                        metr.X = metr.response, 
                        n = n, 
                        p = 1)
    
    dcorvec <-
      sapply(1:p, function(k)
        distcorr(
          X = X[, k],
          Y = distresp,
          affine = FALSE,
          bias.corr = bias.corr,
          type.X = "sample",
          type.Y = "distance",
          metr.X = "euclidean",
          metr.Y = "response",
          use = use
        ))
    
    if  (use == "complete.obs") {
            cc_resp <- complete.cases(response)
            response <- response[cc_resp]
            X <- X[cc_resp,]
    }



    orddcor <- Rfast::Order(dcorvec, descending = TRUE)
    
    if (!is.null(top) | is.null(threshold))
      top <- min(p, ceiling(n / log(n)))
    else if (is.null(top))
      top <- length(which(dcorvec > threshold))

    selected <- orddcor[1:top]
    
    if (return.all == FALSE)
      return(selected)
    else
      return(
        list(
          "order.all" = orddcor,
          "dcor.all" = dcorvec,
          "selected" = selected,
          "dcor.selected" = dcorvec[selected]
        )
      )
    
  }
