#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector DyadUpdate_new(IntegerVector & Y, NumericVector & C)
{
  /* variable declaration */
  int k, pos;
  int n = Y.size();
  int L = ceil(log(n) / log(2));
  NumericVector GAMMA(n), S(pow(2, L + 1));

  for (int ii = 1; ii <= n - 1; ii++) {
    for (int ell = 0; ell <= L - 1; ell++) {
      double k_db = Y[ii - 1] / pow(2.0, ell);
      k = ceil(k_db) - 1;
      pos = k;
      if (ell > 0) {
        for (int scale = ell - 1; scale >= 0; scale--) {
          pos += pow(2, L - scale);
        }
      }
      S[pos] += C[ii - 1];
    }
    for (int ell = 0; ell <= L - 1; ell++) {
      double k_db = (Y[ii] - 1) / pow(2.0, ell);
      k = floor(k_db) - 1;
      if ((k + 1) / 2.0 > floor((k + 1) / 2.0)) {
        pos = k;
        if (ell > 0) {
          for (int scale = ell - 1; scale >= 0; scale--) {
            pos += pow(2, L - scale);
          }
        }
        GAMMA[ii] += S[pos];
      }
    }
  }
  return GAMMA;
}

// [[Rcpp::export]]
NumericVector PartialSum2D_new(NumericVector & X,
                           NumericVector & Y,
                           NumericVector & Z,
                           IntegerVector & IX0,
                           IntegerVector & IY0) {
  int n = X.size();
  double Zdot = 0;
  NumericVector Y_sort(n), Z_sort(n), sX(n + 1), sY(n + 1), gamma(n), gamma1(n);
  IntegerVector IX(n), IY(n);

  for (int i = 0; i < n; i++) {
      int id_x = IX0[i];
      int id_y = IY0[i];
    IX[id_x] = i;
    IY[id_y] = i;
    Y_sort[i] = Y[id_x];
    Z_sort[i] = Z[id_x];
    sX[i + 1] = sX[i] + Z_sort[i];
    sY[i + 1] = sY[id_y] + Z_sort[id_y];
    Zdot += Z_sort[i];
  }
  gamma1 = DyadUpdate_new(IY, Z_sort);
  for (int i = 0; i < n; i++) {
      int id_x = IX[i];
      int id_y = IY[id_x];
    gamma[i] = Zdot - Z_sort[id_x] - 2 * sY[id_y] - 2 * sX[id_x] + 4 * gamma1[id_x];
  }
  return gamma;
}

// // [[Rcpp::export]]
// distcov_fast <- function(X, Y) {
//   n <- length(X)
//   temp <- IX <- IY <- 1:n
//
//   IX0 <- Rfast::Order(X)
//   vX <- X[IX0]
//   IX[IX0] <- temp
//
//   IY0 <- Rfast::Order(Y)
//   vY <- Y[IY0]
//   IY[IY0] <- temp
//
//   sX <- cumsum(vX)
//   sY <- cumsum(vY)
//
//   alphaX <- IX - 1
//   alphaY <- IY - 1
//   betaX <- sX[IX] - vX[IX]
//   betaY <- sY[IY] - vY[IY]
//
//   Xdot <- sum(X)
//   Ydot <- sum(Y)
//
//   aidot <- Xdot + (2 * alphaX - n) * X - 2 * betaX
//   bidot <- Ydot + (2 * alphaY - n) * Y - 2 * betaY
//   Sab <- vector_prod_sum(aidot, bidot) # sum(aidot * bidot)
//
//   adotdot <- 2 * vector_prod_sum(alphaX, X) - 2 * sum(betaX)
//   bdotdot <- 2 * vector_prod_sum(alphaY, Y) - 2 * sum(betaY)
//
//   gamma_1  <- PartialSum2D(X, Y, rep(1, n))
//   gamma_X  <- PartialSum2D(X, Y, X)
//   gamma_Y  <- PartialSum2D(X, Y, Y)
//   gamma_XY <- PartialSum2D(X, Y, X * Y)
//
//   aijbij <- specific_vector_prod_sum(X, Y, gamma_1, gamma_X, gamma_Y, gamma_XY)
//   dCov <- aijbij / n / (n - 3) - 2 * Sab / n / (n - 2) / (n - 3) + adotdot * bdotdot / n / (n - 1) / (n - 2) / (n - 3)
//   return(dCov)
// }
//
//
// distvar_fast <- function(X) {
//   n <- length(X)
//   temp <- IX <- 1:n
//
//   IX0 <- Rfast::Order(X)
//   vX <- X[IX0]
//   IX[IX0] <- temp
//
//
//   sX <- cumsum(vX)
//   alphaX <- IX - 1
//   betaX <- sX[IX] - vX[IX]
//   Xdot <- sum(X)
//
//   aidot <- Xdot + (2 * alphaX - n) * X - 2 * betaX
//   Saa <- vector_prod_sum(aidot, aidot)
//
//   adotdot <- 2 * vector_prod_sum(alphaX, X) - 2 * sum(betaX)
//
//   gamma_1  <- PartialSum2D(X, X, rep(1, n))
//   gamma_X  <- PartialSum2D(X, X, X)
//   gamma_XX <- PartialSum2D(X, X, X * X)
//
//   aijaij <- specific_vector_prod_sum(X, X, gamma_1, gamma_X, gamma_X, gamma_XX)
//   dVar <- aijaij / n / (n - 3) - 2 * Saa / n / (n - 2) / (n - 3) + adotdot * adotdot / n / (n - 1) / (n - 2) / (n - 3)
//   return(dVar)
// }

