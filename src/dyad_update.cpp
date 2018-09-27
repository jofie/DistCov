#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector DyadUpdate(NumericVector & Y, NumericVector & C)
{
  /* variable declaration */
  int k, pos;
  int n = Y.size();
  int L = ceil(log(n) / log(2));
  NumericVector GAMMA(n), S(pow(2, L + 1));
  
  for (int ii = 1; ii <= n - 1; ii++) {
    for (int ell = 0; ell <= L - 1; ell++) {
      double k_db = Y[ii - 1] / pow(2, ell);
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
      double k_db = (Y[ii] - 1) / pow(2, ell);
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


