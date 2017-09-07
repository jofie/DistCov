#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector DyadUpdate(NumericVector & Y, NumericVector & C){
  int n = Y.size();
  NumericVector GAMMA(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      if (Y[j] < Y[i]) {
      	GAMMA[i] += C[j];
      }
    }
  }
  return GAMMA;
}




