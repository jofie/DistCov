#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix add_vector(NumericMatrix & X, NumericVector & Y){
  unsigned int nrow = X.nrow();
  for (unsigned int j=0; j<nrow; j++) {
    X(j, 0) = X(j, 0) + Y[0];
    X(j, 1) = X(j, 1) + Y[1];
  }
  return X;
}

// [[Rcpp::export]]
NumericMatrix match_coords(NumericMatrix X, NumericMatrix Y){
  unsigned int nrow = X.nrow();
  int counter = 0;
  NumericMatrix res(nrow, 2);
  for (unsigned int j = 0; j < nrow; j++) {
    for (unsigned int i = 0; i < nrow; i++) {
      if (X(j, 0) == Y(i, 0) && X(j, 1) == Y(i, 1)) {
	int idx = counter++;
  	res(idx, 0) = j + 1;
  	res(idx, 1) = i + 1;
  	break;
      }
    }
  }
  return res;
}
