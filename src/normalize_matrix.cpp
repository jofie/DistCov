#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix normalize_matrix(NumericMatrix & M, NumericVector & cmM, double & mM) {
  //NumericMatrix res(clone(M));
  NumericVector m(M.nrow(), mM);
  unsigned int ncol = M.ncol();
  for (unsigned int n = 0; n < ncol; ++n) {
    M.row(n) = M.row(n) - cmM;
    M.column(n) = M.column(n) - cmM + m;
  }
  return M;
}

// [[Rcpp::export]]
NumericMatrix hadamard_product(NumericMatrix & X, NumericMatrix & Y){
  unsigned int ncol = X.ncol();
  unsigned int nrow = X.nrow();
  int counter = 0;
  for (unsigned int j=0; j<ncol; j++) {
    for (unsigned int i=0; i<nrow; i++)  {
      X[counter++] *= Y(i, j);
    }
  }
  return X;
}




