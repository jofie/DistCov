#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix normalize_matrix(NumericMatrix & M, NumericVector & cmM, double & mM) {
  //NumericMatrix res(clone(M));
  NumericVector m(M.nrow(), mM);
  unsigned int ncol = M.ncol();
  for (unsigned int n = 0; n < ncol; ++n) {
    M.row(n) = M.row(n) - cmM;
    M.column(n) = M.column(n) - cmM + mM;
  }
  return M;
}



// [[Rcpp::export]]
NumericMatrix hadamard_product(NumericMatrix & X, NumericMatrix & Y){
  unsigned int n = X.size();
  unsigned int counter1 = 0;
  unsigned int counter2 = 0;
  for (unsigned int j=0; j < n; j++) {
      X[counter1++] *= Y[counter2++];
  }
  return X;
}

// [[Rcpp::export]]
NumericVector vector_product(NumericVector & X, NumericVector & Y){
  unsigned int n = X.size();
  int counter = 0;
  for (unsigned int j=0; j<n; j++) {
    X[counter++] *= Y[j];
  }
  return X;
}

// [[Rcpp::export]]
double matrix_prod_sum(const NumericMatrix X, const NumericMatrix Y){
    unsigned int n = X.nrow();
    double res = 0;
    for (unsigned int j = 0; j < n; j++) {
      for (unsigned int i = (j+1); i < n; i++) {
       res += X(i, j) * Y(i, j);
      }
  }
  return 2 * res;
}

// [[Rcpp::export]]
double vector_prod_sum(const NumericVector X, const NumericVector Y){
    unsigned int n = X.size();
    double res = 0;
    for (unsigned int j = 0; j < n; j++) {
        res += X[j] * Y[j];
    }
    return res;
}


// [[Rcpp::export]]
double matrix_prod_sum_sample(const NumericMatrix X, const NumericMatrix Y,const NumericVector s){
    unsigned int n = X.nrow();
    double res = 0;
    for (unsigned int j = 0; j < n; j++) {
        for (unsigned int i = (j+1); i < n; i++) {
            res += X(i, j) * Y(s[i], s[j]);
        }
    }
    return 2 * res;
}




