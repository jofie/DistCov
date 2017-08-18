#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix normalize_matrix(NumericMatrix M, NumericVector cmM, double mM) {
  NumericVector m(M.nrow(), mM);
  for (int n = 0; n < M.nrow(); ++n) {
    M.row(n) = M.row(n) - cmM;
    M.column(n) = M.column(n) - cmM + m;
  }
  return M;
}

