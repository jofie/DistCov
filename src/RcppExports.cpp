// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// DyadUpdate
NumericVector DyadUpdate(NumericVector& Y, NumericVector& C);
RcppExport SEXP _DistCov_DyadUpdate(SEXP YSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(DyadUpdate(Y, C));
    return rcpp_result_gen;
END_RCPP
}
// normalize_matrix
NumericMatrix normalize_matrix(NumericMatrix& M, NumericVector& cmM, double& mM);
RcppExport SEXP _DistCov_normalize_matrix(SEXP MSEXP, SEXP cmMSEXP, SEXP mMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type M(MSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type cmM(cmMSEXP);
    Rcpp::traits::input_parameter< double& >::type mM(mMSEXP);
    rcpp_result_gen = Rcpp::wrap(normalize_matrix(M, cmM, mM));
    return rcpp_result_gen;
END_RCPP
}
// hadamard_product
NumericMatrix hadamard_product(NumericMatrix& X, NumericMatrix& Y);
RcppExport SEXP _DistCov_hadamard_product(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(hadamard_product(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// vector_product
NumericVector vector_product(NumericVector& X, NumericVector& Y);
RcppExport SEXP _DistCov_vector_product(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(vector_product(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// matrix_prod_sum
double matrix_prod_sum(NumericMatrix& X, NumericMatrix& Y);
RcppExport SEXP _DistCov_matrix_prod_sum(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_prod_sum(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// vector_prod_sum
double vector_prod_sum(NumericVector& X, NumericVector& Y);
RcppExport SEXP _DistCov_vector_prod_sum(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(vector_prod_sum(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// specific_vector_prod_sum
double specific_vector_prod_sum(NumericVector& X, NumericVector& Y, NumericVector& gamma_1, NumericVector& gamma_X, NumericVector& gamma_Y, NumericVector& gamma_XY);
RcppExport SEXP _DistCov_specific_vector_prod_sum(SEXP XSEXP, SEXP YSEXP, SEXP gamma_1SEXP, SEXP gamma_XSEXP, SEXP gamma_YSEXP, SEXP gamma_XYSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gamma_1(gamma_1SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gamma_X(gamma_XSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gamma_Y(gamma_YSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type gamma_XY(gamma_XYSEXP);
    rcpp_result_gen = Rcpp::wrap(specific_vector_prod_sum(X, Y, gamma_1, gamma_X, gamma_Y, gamma_XY));
    return rcpp_result_gen;
END_RCPP
}
// matrix_sum
double matrix_sum(NumericMatrix& X);
RcppExport SEXP _DistCov_matrix_sum(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_sum(X));
    return rcpp_result_gen;
END_RCPP
}
// vector_prod_sum_sample
double vector_prod_sum_sample(const NumericVector X, const NumericVector Y, const IntegerVector s);
RcppExport SEXP _DistCov_vector_prod_sum_sample(SEXP XSEXP, SEXP YSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(vector_prod_sum_sample(X, Y, s));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DistCov_DyadUpdate", (DL_FUNC) &_DistCov_DyadUpdate, 2},
    {"_DistCov_normalize_matrix", (DL_FUNC) &_DistCov_normalize_matrix, 3},
    {"_DistCov_hadamard_product", (DL_FUNC) &_DistCov_hadamard_product, 2},
    {"_DistCov_vector_product", (DL_FUNC) &_DistCov_vector_product, 2},
    {"_DistCov_matrix_prod_sum", (DL_FUNC) &_DistCov_matrix_prod_sum, 2},
    {"_DistCov_vector_prod_sum", (DL_FUNC) &_DistCov_vector_prod_sum, 2},
    {"_DistCov_specific_vector_prod_sum", (DL_FUNC) &_DistCov_specific_vector_prod_sum, 6},
    {"_DistCov_matrix_sum", (DL_FUNC) &_DistCov_matrix_sum, 1},
    {"_DistCov_vector_prod_sum_sample", (DL_FUNC) &_DistCov_vector_prod_sum_sample, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_DistCov(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
