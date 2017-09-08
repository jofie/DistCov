#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct Distance : public Worker {

    // input matrix to read from
    const RVector<double> vec;

    // output matrix to write to
    RMatrix<double> rmat;

    // initialize from Rcpp input and output matrixes (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
    Distance(const NumericVector vec, NumericMatrix rmat)
        : vec(vec), rmat(rmat) {}

    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            for (std::size_t j = 0; j < i; j++) {

                // rows we will operate on
                double element1 = vec[i];
                double element2 = vec[j];

                // write to output matrix
                rmat(i, j) = rmat(j, i) = fabs(element1 - element2);
            }
        }
    }
};

// [[Rcpp::export]]
NumericMatrix rcpp_parallel_distance(NumericVector vec) {

    // allocate the matrix we will return
    NumericMatrix rmat(vec.size(), vec.size());

    // create the worker
    Distance Distance(vec, rmat);

    // call it with parallelFor
    parallelFor(0, vec.size(), Distance);

    return rmat;
}

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct Colsums : public Worker {

    // input matrix to read from
    const RMatrix<double> mat;

    // output vector to write to
    RVector<double> rvec;

    // initialize from Rcpp input matrix and output vector (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
    Colsums(const NumericMatrix mat, NumericVector rvec)
        : mat(mat), rvec(rvec) {}

    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            for (unsigned int j = 0; j < mat.nrow(); j++) {

                rvec[i] += mat(j, i);
            }
        }

    }
};

// [[Rcpp::export]]
NumericVector rcpp_parallel_colsums(NumericMatrix mat) {

    // allocate the matrix we will return
    NumericVector rvec(mat.nrow());

    // create the worker
    Colsums Colsums(mat, rvec);

    // call it with parallelFor
    parallelFor(0, mat.nrow(), Colsums);

    return rvec;
}

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct Matrix_prod_sum : public Worker {

    // input matrix to read from
    const RMatrix<double> dist1;
    const RMatrix<double> dist2;

    // output vector to write to
    double rres;

    // initialize from Rcpp input matrix and output vector (the RMatrix class
    // can be automatically converted to from the Rcpp matrix type)
    Matrix_prod_sum(const NumericMatrix dist1, const NumericMatrix dist2,
                    double rres) : dist1(dist1), dist2(dist2), rres(rres) {}

    // function call operator that work for the specified range (begin/end)
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; i++) {
            for (unsigned int j = i + 1; j < dist1.nrow(); j++) {

                rres += 2 * dist1(j, i) * dist2(j, i);
            }
        }

    }
};

// [[Rcpp::export]]
double rcpp_parallel_prod_sum(const NumericMatrix dist1, const NumericMatrix dist2) {

    // allocate the matrix we will return
    double rres = 0;

    // create the worker
    Matrix_prod_sum Matrix_prod_sum(dist1, dist2, rres);

    // call it with parallelFor
    parallelFor(0, dist1.nrow(), Matrix_prod_sum);

    return rres;
}

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct Sum : public Worker
{
    // source vector
    const RVector<double> input;

    // accumulated value
    double value;

    // constructors
    Sum(const NumericVector input) : input(input), value(0) {}
    Sum(const Sum& sum, Split) : input(sum.input), value(0) {}

    // accumulate just the element of the range I've been asked to
    void operator()(std::size_t begin, std::size_t end) {
        value += std::accumulate(input.begin() + begin, input.begin() + end, 0.0);
    }

    // join my value with that of another Sum
    void join(const Sum& rhs) {
        value += rhs.value;
    }
};

// [[Rcpp::export]]
double parallelVectorSum(NumericVector x) {

    // declare the SumBody instance
    Sum sum(x);

    // call parallel_reduce to start the work
    parallelReduce(0, x.length(), sum);

    // return the computed sum
    return sum.value;
}
