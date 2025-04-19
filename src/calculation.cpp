// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP MtxProd(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP expit(const Eigen::MatrixXd& x) {
  Eigen::MatrixXd result = (1 / (1 + (-x.array()).exp()));
  return Rcpp::wrap(result);
}
