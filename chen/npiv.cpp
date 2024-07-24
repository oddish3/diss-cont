#include <RcppArmadillo.h>
#include "npiv.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List npiv(const arma::mat& P, const arma::mat& B, const arma::vec& y) {
  // Rcpp::Rcout << "Debug: P dimensions: " << P.n_rows << "x" << P.n_cols << std::endl;
  // Rcpp::Rcout << "Debug: B dimensions: " << B.n_rows << "x" << B.n_cols << std::endl;
  // Rcpp::Rcout << "Debug: y dimensions: " << y.n_rows << "x" << y.n_cols << std::endl;
  
  // Compute intermediate matrices
  arma::mat BtB = B.t() * B;
  arma::mat BtB_pinv = arma::pinv(BtB);
  arma::mat PtB = P.t() * B;
  
  // Rcpp::Rcout << "Debug: BtB_pinv dimensions: " << BtB_pinv.n_rows << "x" << BtB_pinv.n_cols << std::endl;
  // Rcpp::Rcout << "Debug: PtB dimensions: " << PtB.n_rows << "x" << PtB.n_cols << std::endl;
  
  // Compute Q (corrected to match MATLAB version)
  arma::mat Q = arma::pinv(PtB * BtB_pinv * B.t() * P) * PtB * BtB_pinv;
  
  // Rcpp::Rcout << "Debug: Q dimensions: " << Q.n_rows << "x" << Q.n_cols << std::endl;
  // Rcpp::Rcout << "Debug: Q sum: " << arma::accu(Q) << ", mean: " << arma::mean(arma::vectorise(Q)) << std::endl;
  
  // Compute c
  arma::vec c = Q * B.t() * y;
  
  // Rcpp::Rcout << "Debug: c dimensions: " << c.n_rows << "x" << c.n_cols << std::endl;
  // Rcpp::Rcout << "Debug: c sum: " << arma::accu(c) << ", mean: " << arma::mean(c) << std::endl;
  
  // Compute uhat
  arma::vec uhat = y - P * c;
  
  // Rcpp::Rcout << "Debug: uhat dimensions: " << uhat.n_rows << "x" << uhat.n_cols << std::endl;
  // Rcpp::Rcout << "Debug: uhat sum: " << arma::accu(uhat) << ", mean: " << arma::mean(uhat) << std::endl;
  
  // Scale Q by length of y (to match MATLAB version)
  Q *= y.n_elem;
  
  // Rcpp::Rcout << "Debug: Scaled Q sum: " << arma::accu(Q) << ", mean: " << arma::mean(arma::vectorise(Q)) << std::endl;
  
  // Prepare the return list
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("c") = c,
    Rcpp::Named("uhat") = uhat,
    Rcpp::Named("Q") = Q // Include scaled Q in the return list
  );
  
  return result;
}