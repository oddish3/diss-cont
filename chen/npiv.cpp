// chen/npiv.cpp
#include <RcppArmadillo.h>
#include "npiv.h"
// [[Rcpp::depends(RcppArmadillo)]]

// Add export attribute to the npiv function
// [[Rcpp::export]]
Rcpp::List npiv(const arma::mat& P, const arma::mat& B, const arma::vec& y) {
  // Compute intermediate matrices
  arma::mat BtB = B.t() * B;
  arma::mat BtB_pinv = arma::pinv(BtB);
  arma::mat PtB = P.t() * B;
  
  // Compute Q
  arma::mat Q = arma::pinv(PtB * BtB_pinv * PtB.t()) * PtB * BtB_pinv;
  
  // Compute c
  arma::vec c = Q * B.t() * y;
  
  // Compute uhat
  arma::vec uhat = y - P * c;
  
  // Prepare the return list
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("c") = c,
    Rcpp::Named("uhat") = uhat,
    Rcpp::Named("Q") = Q // Include Q in the return list
  );
  
  return result;
}
