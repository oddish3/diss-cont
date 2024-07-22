#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Function to estimate function coefficient, return residuals
// [[Rcpp::export]]
Rcpp::List npiv(const arma::mat& P, const arma::mat& B, const arma::vec& y) {
  arma::mat BtB = B.t() * B;
  arma::mat PtB = P.t() * B;
  arma::mat BtP = B.t() * P;
  
  arma::mat Q = arma::pinv(PtB * arma::pinv(BtB) * BtP) * PtB * arma::pinv(BtB);
  arma::vec c = Q * B.t() * y;
  arma::vec uhat = y - P * c;
  arma::mat QQ; // Empty matrix for QQ
  
  return Rcpp::List::create(Rcpp::Named("c") = c,
                            Rcpp::Named("u") = uhat,
                            Rcpp::Named("Q") = Q * y.size());
}
