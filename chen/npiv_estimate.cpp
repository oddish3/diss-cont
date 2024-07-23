#include <RcppArmadillo.h>
#include "npiv.h" // Include the header file
// [[Rcpp::depends(RcppArmadillo)]]

// Define the npiv_estimate function
// [[Rcpp::export]]
Rcpp::List npiv_estimate_cpp(arma::uword Ltil, const arma::mat& Px, const arma::mat& PP, const arma::mat& BB, 
                             const arma::uvec& CJ, const arma::uvec& CK, const arma::vec& y, arma::uword n) {
  
  Ltil = Ltil + 1;
  
  // Ensure that Ltil + 1 does not exceed the length of CJ or CK
  if (Ltil >= CJ.n_elem || Ltil >= CK.n_elem) {
    Rcpp::stop("Index Ltil out of bounds for CJ or CK.");
  }
  
  // Adjust the indices to be zero-based
  arma::uword CJ_start = CJ(Ltil - 1);
  arma::uword CJ_end = CJ(Ltil) - 1;
  arma::uword CK_start = CK(Ltil - 1);
  arma::uword CK_end = CK(Ltil) - 1;
  
  // Ensure that the end indices do not exceed the number of columns
  if (CJ_end >= Px.n_cols || CJ_end >= PP.n_cols || CK_end >= BB.n_cols) {
    Rcpp::stop("Index out of bounds for matrix columns.");
  }
  
  arma::mat Px1 = Px.cols(CJ_start, CJ_end);
  arma::mat PP1 = PP.cols(CJ_start, CJ_end);
  arma::mat BB1 = BB.cols(CK_start, CK_end);
  
  // Note: multiply Px1 by Q1 to get Lx1
  Rcpp::List npiv_result = npiv(PP1, BB1, y);
  arma::vec c1 = npiv_result["c"];
  arma::vec u1 = npiv_result["uhat"];
  arma::mat Q1 = npiv_result["Q"];
  
  arma::mat Bu1 = BB1.each_col() % u1;
  arma::mat OL1 = Q1 * (Bu1.t() * Bu1 / n) * Q1.t();
  
  // Variance term
  arma::vec tden(Px.n_rows, arma::fill::zeros);
  for (arma::uword x = 0; x < Px.n_rows; ++x) {
    double s1 = arma::as_scalar(Px1.row(x) * OL1 * Px1.row(x).t());
    tden(x) = std::sqrt(s1);
  }
  
  arma::vec hhat = Px1 * c1;
  arma::vec sigh = tden / std::sqrt(n);
  
  return Rcpp::List::create(Rcpp::Named("hhat") = hhat,
                            Rcpp::Named("sigh") = sigh);
}
