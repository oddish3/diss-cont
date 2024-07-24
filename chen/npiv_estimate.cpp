#include <RcppArmadillo.h>
#include "npiv.h" // Include the header file
// [[Rcpp::depends(RcppArmadillo)]]

// Define the npiv_estimate function
// [[Rcpp::export]]
Rcpp::List npiv_estimate_cpp(arma::uword Ltil, const arma::mat& Px, const arma::mat& PP, const arma::mat& BB, 
                             const arma::uvec& CJ, const arma::uvec& CK, const arma::vec& y, arma::uword n) {
  
  Rcpp::Rcout << "Debug: Initial Ltil = " << Ltil << std::endl;
  Ltil = Ltil + 1;
  Rcpp::Rcout << "Debug: Adjusted Ltil = " << Ltil << std::endl;
  
  // Ensure that Ltil + 1 does not exceed the length of CJ or CK
  if (Ltil >= CJ.n_elem || Ltil >= CK.n_elem) {
    Rcpp::stop("Index Ltil out of bounds for CJ or CK.");
  }
  
  // Adjust the indices to be zero-based
  arma::uword CJ_start = CJ(Ltil - 1);
  arma::uword CJ_end = CJ(Ltil) - 1;
  arma::uword CK_start = CK(Ltil - 1);
  arma::uword CK_end = CK(Ltil) - 1;
  
  Rcpp::Rcout << "Debug: CJ_start = " << CJ_start << ", CJ_end = " << CJ_end << std::endl;
  Rcpp::Rcout << "Debug: CK_start = " << CK_start << ", CK_end = " << CK_end << std::endl;
  
  // Ensure that the end indices do not exceed the number of columns
  if (CJ_end >= Px.n_cols || CJ_end >= PP.n_cols || CK_end >= BB.n_cols) {
    Rcpp::stop("Index out of bounds for matrix columns.");
  }
  
  arma::mat Px1 = Px.cols(CJ_start, CJ_end);
  arma::mat PP1 = PP.cols(CJ_start, CJ_end);
  arma::mat BB1 = BB.cols(CK_start, CK_end);
  
  Rcpp::Rcout << "Debug: Px1 dimensions = " << Px1.n_rows << "x" << Px1.n_cols << std::endl;
  Rcpp::Rcout << "Debug: PP1 dimensions = " << PP1.n_rows << "x" << PP1.n_cols << std::endl;
  Rcpp::Rcout << "Debug: BB1 dimensions = " << BB1.n_rows << "x" << BB1.n_cols << std::endl;
  
  // Note: multiply Px1 by Q1 to get Lx1
  Rcpp::List npiv_result = npiv(PP1, BB1, y);
  arma::vec c1 = npiv_result["c"];
  arma::vec u1 = npiv_result["uhat"];
  arma::mat Q1 = npiv_result["Q"];
  
  Rcpp::Rcout << "Debug: c1 length = " << c1.n_elem << std::endl;
  Rcpp::Rcout << "Debug: u1 length = " << u1.n_elem << std::endl;
  Rcpp::Rcout << "Debug: Q1 dimensions = " << Q1.n_rows << "x" << Q1.n_cols << std::endl;
  
  arma::mat Bu1 = BB1.each_col() % u1;
  arma::mat OL1 = Q1 * (Bu1.t() * Bu1 / n) * Q1.t();
  
  Rcpp::Rcout << "Debug: Bu1 dimensions = " << Bu1.n_rows << "x" << Bu1.n_cols << std::endl;
  Rcpp::Rcout << "Debug: OL1 dimensions = " << OL1.n_rows << "x" << OL1.n_cols << std::endl;
  
  // Variance term
  arma::vec tden(Px.n_rows, arma::fill::zeros);
  for (arma::uword x = 0; x < Px.n_rows; ++x) {
    double s1 = arma::as_scalar(Px1.row(x) * OL1 * Px1.row(x).t());
    tden(x) = std::sqrt(s1);
    if (x % 100 == 0) {
      Rcpp::Rcout << "Debug: Processing row " << x << ", s1 = " << s1 << ", tden(x) = " << tden(x) << std::endl;
    }
  }
  
  arma::vec hhat = Px1 * c1;
  arma::vec sigh = tden / std::sqrt(n);
  
  Rcpp::Rcout << "Debug: hhat length = " << hhat.n_elem << std::endl;
  Rcpp::Rcout << "Debug: sigh length = " << sigh.n_elem << std::endl;
  
  return Rcpp::List::create(Rcpp::Named("hhat") = hhat,
                            Rcpp::Named("sigh") = sigh);
}