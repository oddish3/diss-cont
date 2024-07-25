#include <RcppArmadillo.h>
#include "npiv.h" // Include the header file for npiv function
// [[Rcpp::depends(RcppArmadillo)]]

// Define the npiv_estimate function
// [[Rcpp::export]]
Rcpp::List npiv_estimate_cpp(arma::uword Ltil, const arma::mat& Px, const arma::mat& PP, const arma::mat& BB, 
                             const arma::uvec& CJ, const arma::uvec& CK, const arma::vec& y, arma::uword n) {
  
  // Rcpp::Rcout << "Debug: Initial Ltil = " << Ltil << std::endl;
  Ltil = Ltil + 1;
  // Rcpp::Rcout << "Debug: Adjusted Ltil = " << Ltil << std::endl;
  // 
  // Ensure that Ltil + 1 does not exceed the length of CJ or CK
  if (Ltil >= CJ.n_elem || Ltil >= CK.n_elem) {
    Rcpp::stop("Index Ltil out of bounds for CJ or CK.");
  }
  
  // Adjust the indices to be zero-based
  arma::uword CJ_start = CJ(Ltil - 1);
  arma::uword CJ_end = CJ(Ltil) - 1;
  arma::uword CK_start = CK(Ltil - 1);
  arma::uword CK_end = CK(Ltil) - 1;
  
  // Rcpp::Rcout << "Debug: CJ_start = " << CJ_start << ", CJ_end = " << CJ_end << std::endl;
  // Rcpp::Rcout << "Debug: CK_start = " << CK_start << ", CK_end = " << CK_end << std::endl;
  
  // Ensure that the end indices do not exceed the number of columns
  if (CJ_end >= Px.n_cols || CJ_end >= PP.n_cols || CK_end >= BB.n_cols) {
    Rcpp::stop("Index out of bounds for matrix columns.");
  }
  
  arma::mat Px1 = Px.cols(CJ_start, CJ_end);
  arma::mat PP1 = PP.cols(CJ_start, CJ_end);
  arma::mat BB1 = BB.cols(CK_start, CK_end);
  
  // Calculate and print J
  arma::uword J = Px1.n_cols;
  Rcpp::Rcout << "Value of J (number of columns in Px1): " << J << std::endl;
  
  // Debug Px1
  // Rcpp::Rcout << "Debug: Px1 dimensions = " << Px1.n_rows << "x" << Px1.n_cols << std::endl;
  // Rcpp::Rcout << "Debug: Px1 first few elements:\n" << Px1.submat(0, 0, arma::min(arma::uvec{4, Px1.n_rows-1}), arma::min(arma::uvec{4, Px1.n_cols-1})) << std::endl;
  // Rcpp::Rcout << "Debug: Px1 sum = " << arma::accu(Px1) << ", mean = " << arma::mean(arma::mean(Px1)) << std::endl;
  
  // Rcpp::Rcout << "Debug: PP1 dimensions = " << PP1.n_rows << "x" << PP1.n_cols << std::endl;
  // Rcpp::Rcout << "Debug: BB1 dimensions = " << BB1.n_rows << "x" << BB1.n_cols << std::endl;
  
  // Note: multiply Px1 by Q1 to get Lx1
  Rcpp::List npiv_result = npiv(PP1, BB1, y);
  arma::vec c1 = npiv_result["c"];
  arma::vec u1 = npiv_result["uhat"];
  arma::mat Q1 = npiv_result["Q"];
  
  // Rcpp::Rcout << "Debug: c1 length = " << c1.n_elem << std::endl;
  // Rcpp::Rcout << "Debug: u1 length = " << u1.n_elem << std::endl;
  // Rcpp::Rcout << "Debug: Q1 dimensions = " << Q1.n_rows << "x" << Q1.n_cols << std::endl;
  
  arma::mat Bu1 = BB1.each_col() % u1;
  
  // Rcpp::Rcout << "First few rows of Q1:\n" 
  //             << Q1.submat(0, 0, arma::min(arma::uvec{4, Q1.n_rows-1}), arma::min(arma::uvec{4, Q1.n_cols-1})) 
  //             << std::endl;
  // 
  // Rcpp::Rcout << "First few rows of Bu1:\n" 
  //             << Bu1.submat(0, 0, arma::min(arma::uvec{4, Bu1.n_rows-1}), arma::min(arma::uvec{4, Bu1.n_cols-1})) 
  //             << std::endl;
  // 
  // Rcpp::Rcout << "Value of n: " << n << std::endl;
  
  arma::mat OL1 = Q1 * (Bu1.t() * Bu1 / n) * Q1.t();
  
  // Rcpp::Rcout << "First few rows of OL1:\n" 
  //             << OL1.submat(0, 0, arma::min(arma::uvec{4, OL1.n_rows-1}), arma::min(arma::uvec{4, OL1.n_cols-1})) 
  //             << std::endl;
  
  
  
  // Rcpp::Rcout << "Debug: Bu1 dimensions = " << Bu1.n_rows << "x" << Bu1.n_cols << std::endl;
  // Rcpp::Rcout << "Debug: Bu1 first few elements:\n" << Bu1.submat(0, 0, arma::min(arma::uvec{4, Bu1.n_rows-1}), arma::min(arma::uvec{4, Bu1.n_cols-1})) << std::endl;
  
  // Debug OL1
  // Rcpp::Rcout << "Debug: OL1 dimensions = " << OL1.n_rows << "x" << OL1.n_cols << std::endl;
  // Rcpp::Rcout << "Debug: OL1 first few elements:\n" << OL1.submat(0, 0, arma::min(arma::uvec{4, OL1.n_rows-1}), arma::min(arma::uvec{4, OL1.n_cols-1})) << std::endl;
  // Rcpp::Rcout << "Debug: OL1 sum = " << arma::accu(OL1) << ", mean = " << arma::mean(arma::mean(OL1)) << std::endl;
  
  
  // Variance term
  arma::vec tden(Px.n_rows, arma::fill::zeros);
  
  // Debug: Print the first few rows of Px1
  // Rcpp::Rcout << "Debug: Px1 first few rows:\n" << Px1.rows(0, arma::min(arma::uvec{4, Px1.n_rows-1})) << std::endl;
  
  for (arma::uword x = 0; x < Px.n_rows; ++x) {
    arma::rowvec Px1_row = Px1.row(x);
    double s1 = arma::as_scalar(Px1_row * OL1 * Px1_row.t());
    
    // Check for negative values before sqrt
    if (s1 < 0) {
      Rcpp::Rcout << "Warning: Negative value encountered for s1 at row " << x << ": " << s1 << std::endl;
      // Handle negative value (e.g., set to 0 or small positive number)
      s1 = std::max(s1, 0.0);
    }
    
    tden(x) = std::sqrt(s1);
    
    // Print debug info for every 1000th row and the first 5 rows
    if (x % 1000 == 0 || x < 5) {
      // Rcpp::Rcout << "Debug: Row " << x << ":" << std::endl;
      // Rcpp::Rcout << "  Px1_row: " << Px1_row.subvec(0, arma::min(arma::uvec{4, Px1_row.n_elem-1})).t() << std::endl;
      // Rcpp::Rcout << "  s1: " << s1 << std::endl;
      // Rcpp::Rcout << "  tden(x): " << tden(x) << std::endl;
    }
  }
  
  // Debug: Print summary statistics for tden
  // Rcpp::Rcout << "Debug: tden summary:" << std::endl;
  // Rcpp::Rcout << "  Min: " << arma::min(tden) << std::endl;
  // Rcpp::Rcout << "  Max: " << arma::max(tden) << std::endl;
  // Rcpp::Rcout << "  Mean: " << arma::mean(tden) << std::endl;
  // Rcpp::Rcout << "  Median: " << arma::median(tden) << std::endl;
  
  arma::vec hhat = Px1 * c1;
  arma::vec sigh = tden / std::sqrt(n);
  
  // Debug: Print summary statistics for sigh
  // Rcpp::Rcout << "Debug: sigh summary:" << std::endl;
  // Rcpp::Rcout << "  Min: " << arma::min(sigh) << std::endl;
  // Rcpp::Rcout << "  Max: " << arma::max(sigh) << std::endl;
  // Rcpp::Rcout << "  Mean: " << arma::mean(sigh) << std::endl;
  // Rcpp::Rcout << "  Median: " << arma::median(sigh) << std::endl;
  // Rcpp::Rcout << "  n: " << n << std::endl;
  
  return Rcpp::List::create(Rcpp::Named("hhat") = hhat,
                            Rcpp::Named("sigh") = sigh,
                            Rcpp::Named("basis_functions") = Px1);
}