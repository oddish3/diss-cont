#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Helper function to print matrix summary
void print_matrix_summary(const arma::mat& matrix, const std::string& name) {
  // Rcpp::Rcout << "Matrix " << name << " dimensions: " << matrix.n_rows << "x" << matrix.n_cols << std::endl;
  // Rcpp::Rcout << "First few elements of " << name << ":" << std::endl;
  for (uword i = 0; i < std::min(10u, matrix.n_rows); ++i) {
    for (uword j = 0; j < std::min(5u, matrix.n_cols); ++j) {
      Rcpp::Rcout << std::setprecision(10) << std::setw(15) << matrix(i, j) << " ";
    }
    Rcpp::Rcout << std::endl;
  }
  Rcpp::Rcout << std::endl;
}

// [[Rcpp::export]]
double shat(const arma::mat& P, const arma::mat& B) {
  // print_matrix_summary(P, "P");
  // print_matrix_summary(B, "B");
  
  mat Gp = P.t() * P;
  mat Gb = B.t() * B;
  mat S = B.t() * P;
  
  double s;
  
  vec eigvals = eig_sym(Gb);
  // Rcpp::Rcout << "Minimum eigenvalue of Gb: " << eigvals(0) << std::endl;
  
  if (eigvals(0) > 0) {
    mat sqrtGb = sqrtmat_sympd(Gb);
    mat sqrtGp = sqrtmat_sympd(Gp);
    mat temp = solve(sqrtGb, S) * inv(sqrtGp);
    vec ss = svd(temp);
    s = ss(ss.n_elem - 1);  // Smallest singular value
    // Rcpp::Rcout << "Smallest singular value: " << s << std::endl;
  } else {
    s = 1e-20;
    // Rcpp::Rcout << "Using default s value: " << s << std::endl;
  }
  
  return s;
}

// [[Rcpp::export]]
Rcpp::List jhat(const arma::mat& PP, const arma::mat& BB, 
                const arma::vec& CJ, const arma::vec& CK, 
                const arma::vec& TJ, double M, int n, int nL) {
  arma::vec lb(nL + 1);
  arma::vec ub(nL + 1);
  
  for (int ll = 1; ll <= nL + 1; ++ll) {
    double s;
    try {
      // Rcpp::Rcout << "Cj " << std::fixed << std::setprecision(10) << CJ(ll-1) << std::endl;
      arma::mat P_sub = PP.cols(CJ(ll-1), CJ(ll) - 1);
      arma::mat B_sub = BB.cols(CK(ll-1), CK(ll) - 1);
      
      // Rcpp::Rcout << "Iteration " << ll << ":" << std::endl;
      // print_matrix_summary(P_sub, "P_sub");
      // print_matrix_summary(B_sub, "B_sub");
      
      s = shat(P_sub, B_sub);
    } catch (std::exception& e) {
      Rcpp::Rcout << "Exception caught: " << e.what() << std::endl;
      s = 1e-20;
    }
    
    double J = TJ(ll-1);
    lb(ll-1) = J * std::sqrt(std::log(J)) * std::max(0.0, 1.0 / s);
    
    // Rcpp::Rcout << "J: " << std::fixed << std::setprecision(10) << J
    //             << ", lb(" << ll << "): " << std::setprecision(10) << lb(ll-1) << std::endl;
  }
  
  // Rcpp::Rcout << "Before setting ub:" << std::endl;
  // Rcpp::Rcout << "nL: " << nL << std::endl;
  // Rcpp::Rcout << "ub size: " << ub.size() << std::endl;
  // Rcpp::Rcout << "lb size: " << lb.size() << std::endl;
  
  ub.head(nL) = lb.subvec(1, nL);
  ub(nL) = arma::datum::inf;
  
  // Rcpp::Rcout << "After setting ub:" << std::endl;
  // Rcpp::Rcout << "ub: " << ub.t() << std::endl;
  
  double threshold = 2 * M * std::sqrt(n);
  // Rcpp::Rcout << "Threshold: " << threshold << std::endl;
  
  // Rcpp::Rcout << "lb <= threshold: " << arma::conv_to<arma::rowvec>::from(lb <= threshold) << std::endl;
  // Rcpp::Rcout << "threshold <= ub: " << arma::conv_to<arma::rowvec>::from(threshold <= ub) << std::endl;
  
  arma::uvec L = arma::find(lb <= threshold && threshold <= ub);
  
  // Rcpp::Rcout << "L after first find: " << L.t() << std::endl;
  
  int LL;
  int flag;
  
  if (L.n_elem > 0) {
    LL = L(0);
    flag = 0;
    // Rcpp::Rcout << "L.n_elem > 0, LL = " << LL << ", flag = " << flag << std::endl;
  } else {
    L = arma::find(lb <= threshold);
    // Rcpp::Rcout << "L after second find: " << L.t() << std::endl;
    if (L.n_elem > 0) {
      LL = L(L.n_elem - 1);
      flag = 1;
      // Rcpp::Rcout << "L.n_elem > 0 (second check), LL = " << LL << ", flag = " << flag << std::endl;
    } else {
      LL = 0;
      flag = 2;
      // Rcpp::Rcout << "L.n_elem == 0, LL = " << LL << ", flag = " << flag << std::endl;
    }
  }
  
  LL = std::max(LL, 1);
  // Rcpp::Rcout << "Final LL: " << LL << std::endl;
  
  return Rcpp::List::create(
    Rcpp::Named("LL") = LL,
    Rcpp::Named("flag") = flag
  );
}