#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double shat_optimized_cpp(const arma::mat& P, const arma::mat& B) {
  arma::mat Gp = P.t() * P;
  arma::mat Gb = B.t() * B;
  arma::mat S = B.t() * P;
  
  // Compute the smallest eigenvalue of Gb
  arma::vec eigval;
  bool success = arma::eig_sym(eigval, Gb);
  
  if (!success) {
    Rcpp::warning("Eigenvalue computation failed. Returning 1e-20.");
    return 1e-20;
  }
  
  double min_eigen_Gb = eigval(0);
  
  if (min_eigen_Gb > 0) {
    // Compute the Cholesky decomposition
    arma::mat L_Gb = arma::chol(Gb, "lower");
    arma::mat L_Gp = arma::chol(Gp, "lower");
    
    // Solve the system using Cholesky factors
    arma::mat intermediate = arma::solve(arma::trimatl(L_Gb), S);
    intermediate = arma::solve(arma::trimatu(L_Gp.t()), intermediate.t()).t();
    
    // Compute the largest singular value of the intermediate result
    arma::vec s = arma::svd(intermediate);
    return s(0);
  } else {
    return 1e-20;
  }
}
