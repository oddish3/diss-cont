#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Function to compute the smallest eigenvalue and the singular value
// [[Rcpp::export]]
double shat_optimized_cpp(const arma::mat& P, const arma::mat& B) {
  arma::mat Gp = P.t() * P;
  arma::mat Gb = B.t() * B;
  arma::mat S = B.t() * P;
  
  // Compute the smallest eigenvalue of Gb
  arma::vec eigval = arma::eig_sym(Gb);
  double min_eigen_Gb = eigval(0);
  
  double s;
  if (min_eigen_Gb > 0) {
    // Compute the Cholesky decomposition
    arma::mat L_Gb = arma::chol(Gb);
    arma::mat L_Gp = arma::chol(Gp);
    
    // Solve the system using Cholesky factors
    arma::mat intermediate = arma::solve(arma::trimatl(L_Gb.t()), S);
    intermediate = arma::solve(arma::trimatl(L_Gb), intermediate);
    intermediate = arma::solve(arma::trimatu(L_Gp.t()), intermediate);
    intermediate = arma::solve(arma::trimatu(L_Gp), intermediate);
    
    // Compute the SVD of the intermediate result
    arma::vec singular_values = arma::svd(intermediate);
    s = singular_values(0);
  } else {
    s = 1e-20;
  }
  return s;
}
