#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Include the shat function we optimized earlier
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

// [[Rcpp::export]]
Rcpp::List jhat(const arma::mat& PP, const arma::mat& BB, 
                const arma::uvec& CJ, const arma::uvec& CK, 
                const arma::vec& TJ, double M, int n, int nL) {
  arma::vec lb(nL + 1);
  arma::vec ub(nL + 1);
  
  for (int ll = 0; ll < nL + 1; ++ll) {
    double s;
    try {
      arma::mat P_sub = PP.cols(CJ(ll), CJ(ll+1) - 1);
      arma::mat B_sub = BB.cols(CK(ll), CK(ll+1) - 1);
      s = shat_optimized_cpp(P_sub, B_sub);
      
      // Print the s value
      Rcpp::Rcout << "s: " << s << std::endl;
      
      // Print the head of P_sub
      Rcpp::Rcout << "Head of P_sub:" << std::endl;
      for (size_t i = 0; i < std::min((size_t)5, (size_t)P_sub.n_rows); ++i) {
        for (size_t j = 0; j < std::min((size_t)5, (size_t)P_sub.n_cols); ++j) {
          Rcpp::Rcout << P_sub(i, j) << " ";
        }
        Rcpp::Rcout << std::endl;
      }
      
      // Print the head of B_sub
      Rcpp::Rcout << "Head of B_sub:" << std::endl;
      for (size_t i = 0; i < std::min((size_t)5, (size_t)B_sub.n_rows); ++i) {
        for (size_t j = 0; j < std::min((size_t)5, (size_t)B_sub.n_cols); ++j) {
          Rcpp::Rcout << B_sub(i, j) << " ";
        }
        Rcpp::Rcout << std::endl;
      }
      
    } catch (...) {
      s = 1e-20;
    }
    
    
    double J = TJ(ll);
    lb(ll) = J * std::sqrt(std::log(J)) * std::max(0.0, 1.0 / s);
    
    // Print the lower bound calculation
    Rcpp::Rcout << "J: " << J << ", lb(" << ll << "): " << lb(ll) << std::endl;
  }
  
  ub.head(nL) = lb.subvec(1, nL);
  ub(nL) = arma::datum::inf; // Set the last element to infinity
  
  // Print the upper bounds
  Rcpp::Rcout << "ub: " << ub << std::endl;
  
  double threshold = 2 * M * std::sqrt(n);
  arma::uvec L = arma::find(lb <= threshold && threshold <= ub);
  
  // Print the threshold and indices found
  // Rcpp::Rcout << "threshold: " << threshold << ", L: " << L << std::endl;
  
  int LL;
  int flag;
  
  if (L.n_elem > 0) {
    LL = L(0);
    flag = 0;
  } else {
    L = arma::find(lb <= threshold);
    if (L.n_elem > 0) {
      LL = L(L.n_elem - 1);
      flag = 1;
    } else {
      LL = 0;
      flag = 2;
    }
  }
  
  LL = std::max(LL, 0);
  
  // Print final results
  Rcpp::Rcout << "LL: " << LL << ", flag: " << flag << std::endl;
  
  return Rcpp::List::create(
    Rcpp::Named("LL") = LL,
    Rcpp::Named("flag") = flag
  );
}
