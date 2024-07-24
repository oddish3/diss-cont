#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double shat(const arma::mat& P, const arma::mat& B) {
  mat Gp = P.t() * P;
  mat Gb = B.t() * B;
  mat S = B.t() * P;
  
  double s;
  
  vec eigvals = eig_sym(Gb);
  if (eigvals(0) > 0) {
    mat sqrtGb = sqrtmat_sympd(Gb);
    mat sqrtGp = sqrtmat_sympd(Gp);
    mat temp = inv(sqrtGb) * S * inv(sqrtGp);
    vec ss = svd(temp);
    s = ss(ss.n_elem - 1);  // Smallest singular value
  } else {
    s = 1e-20;
  }
  
  return s;
}