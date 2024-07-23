#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <random>
#include <algorithm>
#include "npiv.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// Quantile function
vec quantile(const mat& v, const vec& q) {
  vec sorted_v = sort(v);
  vec result(q.n_elem);
  for (uword i = 0; i < q.n_elem; ++i) {
    double idx = q(i) * (v.n_elem - 1);
    uword lower = floor(idx);
    uword upper = ceil(idx);
    double weight = idx - lower;
    result(i) = (1 - weight) * sorted_v(lower) + weight * sorted_v(upper);
  }
  return result;
}

// [[Rcpp::export]]
vec ucb_cc(int L, mat Px, mat PP, mat BB, ivec CJ, ivec CK, vec y, int n, int nb, int type, vec alpha) {
  // Random number generation
  mat omega = randn(n, nb);
  
  // Step 1: compute critical value
  vec z(nb);
  int i = L;
  
  // Ensure L is valid
  if (i >= CJ.size() - 1) {
    stop("Invalid value for L: L exceeds or equals CJ size - 1.");
  }
  
  // Check if indices are valid
  if (CJ(i) >= Px.n_cols || CJ(i + 1) > Px.n_cols) {
    stop("Index out of bounds: CJ values are too large for Px columns.");
  }
  
  if (CK(i) >= BB.n_cols || CK(i + 1) > BB.n_cols) {
    stop("Index out of bounds: CK values are too large for BB columns.");
  }
  
  mat Px1 = Px.cols(CJ(i), CJ(i + 1) - 1);
  mat PP1 = PP.cols(CJ(i), CJ(i + 1) - 1);
  mat BB1 = BB.cols(CK(i), CK(i + 1) - 1);
  
  Rcpp::List npiv_result;
  try {
    npiv_result = npiv(PP1, BB1, y);
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  
  // Check if 'uhat' and 'Q' exist in the npiv_result
  if (!npiv_result.containsElementNamed("uhat") || !npiv_result.containsElementNamed("Q")) {
    stop("npiv function did not return expected elements 'uhat' and 'Q'");
  }
  
  vec u1 = as<vec>(npiv_result["uhat"]);
  mat Q1 = as<mat>(npiv_result["Q"]);
  
  mat Bu1 = BB1.each_col() % u1;
  
  mat OL1 = Q1 * (Bu1.t() * Bu1 / n) * Q1.t();
  
  // Variance term
  vec tden(Px.n_rows);
  for (uword x = 0; x < Px.n_rows; ++x) {
    double s1 = as_scalar(Px1.row(x) * OL1 * Px1.row(x).t());
    tden(x) = sqrt(s1);
  }
  
  // Bootstrap
  for (int b = 0; b < nb; ++b) {
    vec Buw1 = Bu1.t() * omega.col(b) / sqrt(n);
    
    // Compute bootstrapped sup-t-stat at (J, J2)
    vec tnum = Px1 * Q1 * Buw1;
    if (type == 0) {
      z(b) = max(abs(tnum / tden));
    } else if (type == -1) {
      z(b) = max(tnum / tden);
    } else if (type == 1) {
      z(b) = min(tnum / tden);
    }
  }
  
  // Critical value
  vec cv;
  mat z_matrix = conv_to<mat>::from(z);  // Convert z to a matrix
  
  if (type == 0 || type == -1) {
    cv = quantile(z_matrix, 1 - alpha);
  } else if (type == 1) {
    cv = -quantile(z_matrix, alpha);
  } else {
    stop("Invalid type value. Must be -1, 0, or 1.");
  }
  
  return cv;
}
