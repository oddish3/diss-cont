// chen_ucb_cv.cpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <numeric>
#include <iostream>
// Assuming npiv.h is available and properly included
#include "npiv.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

// Debug macro
#define DEBUG_PRINT(x) std::cout << #x << " = " << x << std::endl

double quantile(const vector<double>& data, double alpha) {
  if (data.empty()) {
    std::cerr << "Error: Empty data in quantile function" << std::endl;
    return NA_REAL;
  }
  vector<double> sorted_data = data;
  sort(sorted_data.begin(), sorted_data.end());
  size_t idx = static_cast<size_t>(alpha * (sorted_data.size() - 1));
  return sorted_data[idx];
}

// [[Rcpp::export]]
vec ucb_cv(const int Ltil, const int Lhat, const mat& Px, const mat& PP, const mat& BB, 
           const ivec& CJ, const ivec& CK, const vec& y, const int n, const int nb, 
           const int type, const vec& alpha) {
  
  DEBUG_PRINT(Ltil);
  DEBUG_PRINT(Lhat);
  DEBUG_PRINT(n);
  DEBUG_PRINT(nb);
  DEBUG_PRINT(type);
  DEBUG_PRINT(alpha);
  
  int Lmax = Lhat + 1;
  mat omega = randn<mat>(n, nb);
  mat ZZ = zeros<mat>(max(Ltil, Lmax - 2), nb);
  
  DEBUG_PRINT(Lmax);
  DEBUG_PRINT(omega.n_rows);
  DEBUG_PRINT(omega.n_cols);
  DEBUG_PRINT(ZZ.n_rows);
  DEBUG_PRINT(ZZ.n_cols);
  
  for (int i = 0; i < max(Ltil, Lmax - 2); ++i) {
    DEBUG_PRINT(i);
    
    // Check if indices are within bounds
    if (i + 1 >= CJ.n_elem || i + 1 >= CK.n_elem) {
      std::cerr << "Error: Index out of bounds for CJ or CK" << std::endl;
      return vec(alpha.n_elem, fill::value(NA_REAL));
    }
    
    // Precompute that which can be pre-computed
    uword start_col_J = CJ(i);
    uword end_col_J = CJ(i + 1) - 1;
    uword start_col_K = CK(i);
    uword end_col_K = CK(i + 1) - 1;
    
    // DEBUG_PRINT(start_col_J);
    // DEBUG_PRINT(end_col_J);
    // DEBUG_PRINT(start_col_K);
    // DEBUG_PRINT(end_col_K);
    
    // Check if column indices are within bounds
    if (start_col_J >= Px.n_cols || end_col_J >= Px.n_cols || 
        start_col_J >= PP.n_cols || end_col_J >= PP.n_cols || 
        start_col_K >= BB.n_cols || end_col_K >= BB.n_cols) {
      std::cerr << "Error: Column indices out of bounds for Px, PP, or BB" << std::endl;
      return vec(alpha.n_elem, fill::value(NA_REAL));
    }
    
    mat Px1 = Px.cols(start_col_J, end_col_J);
    mat PP1 = PP.cols(start_col_J, end_col_J);
    mat BB1 = BB.cols(start_col_K, end_col_K);
    
    // DEBUG_PRINT(Px1.n_rows);
    // DEBUG_PRINT(Px1.n_cols);
    // DEBUG_PRINT(PP1.n_rows);
    // DEBUG_PRINT(PP1.n_cols);
    // DEBUG_PRINT(BB1.n_rows);
    // DEBUG_PRINT(BB1.n_cols);
    
    // Note: multiply Px1 by Q1 to get Lx1
    List npiv_result = npiv(PP1, BB1, y);
    vec u1 = npiv_result["uhat"];
    mat Q1 = npiv_result["Q"];
    
    // DEBUG_PRINT(u1.n_elem);
    // DEBUG_PRINT(Q1.n_rows);
    // DEBUG_PRINT(Q1.n_cols);
    
    mat Bu1 = BB1.each_col() % u1;
    
    mat OL1 = Q1 * (Bu1.t() * Bu1 / n) * Q1.t();
    
    // DEBUG_PRINT(OL1.n_rows);
    // DEBUG_PRINT(OL1.n_cols);
    
    // Variance term
    vec tden = vec(Px.n_rows, fill::zeros);
    for (uword x = 0; x < Px.n_rows; ++x) {
      double s1 = as_scalar(Px1.row(x) * OL1 * Px1.row(x).t());
      tden(x) = sqrt(s1);
    }
    
    // DEBUG_PRINT(tden.n_elem);
    
    // Bootstrap
    for (int b = 0; b < nb; ++b) {
      vec Buw1 = Bu1.t() * omega.col(b) / sqrt(n);
      
      // Compute bootstrapped sup-t-stat at (J, J2)
      vec tnum = Px1 * Q1 * Buw1;
      vec ratio = tnum / tden;
      
      if (type == 0) {
        ZZ(i, b) = max(abs(ratio));
      } else if (type == -1) {
        ZZ(i, b) = max(ratio);
      } else if (type == 1) {
        ZZ(i, b) = min(ratio);
      }
    }
  }
  
  // Critical value
  vector<double> z;
  for (int b = 0; b < nb; ++b) {
    if (type == 0 || type == -1) {
      z.push_back(max(ZZ.col(b)));
    } else if (type == 1) {
      z.push_back(min(ZZ.col(b)));
    }
  }
  
  // Calculate critical values for each alpha
  vec cv(alpha.n_elem);
  for (uword i = 0; i < alpha.n_elem; ++i) {
    if (type == 0 || type == -1) {
      cv(i) = quantile(z, 1 - alpha(i));
    } else if (type == 1) {
      cv(i) = -quantile(z, alpha(i));
    }
    // DEBUG_PRINT(cv(i));
  }
  
  return cv;
}