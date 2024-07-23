// chen_ucb_cv.cpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <numeric>
// Assuming npiv.h is available and properly included
#include "npiv.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

double quantile(const vector<double>& data, double alpha) {
  vector<double> sorted_data = data;
  sort(sorted_data.begin(), sorted_data.end());
  size_t idx = static_cast<size_t>(alpha * sorted_data.size());
  return sorted_data[idx];
}

// [[Rcpp::export]]
double ucb_cv(const int Ltil, const int Lhat, const mat& Px, const mat& PP, const mat& BB, 
              const ivec& CJ, const ivec& CK, const vec& y, const int n, const int nb, 
              const int type, const double alpha) {
  
  int Lmax = Lhat + 1;
  mat omega = randn<mat>(n, nb);
  mat ZZ = zeros<mat>(max(Ltil, Lmax - 2), nb);
  
  for (int i = 0; i < max(Ltil, Lmax - 2); ++i) {
    
    // Precompute that which can be pre-computed
    mat Px1 = Px.cols(CJ[i], CJ[i + 1] - 1);
    mat PP1 = PP.cols(CJ[i], CJ[i + 1] - 1);
    mat BB1 = BB.cols(CK[i], CK[i + 1] - 1);
    
    // Note: multiply Px1 by Q1 to get Lx1
    List npiv_result = npiv(PP1, BB1, y);
    vec u1 = npiv_result["uhat"];
    mat Q1 = npiv_result["Q"];
    
    mat Bu1 = BB1.each_col() % u1;
    
    mat OL1 = Q1 * (Bu1.t() * Bu1 / n) * Q1.t();
    
    // Variance term
    vec tden = vec(Px.n_rows, fill::zeros);
    for (int x = 0; x < Px.n_rows; ++x) {
      double s1 = as_scalar(Px1.row(x) * OL1 * Px1.row(x).t());
      tden[x] = sqrt(s1);
    }
    
    // Bootstrap
    for (int b = 0; b < nb; ++b) {
      vec Buw1 = Bu1.t() * omega.col(b) / sqrt(n);
      
      // Compute bootstrapped sup-t-stat at (J, J2)
      vec tnum = Px1 * Q1 * Buw1;
      if (type == 0) {
        ZZ(i, b) = max(abs(tnum / tden));
      } else if (type == -1) {
        ZZ(i, b) = max(tnum / tden);
      } else if (type == 1) {
        ZZ(i, b) = min(tnum / tden);
      }
    }
  }
  
  // Critical value
  vector<double> z;
  if (type == 0 || type == -1) {
    for (int b = 0; b < nb; ++b) {
      z.push_back(ZZ.col(b).max());
    }
    return quantile(z, 1 - alpha);
  } else if (type == 1) {
    for (int b = 0; b < nb; ++b) {
      z.push_back(ZZ.col(b).min());
    }
    return -quantile(z, alpha);
  }
  return 0.0;
}
