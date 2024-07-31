#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

using namespace Rcpp;

// [[Rcpp::export]]
List bspline(NumericVector x, int l, int r) {
  int N = x.size();
  int m = std::pow(2, l) - 1;
  r = r + 1;
  
  // Define the augmented knot set
  std::vector<double> kts;
  if (l == 0) {
    kts.resize(2 * (r - 1));
    std::fill(kts.begin(), kts.begin() + (r - 1), 0.0);
    std::fill(kts.begin() + (r - 1), kts.end(), 1.0);
  } else if (l >= 1) {
    kts.resize(2 * (r - 2) + m + 2);
    std::fill(kts.begin(), kts.begin() + (r - 2), 0.0);
    for (int i = 0; i <= m; ++i) {
      kts[r - 2 + i] = i / static_cast<double>(std::pow(2, l));
    }
    std::fill(kts.begin() + r - 2 + m + 1, kts.end(), 1.0);
  }
  
  // Initialize for recursion
  std::vector<std::vector<std::vector<double>>> BB(N, std::vector<std::vector<double>>(m + 2 * r - 2, std::vector<double>(r - 1, 0.0)));
  for (int i = 0; i < N; ++i) {
    for (int j = r - 1; j <= r + m - 1; ++j) {
      if (x[i] >= kts[j] && x[i] < kts[j + 1]) {
        BB[i][j - (r - 1)][0] = 1.0;
        break;
      }
    }
  }
  
  // Recursion
  for (int j = 1; j < r - 1; ++j) {
    for (int i = 0; i < m + 2 * r - 2 - j; ++i) {
      std::vector<double> a1(N, 0.0);
      if (kts[i + j] - kts[i] != 0) {
        for (int n = 0; n < N; ++n) {
          a1[n] = (x[n] - kts[i]) / (kts[i + j] - kts[i]);
        }
      }
      
      std::vector<double> a2(N, 0.0);
      if (kts[i + j + 1] - kts[i + 1] != 0) {
        for (int n = 0; n < N; ++n) {
          a2[n] = (x[n] - kts[i + 1]) / (kts[i + j + 1] - kts[i + 1]);
        }
      }
      
      for (int n = 0; n < N; ++n) {
        BB[n][i][j] = a1[n] * BB[n][i][j - 1] + (1 - a2[n]) * BB[n][i + 1][j - 1];
      }
    }
  }
  
  NumericMatrix XX(N, std::pow(2, l) + r - 2);
  for (int n = 0; n < N; ++n) {
    for (int i = 0; i < std::pow(2, l) + r - 2; ++i) {
      XX(n, i) = BB[n][i][r - 2 - 1]; // adjust indexing to match MATLAB's (r-2) indexing
    }
  }
  
  // Calculate derivative
  NumericMatrix DX(N, m + r - 1);
  for (int i = 0; i < m + r - 1; ++i) {
    std::vector<double> a1(N, 0.0);
    if (kts[i + r - 2] - kts[i] != 0) {
      std::fill(a1.begin(), a1.end(), 1.0 / (kts[i + r - 2] - kts[i]));
    }
    
    std::vector<double> a2(N, 0.0);
    if (kts[i + r - 1] - kts[i + 1] != 0) {
      std::fill(a2.begin(), a2.end(), 1.0 / (kts[i + r - 1] - kts[i + 1]));
    }
    
    if (i < m + r - 2) {
      for (int n = 0; n < N; ++n) {
        DX(n, i) = (r - 2) * (a1[n] * BB[n][i][r - 3] - a2[n] * BB[n][i + 1][r - 3]);
      }
    } else {
      for (int n = 0; n < N; ++n) {
        DX(n, i) = (r - 2) * (a1[n] * BB[n][i][r - 3]);
      }
    }
  }
  
  return List::create(Named("XX") = XX, Named("DX") = DX);
}
