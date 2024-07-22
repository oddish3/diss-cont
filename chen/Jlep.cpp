#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
Rcpp::List npiv(const arma::mat& P, const arma::mat& B, const arma::vec& y); // Declaration of the npiv function

// [[Rcpp::export]]
Rcpp::List Jlep(int Lhat, const arma::mat& Px, const arma::mat& PP, const arma::mat& BB, 
                const arma::ivec& CJ, const arma::ivec& CK, const arma::ivec& TJ, 
                const arma::vec& y, int n, int nb) {
  
  int Lmax = Lhat + 1;
  int Jmax = TJ[Lhat + 1];
  arma::mat omega = arma::randn(n, nb);
  
  arma::cube ZZ(Lmax, Lmax, nb, arma::fill::zeros);
  arma::mat HH(Lmax, Lmax, arma::fill::zeros);
  
  // Pre-compute npiv results
  std::vector<Rcpp::List> npiv_results(Lmax);
  for (int idx = 0; idx < Lmax; ++idx) {
    arma::mat PP_sub = PP.cols(CJ[idx], CJ[idx + 1] - 1);
    arma::mat BB_sub = BB.cols(CK[idx], CK[idx + 1] - 1);
    npiv_results[idx] = npiv(PP_sub, BB_sub, y);
  }
  
  // Main computation
  for (int i = 0; i < Lmax - 1; ++i) {
    for (int j = i + 1; j < Lmax; ++j) {
      arma::mat Px1 = Px.cols(CJ[i], CJ[i + 1] - 1);
      arma::mat Px2 = Px.cols(CJ[j], CJ[j + 1] - 1);
      
      Rcpp::List npiv1 = npiv_results[i];
      Rcpp::List npiv2 = npiv_results[j];
      
      arma::vec c1 = npiv1["c"];
      arma::vec u1 = npiv1["u"];
      arma::mat Q1 = npiv1["Q"];
      arma::vec c2 = npiv2["c"];
      arma::vec u2 = npiv2["u"];
      arma::mat Q2 = npiv2["Q"];
      
      arma::mat BB1 = BB.cols(CK[i], CK[i + 1] - 1);
      arma::mat BB2 = BB.cols(CK[j], CK[j + 1] - 1);
      
      arma::mat Bu1 = BB1.each_col() % u1;
      arma::mat Bu2 = BB2.each_col() % u2;
      
      arma::mat OL1 = Q1 * (Bu1.t() * Bu1 / n) * Q1.t();
      arma::mat OL2 = Q2 * (Bu2.t() * Bu2 / n) * Q2.t();
      arma::mat OL12 = Q1 * (Bu1.t() * Bu2 / n) * Q2.t();
      
      arma::vec tden = arma::sqrt(arma::sum((Px1 * OL1) % Px1, 1) + arma::sum((Px2 * OL2) % Px2, 1) - 2 * arma::sum((Px1 * OL12) % Px2, 1));
      arma::vec tnum = arma::sqrt(n) * (Px1 * c1 - Px2 * c2);
      HH(i, j) = arma::max(arma::abs(tnum / tden));
      
      for (int b = 0; b < nb; ++b) {
        arma::vec Buw1 = Bu1.t() * omega.col(b) / sqrt(n);
        arma::vec Buw2 = Bu2.t() * omega.col(b) / sqrt(n);
        arma::vec tnum_boot = Px1 * Q1 * Buw1 - Px2 * Q2 * Buw2;
        ZZ(i, j, b) = arma::max(arma::abs(tnum_boot / tden));
      }
    }
  }
  
  // Critical value
  arma::vec z = arma::max(ZZ, 2);
  double theta = arma::quantile(z, std::max(0.5, 1 - sqrt(log(Jmax) / Jmax)));
  
  // Step 2: cut-off rule
  arma::uvec LL_indices = arma::find(arma::max(HH, 1) <= 1.1 * theta);
  int LL = LL_indices.is_empty() ? -1 : LL_indices[0] - 1;
  
  return Rcpp::List::create(Rcpp::Named("LL") = LL,
                            Rcpp::Named("theta") = theta);
}
