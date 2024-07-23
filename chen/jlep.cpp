#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Declare the npiv function
// [[Rcpp::export]]
Rcpp::List npiv(const arma::mat& P, const arma::mat& B, const arma::vec& y) {
  // Compute intermediate matrices
  arma::mat BtB = B.t() * B;
  arma::mat BtB_pinv = arma::pinv(BtB);
  arma::mat PtB = P.t() * B;
  
  // Compute Q
  arma::mat Q = arma::pinv(PtB * BtB_pinv * PtB.t()) * PtB * BtB_pinv;
  
  // Compute c
  arma::vec c = Q * B.t() * y;
  
  // Compute uhat
  arma::vec uhat = y - P * c;
  
  // Compute QQ
  arma::mat QQ = Q * y.n_elem;
  
  // Prepare the return list
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("c") = c,
    Rcpp::Named("uhat") = uhat,
    Rcpp::Named("Q") = Q,
    Rcpp::Named("QQ") = QQ // Include QQ in the return list
  );
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List jlep(arma::uword Lhat, const arma::mat& Px, const arma::mat& PP, const arma::mat& BB, 
                const arma::uvec& CJ, const arma::uvec& CK, const arma::vec& TJ, 
                const arma::vec& y, arma::uword n, arma::uword nb) {
  
  arma::uword Lmax = Lhat + 1;
  double Jmax = TJ(Lhat);
  arma::mat omega = arma::randn(n, nb);
  
  arma::cube ZZ(Lmax, Lmax, nb, arma::fill::zeros);
  arma::mat HH(Lmax, Lmax, arma::fill::zeros);
  
  for (arma::uword i = 0; i < Lmax - 1; ++i) {
    for (arma::uword j = i + 1; j < Lmax; ++j) {
      arma::mat Px1 = Px.cols(CJ(i), CJ(i + 1) - 1);
      arma::mat PP1 = PP.cols(CJ(i), CJ(i + 1) - 1);
      arma::mat BB1 = BB.cols(CK(i), CK(i + 1) - 1);
      arma::mat Px2 = Px.cols(CJ(j), CJ(j + 1) - 1);
      arma::mat PP2 = PP.cols(CJ(j), CJ(j + 1) - 1);
      arma::mat BB2 = BB.cols(CK(j), CK(j + 1) - 1);
      
      Rcpp::List npiv1 = npiv(PP1, BB1, y);
      Rcpp::List npiv2 = npiv(PP2, BB2, y);
      
      arma::vec c1 = npiv1["c"];
      arma::vec u1 = npiv1["uhat"];
      arma::mat Q1 = npiv1["QQ"];
      arma::vec c2 = npiv2["c"];
      arma::vec u2 = npiv2["uhat"];
      arma::mat Q2 = npiv2["QQ"];
      
      arma::mat Bu1 = BB1.each_col() % u1;
      arma::mat Bu2 = BB2.each_col() % u2;
      
      arma::mat OL1 = Q1 * (Bu1.t() * Bu1 / n) * Q1.t();
      arma::mat OL2 = Q2 * (Bu2.t() * Bu2 / n) * Q2.t();
      arma::mat OL12 = Q1 * (Bu1.t() * Bu2 / n) * Q2.t();
      
      arma::vec tden(Px.n_rows);
      for (arma::uword x = 0; x < Px.n_rows; ++x) {
        double s1 = arma::as_scalar(Px1.row(x) * OL1 * Px1.row(x).t());
        double s2 = arma::as_scalar(Px2.row(x) * OL2 * Px2.row(x).t());
        double s12 = arma::as_scalar(Px1.row(x) * OL12 * Px2.row(x).t());
        tden(x) = std::sqrt(s1 + s2 - 2 * s12);
      }
      
      arma::vec tnum = std::sqrt(n) * (Px1 * c1 - Px2 * c2);
      HH(i, j) = arma::max(arma::abs(tnum / tden));
      
      for (arma::uword b = 0; b < nb; ++b) {
        arma::vec Buw1 = Bu1.t() * omega.col(b) / std::sqrt(n);
        arma::vec Buw2 = Bu2.t() * omega.col(b) / std::sqrt(n);
        arma::vec tnum_boot = Px1 * Q1 * Buw1 - Px2 * Q2 * Buw2;
        ZZ(i, j, b) = arma::max(arma::abs(tnum_boot / tden));
      }
    }
  }
  
  arma::vec z(nb);
  for (arma::uword b = 0; b < nb; ++b) {
    z(b) = ZZ.slice(b).max();
  }
  arma::vec z_sorted = arma::sort(z);
  double quantile_index = std::max(0.5, 1.0 - std::sqrt(std::log(Jmax) / Jmax));
  double theta = z_sorted(static_cast<arma::uword>(quantile_index * (nb - 1)));
  
  
  arma::uvec indices = arma::find(arma::max(HH, 1) <= 1.1 * theta, 1);
  arma::uword LL = indices.is_empty() ? 0 : indices.at(0);
  if (LL > 0) LL -= 1;
  
  
  return Rcpp::List::create(Rcpp::Named("LL") = LL,
                            Rcpp::Named("theta") = theta);
}