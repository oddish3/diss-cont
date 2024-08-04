#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List acr_inference(const vec& x, const vec& y, const mat& Px, const mat& Dx, const mat& PP, 
                   const vec& hhat, const vec& dhat, const vec& Xx_sub) {
  int n = x.n_elem;
  int n_basis = Px.n_cols;
  
  // Calculate E_n[ΔY|D=0]
  vec D0_indices = find(x == 0);
  double E_DY_D0 = mean(y(D0_indices));
  
  // Calculate ATE_K^(D)
  vec ATE_K_D = Px * PP.i() * PP.t() * y;
  
  // Identify positive doses
  uvec positive_indices = find(x > 0);
  vec positive_doses = x(positive_indices);
  int n_positive = positive_indices.n_elem;
  
  // Calculate ACR^o
  vec ACR_K_D = dhat(round(interp1(Xx_sub, regspace(0, Xx_sub.n_elem - 1), positive_doses)));
  double ACR_o = mean(ACR_K_D);
  
  // Calculate η_acr^o(W_i)
  mat psi_K_D = Px.rows(round(interp1(Xx_sub, regspace(0, Xx_sub.n_elem - 1), positive_doses)));
  mat d_psi_K_D = Dx.rows(round(interp1(Xx_sub, regspace(0, Xx_sub.n_elem - 1), positive_doses)));
  
  vec u_i = y(positive_indices) - E_DY_D0 - ATE_K_D(positive_indices);
  
  mat E_psi_K_D_psi_K_D_inv = inv_sympd(psi_K_D.t() * psi_K_D / n_positive);
  vec E_d_psi_K_D = mean(d_psi_K_D, 0).t();
  
  vec eta = ACR_K_D - ACR_o + 
    (E_d_psi_K_D.t() * E_psi_K_D_psi_K_D_inv * psi_K_D.t() * u_i / n_positive);
  
  // Calculate σ^2_ACR^o
  double sigma_sq_ACR_o = mean(square(eta));
  
  // Standard error for ACR^o
  double se_ACR_o = sqrt(sigma_sq_ACR_o / n_positive);
  
  // Calculate confidence intervals (95%)
  double z_score = R::qnorm(0.975, 0.0, 1.0);
  double ci_lower = ACR_o - z_score * se_ACR_o;
  double ci_upper = ACR_o + z_score * se_ACR_o;
  
  return List::create(
    Named("ACR_o") = ACR_o,
    Named("se_ACR_o") = se_ACR_o,
    Named("ci_lower") = ci_lower,
    Named("ci_upper") = ci_upper,
    Named("sigma_sq_ACR_o") = sigma_sq_ACR_o
  );
}