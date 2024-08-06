#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Custom interpolation function
vec custom_interp1(const vec& x, const vec& y, const vec& x_new) {
  vec y_new(x_new.n_elem);
  for (uword i = 0; i < x_new.n_elem; ++i) {
    uword idx = as_scalar(find(x <= x_new(i), 1, "last"));
    if (idx == x.n_elem - 1) {
      y_new(i) = y(idx);
    } else {
      double t = (x_new(i) - x(idx)) / (x(idx + 1) - x(idx));
      y_new(i) = y(idx) + t * (y(idx + 1) - y(idx));
    }
  }
  return y_new;
}

// [[Rcpp::export]]
List acr_inference(const vec& x, const vec& y, const mat& Px, const mat& Dx, const mat& PP, 
                   const vec& hhat, const vec& dhat, const vec& Xx_sub) {
  Rcout << "Function started. Input dimensions:" << std::endl;
  Rcout << "x: " << x.n_elem << ", y: " << y.n_elem << std::endl;
  Rcout << "Px: " << Px.n_rows << "x" << Px.n_cols << std::endl;
  Rcout << "Dx: " << Dx.n_rows << "x" << Dx.n_cols << std::endl;
  Rcout << "PP: " << PP.n_rows << "x" << PP.n_cols << std::endl;
  Rcout << "hhat: " << hhat.n_elem << ", dhat: " << dhat.n_elem << ", Xx_sub: " << Xx_sub.n_elem << std::endl;
  
  int n = x.n_elem;
  int n_basis = Px.n_cols;
  
  // Calculate E_n[ΔY|D=0]
  uvec D0_indices = find(x == 0);
  double E_DY_D0 = mean(y.elem(D0_indices));
  Rcout << "E_DY_D0 calculated: " << E_DY_D0 << std::endl;
  
  // Calculate ATE_K^(D)
  Rcout << "Calculating ATE_K_D. Dimensions:" << std::endl;
  Rcout << "Px: " << Px.n_rows << "x" << Px.n_cols << std::endl;
  Rcout << "PP: " << PP.n_rows << "x" << PP.n_cols << std::endl;
  Rcout << "y: " << y.n_elem << std::endl;
  
  vec ATE_K_D = Px * pinv(PP) * PP.t() * y;
  Rcout << "ATE_K_D calculated. Dimension: " << ATE_K_D.n_elem << std::endl;
  
  // Identify positive doses
  uvec positive_indices = find(x > 0);
  vec positive_doses = x.elem(positive_indices);
  int n_positive = positive_indices.n_elem;
  Rcout << "Positive doses identified. Count: " << n_positive << std::endl;
  
  // Calculate ACR^o
  vec ACR_K_D = dhat.elem(conv_to<uvec>::from(round(custom_interp1(Xx_sub, regspace<vec>(0, Xx_sub.n_elem - 1), positive_doses))));
  double ACR_o = mean(ACR_K_D);
  Rcout << "ACR_o calculated: " << ACR_o << std::endl;
  
  // Calculate η_acr^o(W_i)
  Rcout << "Calculating η_acr^o(W_i). Dimensions:" << std::endl;
  mat psi_K_D = Px.rows(conv_to<uvec>::from(round(custom_interp1(Xx_sub, regspace<vec>(0, Xx_sub.n_elem - 1), positive_doses))));
  Rcout << "psi_K_D: " << psi_K_D.n_rows << "x" << psi_K_D.n_cols << std::endl;
  
  mat d_psi_K_D = Dx.rows(conv_to<uvec>::from(round(custom_interp1(Xx_sub, regspace<vec>(0, Xx_sub.n_elem - 1), positive_doses))));
  Rcout << "d_psi_K_D: " << d_psi_K_D.n_rows << "x" << d_psi_K_D.n_cols << std::endl;
  
  vec u_i = y.elem(positive_indices) - E_DY_D0 - ATE_K_D.elem(positive_indices);
  Rcout << "u_i calculated. Dimension: " << u_i.n_elem << std::endl;
  
  mat E_psi_K_D_psi_K_D_pinv = pinv(psi_K_D.t() * psi_K_D / n_positive);
  Rcout << "E_psi_K_D_psi_K_D_pinv: " << E_psi_K_D_psi_K_D_pinv.n_rows << "x" << E_psi_K_D_psi_K_D_pinv.n_cols << std::endl;
  
  vec E_d_psi_K_D = mean(d_psi_K_D, 0).t();
  Rcout << "E_d_psi_K_D: " << E_d_psi_K_D.n_elem << std::endl;
  
  vec eta = ACR_K_D - ACR_o + 
    (E_d_psi_K_D.t() * E_psi_K_D_psi_K_D_pinv * psi_K_D.t() * u_i / n_positive);
  Rcout << "eta calculated. Dimension: " << eta.n_elem << std::endl;
  
  // Calculate σ^2_ACR^o
  double sigma_sq_ACR_o = mean(square(eta));
  Rcout << "sigma_sq_ACR_o calculated: " << sigma_sq_ACR_o << std::endl;
  
  // Standard error for ACR^o
  double se_ACR_o = sqrt(sigma_sq_ACR_o / n_positive);
  Rcout << "se_ACR_o calculated: " << se_ACR_o << std::endl;
  
  // Calculate confidence intervals (95%)
  double z_score = R::qnorm(0.975, 0.0, 1.0, 1, 0);
  double ci_lower = ACR_o - z_score * se_ACR_o;
  double ci_upper = ACR_o + z_score * se_ACR_o;
  Rcout << "Confidence interval calculated: [" << ci_lower << ", " << ci_upper << "]" << std::endl;
  
  return List::create(
    Named("ACR_o") = ACR_o,
    Named("se_ACR_o") = se_ACR_o,
    Named("ci_lower") = ci_lower,
    Named("ci_upper") = ci_upper,
    Named("sigma_sq_ACR_o") = sigma_sq_ACR_o
  );
}