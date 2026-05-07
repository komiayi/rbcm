#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector log_lik_med_arma(vec params,
                                     vec M1,
                                     vec M2,
                                     mat Z1,
                                     mat Z2,
                                     vec T_val,
                                     uvec idx0,
                                     uvec idx1,
                                     std::string methods = "CC",
                                     bool return_sum = true) {
  
  int n  = T_val.n_elem;
  int p1 = Z1.n_cols;
  int p2 = Z2.n_cols;
  
  vec s1_vec(n), s2_vec(n), rho_vec(n);
  
  // 1. Coefficients
  vec alpha1 = params.subvec(0, p1 - 1);
  vec alpha2 = params.subvec(p1, p1 + p2 - 1);
  
  // 2. Résidus (vectorisé)
  vec res1 = M1 - Z1 * alpha1;
  vec res2 = M2 - Z2 * alpha2;
  
  vec f(n);
  if(methods == "CC") {
    
    double s1  = std::exp(params[p1 + p2]);
    double s2  = std::exp(params[p1 + p2 + 1]);
    double rho = std::tanh(params[p1 + p2 + 2]);
    
    // calcul vectorisé
    vec term = square(res1)/ (s1*s1)
      + square(res2)/ (s2*s2)
      - 2*rho*(res1 % res2)/(s1*s2);
      
      term = term / (1 - rho*rho);
      
      f = -std::log(2*M_PI)
        - std::log(s1)
        - std::log(s2)
        - 0.5*std::log(1 - rho*rho)
        - 0.5 * term;
        
  } else {
    
    double s1_1 = std::exp(params[p1 + p2]);
    double s1_0 = std::exp(params[p1 + p2 + 1]);
    double s2_1 = std::exp(params[p1 + p2 + 2]);
    double s2_0 = std::exp(params[p1 + p2 + 3]);
    
    double rho0 = std::tanh(params[p1 + p2 + 4]);
    double rho1 = std::tanh(params[p1 + p2 + 5]);
    
    // ⚠️ idx venant de R → convertir en 0-based
    uvec id1 = idx1 - 1;
    uvec id0 = idx0 - 1;
    
    s1_vec.elem(id1).fill(s1_1);
    s1_vec.elem(id0).fill(s1_0);
    
    s2_vec.elem(id1).fill(s2_1);
    s2_vec.elem(id0).fill(s2_0);
    
    rho_vec.elem(id1).fill(rho1);
    rho_vec.elem(id0).fill(rho0);
    
    vec am1 = res1 / s1_vec;
    vec am2 = res2 / s2_vec;
    
    vec A40 = square(am1) + square(am2) - 2 * rho_vec % am1 % am2;
    vec A41 = A40 / (1 - square(rho_vec));
    
    f = -std::log(2*M_PI)
      - log(s1_vec)
      - log(s2_vec)
      - 0.5 * log(1 - square(rho_vec))
      - 0.5 * A41;
  }
  if(return_sum){
    return Rcpp::NumericVector::create(-sum(f)); 
  } else {
    return  Rcpp::wrap(-f);       
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector log_lik_outc_arma(vec params,
                                      vec Y,
                                      mat Y1,
                                      bool return_sum = true) {
  
  int p = params.n_elem;
  
  // coefficients beta
  vec beta = params.subvec(0, p - 2);
  
  // résidus
  vec res = Y - Y1 * beta;
  
  // écart-type
  double s = std::exp(params[p - 1]);
  
  // log-vraisemblance vectorisée
  vec f = -0.5 * std::log(2 * M_PI)
    - std::log(s)
    - 0.5 * square(res) / (s * s);
    
    if(return_sum){
      return Rcpp::NumericVector::create(-sum(f)); 
    } else {
      return  Rcpp::wrap(-f);       
    }
}


// [[Rcpp::export]]
Rcpp::List get_initials_arma(arma::mat Z, arma::vec M) {
  arma::vec coeffs = arma::solve(Z, M); // Résolution QR/least squares
  coeffs.replace(arma::datum::nan, 0.0); // Remplace NA par 0
  
  arma::vec residuals = M - Z * coeffs;
  
  return Rcpp::List::create(
    Rcpp::Named("alphas") = coeffs,
    Rcpp::Named("res") = residuals
  );
}
