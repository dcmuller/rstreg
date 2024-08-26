#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat weibull_gr(const arma::vec& theta, const arma::mat& X, Nullable<arma::mat> Z,
               const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
               const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec& offset) {

  int n = d.n_elem;       // Number of observations
  int pX = X.n_cols;      // Number of parameters for X
  int pZ = Z.isNotNull() ? as<mat>(Z).n_cols : 0; // Number of parameters for Z (if applicable)

  // Check that the length of theta matches the sum of columns of X and Z
  if (theta.n_elem != (pX + pZ)) {
    stop("Length of theta must be equal to the sum of the number of columns in X and Z.");
  }

  // Extract parameters from theta
  vec beta = theta.head(pX);
  vec gamma = pZ > 0 ? theta.tail(pZ) : vec();


  // Compute lambda and log_lambda
  vec log_lambda = X * beta;
  vec lambda = exp(log_lambda);

  // Compute log_p and p
  vec log_p;
  if (pfixed.isNotNull()) {
    log_p = log(as<vec>(pfixed));
  } else {
    if (pZ > 0) {
      log_p = as<mat>(Z) * gamma;
    } else {
      stop("Z matrix must be provided when pfixed is NULL.");
    }
  }
  vec p = exp(log_p);

  // Gradient computation
  mat gr(n, theta.n_elem, fill::zeros);

  // Gradient for beta terms
  for (int i = 0; i < n; ++i) {
    double lambda_i = lambda[i];
    double tt_i = tt[i];
    double tt0_i = tt0[i];
    double w_i = w[i];
    double p_i = p[i];
    double d_i = d[i];
    vec X_i = X.row(i).t();
    double tmp1=0;
    if (tt0_i > 0) {
      tmp1 = exp(p_i * log(tt0_i));
    }
    double tmp2 = exp(p_i * log(tt_i));
    double tmpval1 = w_i *(d_i + lambda_i * (tmp1 - tmp2));
    vec tmpvec1(X_i.n_elem, fill::zeros);
      for(int j=0; j<tmpvec1.n_elem; ++j)
        tmpvec1[j] += tmpval1;

    gr.row(i).head(pX) = (X_i % tmpvec1).t();

    if (pZ > 0) {
      vec Z_i = as<mat>(Z).row(i).t();
      double tmp1 = 0;
      if (tt0_i > 0)
        tmp1 = exp(p_i * log(tt0_i)) * log(tt0_i);
      double tmp2 = exp(p_i * log(tt_i)) * log(tt_i);
      double tmpval2 = w_i*(d_i*(1.0 + p_i*log(tt_i)) + lambda_i * p_i * (tmp1 - tmp2));
      vec tmpvec2(Z_i.n_elem, fill::zeros);
      for(int j=0; j<tmpvec2.n_elem; ++j)
        tmpvec2[j] += tmpval2;

      gr.row(i).tail(pZ) = (Z_i % tmpvec2).t();
    }
  }
  return gr;
}

