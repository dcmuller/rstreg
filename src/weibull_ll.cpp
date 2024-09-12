#include <RcppArmadillo.h>
#include "streg.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
arma::vec weibull_ll(const arma::vec& theta, const arma::mat& X, Nullable<arma::mat> Z,
                          const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
                          const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec & offset) {

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
  vec log_lambda = X * beta + offset;
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

  // Compute log-likelihood
  vec ll = d % (log_lambda + log_p + (p - 1) % log(tt)) - lambda % exp(p % log(tt));
  for (int i = 0; i < n; ++i) {
    if (tt0[i] > 0) {
      ll[i] += lambda[i] * exp(p[i] * log(tt0[i]));
    }
    ll[i] *= w[i];
  }
  return ll;
}

