#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat weibull_hess(const arma::vec& theta, const arma::mat& X, Nullable<arma::mat> Z,
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

  vec log_tt = log(tt);
  vec log_tt0 = arma::vec(n, fill::zeros);
  for (int i=0; i<n; ++i)
    log_tt0[i] = tt0[i]>0 ? log(tt0[i]) : 0;

  // Hessian computation
  mat H(theta.n_elem, theta.n_elem, fill::zeros);

  // Hessian for beta terms
  mat Hbb = zeros<mat>(pX, pX);
  for (int i = 0; i < n; ++i) {
    double log_lambda_i = log_lambda[i];
    double lambda_i = lambda[i];
    double log_tt_i = log(tt[i]);
    double log_tt0_i = log(tt0[i]);
    double tt0_i = tt0[i];
    double w_i = w[i];
    double p_i = p[i];
    vec X_i = X.row(i).t();

    double termt0 = tt0_i>0 ? exp(p_i*log_tt0_i + log_lambda_i) : 0;
    double termt = exp(p_i*log_tt_i + log_lambda_i);
    double tmpval_bb = w_i * (termt0 - termt);

    vec tmpvec_bb(X_i.n_elem, fill::zeros);
    for(int j=0; j<tmpvec_bb.n_elem; ++j)
      tmpvec_bb[j] += tmpval_bb;

    mat tmp = (X_i % tmpvec_bb) * X_i.t();
    Hbb += tmp;
  }
  H.submat(0, 0, pX - 1, pX - 1) = Hbb;

  if (pZ > 0) {
    // Hessian for beta and gamma cross terms, and gamma^2 terms
    mat Hpb = zeros<mat>(pX, pZ);
    mat Hpp = zeros<mat>(pZ, pZ);
    for (int i = 0; i < n; ++i) {
      double log_lambda_i = log_lambda[i];
      double lambda_i = lambda[i];
      double log_tt_i = log_tt[i];
      double log_tt0_i = log_tt0[i];
      double tt0_i = tt0[i];
      double w_i = w[i];
      double p_i = p[i];
      double d_i = d[i];
      vec X_i = X.row(i).t();
      vec Z_i = as<mat>(Z).row(i).t();

    double termt0 = tt0_i>0 ? exp(p_i*log_tt0_i + log_lambda_i) : 0;
    double termt = exp(p_i*log_tt_i + log_lambda_i);
      double tmpval_pb = w_i * p_i * (termt0*log_tt0_i - termt*log_tt_i);
      vec tmpvec_pb(X_i.n_elem, fill::zeros);
      for(int j=0; j<tmpvec_pb.n_elem; ++j)
        tmpvec_pb[j] += tmpval_pb;

      mat tmp_pb = (X_i % tmpvec_pb) * Z_i.t() ;
      Hpb += tmp_pb;

      double tmpval_pp = w_i*p_i*(d_i*log_tt_i +
                                  (termt0*log_tt0_i*(1+p_i*log_tt0_i) -
                                  termt*log_tt_i*(1+p_i*log_tt_i)));

      vec tmpvec_pp(Z_i.n_elem, fill::zeros);
      for(int j=0; j<tmpvec_pp.n_elem; ++j)
        tmpvec_pp[j] += tmpval_pp;

      mat tmp_pp = (Z_i % tmpvec_pp) * Z_i.t();
      Hpp += tmp_pp;
    }
    H.submat(0, pX, pX-1, theta.n_elem - 1) = Hpb;
    H.submat(pX, 0, theta.n_elem - 1, pX - 1) = Hpb.t();
    H.submat(pX, pX, theta.n_elem - 1, theta.n_elem - 1) = Hpp;
  }
  return H;
}

