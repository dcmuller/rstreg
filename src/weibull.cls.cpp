#include <cmath>
#include <cstddef>

#include <algorithm>

#include <RcppArmadillo.h>

#include "roptim.h"

using namespace roptim;
using namespace Rcpp;
using namespace arma;

class Weibull : public Functor {
public:
  arma::mat X;
  Nullable<arma::mat> Z;
  arma::vec tt0;
  arma::vec tt;
  arma::vec d;
  arma::vec pfixed_val;
  arma::vec w;
  arma::vec offset;
  arma::vec init_theta;
  bool p_is_fixed;

  // constructor
  Weibull(const arma::mat& X, Nullable<arma::mat> Z,
          const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
          const Nullable<arma::vec>& pfixed,
          const arma::vec& w, const arma::vec & offset)
    : X(X), Z(Z), tt0(tt0), tt(tt), d(d), p_is_fixed(pfixed.isNotNull()),
      w(w), offset(offset) {
    if (p_is_fixed) {
      pfixed_val = as<arma::vec>(pfixed);
    } else {
      pfixed_val = arma::ones<arma::vec>(d.n_elem);
    }
  }

  double operator()(const arma::vec& theta) override {
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
  if (p_is_fixed) {
    log_p = log(pfixed_val);
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
  return sum(ll);
  }
  void Gradient(const arma::vec &theta, arma::vec &gr) override {
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
    if (p_is_fixed) {
      log_p = log(pfixed_val);
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

    // Gradient computation
    gr = arma::zeros<arma::vec>(theta.n_elem);

    // Gradient for beta terms
    for (int i = 0; i < n; ++i) {
      double lambda_i = lambda[i];
      double tt0_i = tt0[i];
      double log_tt_i = log_tt[i];
      double log_tt0_i = log_tt0[i];
      double w_i = w[i];
      double p_i = p[i];
      double d_i = d[i];
      vec X_i = X.row(i).t();
      double term1 = tt0_i > 0 ? exp(p_i*log_tt0_i) : 0;
      double tmpval_b = w_i * (d_i + lambda_i * (term1 - exp(p_i*log_tt_i)));
      vec tmpvec_b(X_i.n_elem, fill::zeros);
      for(int j=0; j<tmpvec_b.n_elem; ++j)
        tmpvec_b[j] += tmpval_b;

      gr.head(pX) += (X_i % tmpvec_b);

      // gradient for gamma terms
      if (pZ > 0) {
        vec Z_i = as<mat>(Z).row(i).t();
        double tmpval_p = w_i * (d_i*(1 + p_i*log_tt_i) +
                                 lambda_i*p_i*(exp(p_i*log_tt0_i)*log_tt0_i - exp(p_i*log_tt_i)*log_tt_i));
        vec tmpvec_p(Z_i.n_elem, fill::zeros);
        for(int j=0; j<tmpvec_p.n_elem; ++j)
          tmpvec_p[j] += tmpval_p;

        gr.tail(pZ) += (Z_i % tmpvec_p);
      }
    }
  }
  void Hessian(const arma::vec &theta, arma::mat &H) override {
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
    if (p_is_fixed) {
      log_p = log(pfixed_val);
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

    // Hessian for beta terms
    mat Hbb = zeros<mat>(pX, pX);
    for (int i = 0; i < n; ++i) {
      double lambda_i = lambda[i];
      double log_tt_i = log(tt[i]);
      double log_tt0_i = log(tt0[i]);
      double tt0_i = tt0[i];
      double w_i = w[i];
      double p_i = p[i];
      vec X_i = X.row(i).t();

      double term1 = tt0_i>0 ? exp(p_i*log_tt0_i) : 0;
      double tmpval_bb = w_i * lambda_i * (term1 - exp(p_i*log_tt_i));

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
        double lambda_i = lambda[i];
        double log_tt_i = log_tt[i];
        double log_tt0_i = log_tt0[i];
        double tt0_i = tt0[i];
        double w_i = w[i];
        double p_i = p[i];
        double d_i = d[i];
        vec X_i = X.row(i).t();
        vec Z_i = as<mat>(Z).row(i).t();

        double term1 = tt0_i > 0 ? exp(p_i*log_tt0_i)*log_tt0_i : 0;
        double tmpval_pb = w_i * lambda_i * p_i * (term1 - exp(p_i*log_tt_i)*log_tt_i);
        vec tmpvec_pb(X_i.n_elem, fill::zeros);
        for(int j=0; j<tmpvec_pb.n_elem; ++j)
          tmpvec_pb[j] += tmpval_pb;

        mat tmp_pb = (X_i % tmpvec_pb) * Z_i.t() ;
        Hpb += tmp_pb;

        double tmpval_pp = w_i*(d_i*p_i*log_tt_i + lambda_i *
                                (term1*p_i*log_tt0_i*(1+p_i*log_tt0_i) -
                                exp(p_i*log_tt_i)*p_i*log_tt_i*(1+p_i*log_tt_i)));

        vec tmpvec_pp(Z_i.n_elem, fill::zeros);
        for(int j=0; j<tmpvec_pp.n_elem; ++j)
          tmpvec_pp[j] += tmpval_pp;

        mat tmp_pp = (Z_i % tmpvec_pp) * Z_i.t();
        Hpp += tmp_pp;
      }
      H.submat(0, pX, pX-1, theta.n_elem - 1) = Hpb;
      H.submat(pX, 0, pX + pZ - 1, pX - 1) = Hpb.t();
      H.submat(pX, pX, theta.n_elem - 1, theta.n_elem - 1) = Hpp;
    }
  }
};

// [[Rcpp::export]]
Rcpp::List weibull(const arma::mat& X, Nullable<arma::mat> Z,
                  const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
                  const Nullable<arma::vec>& pfixed,
                  const arma::vec& w, const arma::vec & offset, bool print=TRUE,
                  const arma::vec init_theta=arma::vec(), const Rcpp::String method="BFGS",
                  const Nullable<Rcpp::List>& control=NULL) {
  Weibull ll(X, Z, tt0, tt, d, pfixed, w, offset);
  Roptim<Weibull> opt(method);

  // default control options
  opt.control.fnscale = -1;
  opt.control.trace = 0;
  opt.control.maxit = 500;
  opt.control.reltol = 1e-12;
  opt.control.factr = 1e12;

  if(control.isNotNull()) {
    Rcpp::List cont = as<Rcpp::List>(control);
    if (cont.containsElementNamed("trace"))
      opt.control.trace = cont.offset("trace");
    if (cont.containsElementNamed("fnscale"))
      opt.control.fnscale = cont.offset("fnscale");
    if (cont.containsElementNamed("parscale"))
      opt.control.parscale = cont.offset("parscale");
    if (cont.containsElementNamed("ndeps"))
      opt.control.ndeps = cont.offset("ndeps");
    if (cont.containsElementNamed("maxit"))
      opt.control.maxit = cont.offset("maxit");
    if (cont.containsElementNamed("abstol"))
      opt.control.abstol = cont.offset("abstol");
    if (cont.containsElementNamed("reltol"))
      opt.control.reltol = cont.offset("reltol");
    if (cont.containsElementNamed("alpha"))
      opt.control.alpha = cont.offset("alpha");
    if (cont.containsElementNamed("beta"))
      opt.control.beta = cont.offset("beta");
    if (cont.containsElementNamed("gamma"))
      opt.control.gamma = cont.offset("gamma");
    if (cont.containsElementNamed("REPORT"))
      opt.control.REPORT = cont.offset("REPORT");
    if (cont.containsElementNamed("warn_1d_NelderMead"))
      opt.control.warn_1d_NelderMead = cont.offset("warn_1d_NelderMead");
    if (cont.containsElementNamed("lmm"))
      opt.control.lmm = cont.offset("lmm");
    if (cont.containsElementNamed("factr"))
      opt.control.factr= cont.offset("factr");
    if (cont.containsElementNamed("pgtol"))
      opt.control.pgtol = cont.offset("pgtol");
  }

  opt.control.trace = print ? 1 : 0;
  arma::vec init_params = init_theta;
  opt.minimize(ll, init_params);

  if (print) {
    Rcpp::Rcout << "-------------------------" << std::endl;
    opt.print();
  }

  return Rcpp::List::create(
    Rcpp::Named("par") = opt.par(),
    Rcpp::Named("value") = opt.value(),
    Rcpp::Named("fncount") = opt.fncount(),
    Rcpp::Named("grcount") = opt.grcount(),
    Rcpp::Named("convergence") = opt.convergence(),
    Rcpp::Named("message") = opt.message(),
    Rcpp::Named("hessian") = opt.hessian(),
    Rcpp::Named("control") = List::create(Named("fnscale") = opt.control.fnscale,
                                          Named("trace") = opt.control.trace,
                                          Named("maxit") = opt.control.maxit,
                                          Named("abstol") = opt.control.abstol,
                                          Named("reltol") = opt.control.reltol,
                                          Named("alpha") = opt.control.alpha,
                                          Named("beta") = opt.control.beta,
                                          Named("gamma") = opt.control.gamma,
                                          Named("REPORT") = opt.control.REPORT,
                                          Named("warn1dNelderMead") = opt.control.warn_1d_NelderMead,
                                          Named("lmm") = opt.control.lmm,
                                          Named("lmm") = opt.control.lmm,
                                          Named("factr") = opt.control.factr,
                                          Named("pgtol") = opt.control.pgtol));
}
