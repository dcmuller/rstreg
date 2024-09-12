#include <cmath>
#include <cstddef>
#include <RcppArmadillo.h>
#include "streg.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
List weibull_max_nr(const int maxiter, const arma::mat& X, const Nullable<arma::mat> Z,
        const arma::vec& tt0, const arma::vec & tt, const arma::vec& d, const Nullable<arma::vec> pfixed,
        const arma::vec& wt, const arma::vec & offset, const arma::vec& theta,
              const double eps, const double tol) {
//              XPtr<LikelihoodFunc> likelihood_func,
//              XPtr<GradFunc> grad_func,
//              XPtr<HessFunc> hess_func) {

  const int n = d.n_elem;       // Number of observations
  const int pX = X.n_cols;      // Number of parameters for X
  const int pZ = Z.isNotNull() ? as<mat>(Z).n_cols : 0; // Number of parameters for Z (if applicable)

  // Check that the length of theta matches the sum of columns of X and Z
  if (theta.n_elem != (pX + pZ)) {
    stop("Length of theta must be equal to the sum of the number of columns in X and Z.");
  }
  const int nvar2 = theta.n_elem;

  arma::vec newtheta(theta.n_elem);
  arma::vec par(theta.n_elem); // working parameter vector
  arma::mat imat(nvar2, nvar2);
  arma::mat JJ(nvar2, nvar2);
  arma::mat grobs(n, nvar2);
  arma::vec u(nvar2);
  arma::vec usave(nvar2);
  arma::vec loglik(1);
  bool flag;
  int halving = 0;
  int finaliter = 0;
  arma::vec step(nvar2);
  double newlk = 0.0;

  newtheta = theta;
  par = theta;
  // Call the likelihood function for the initial values
  weibull_inplace(loglik[0], grobs, u, imat,
                  theta, X, Z, tt0, tt, d,
                  pfixed, wt, offset);
  imat = -imat;
  flag = imat.is_sympd(tol);
  if (!flag) {
    JJ = grobs.t() * grobs;
    step = solve(JJ, u);
  }
  else
    step = solve(imat, u);

  // Update estimates using Newton-Raphson step
  newtheta = theta + step;

  weibull_inplace(newlk, grobs, u, imat,
                  newtheta, X, Z, tt0, tt, d,
                  pfixed, wt, offset);
  imat = -imat;

  // main loop
  for (int iter = 1; iter <= maxiter; iter++) {
    if (!std::isnan(newlk) && !std::isinf(newlk) &&
        fabs(1 - (loglik[0] / newlk)) <= eps && halving == 0) {
      loglik[0] = newlk;
      flag = imat.is_sympd(tol);
      finaliter = iter-1;
      par = newtheta;
//      return List::create(Named("coef") = par,
//                          Named("iter") = finaliter,
//                          Named("hessian") = -imat,
//                          Named("loglik") = loglik,
//                          Named("flag") = flag,
//                          Named("u") = u);
      break;
    }

    if (std::isnan(newlk) || std::isinf(newlk) || newlk < loglik[0]) {
      // Step halving
      for (int j = 0; j < 10 && newlk < loglik[0]; j++) {
        halving++;
        newtheta = (newtheta + par) / 2;
        if (halving==1) {  /* only the first time */
          for (int i=pX; i<theta.n_elem; i++) {
            if ((par[i]-newtheta[i])> 1.1)
              newtheta[i] = par[i] - 1.1;
          }
        }
        weibull_inplace(newlk, grobs, u, imat,
                        newtheta, X, Z, tt0, tt, d,
                        pfixed, wt, offset);
        imat = -imat;
      }
    } else {
      // Standard Newton-Raphson step
      halving = 0;
      flag = imat.is_sympd(tol);
      if (!flag) {
        JJ = grobs.t() * grobs;
        step = solve(JJ, u);
      }
      else
        step = solve(imat, u);
      par = newtheta;
      newtheta += step;
    }
    loglik[0] = newlk;
    weibull_inplace(newlk, grobs, u, imat,
                    newtheta, X, Z, tt0, tt, d,
                    pfixed, wt, offset);
    imat = -imat;
  }

  // If max iterations reached without convergence
  if (finaliter==0) {
    finaliter = maxiter;
    flag = false;
  }
  return List::create(Named("par") = par,
                      Named("value") = loglik,
                      Named("iter") = finaliter,
                      Named("hessian") = -imat,
                      Named("u") = u,
                      Named("convergence") = flag ? 0 : 1);
}
