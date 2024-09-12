#ifndef STREG_H
#define STREG_H
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

void weibull_inplace(double &llval, arma::mat &grObsval, arma::vec &grval, arma::mat &hessval,
                     const arma::vec& theta, const arma::mat& X, const Nullable<arma::mat>& Z,
                     const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
                     const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec& offset);
arma::vec weibull_ll(const arma::vec& theta, const arma::mat& X, Nullable<arma::mat>& Z,
                     const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
                     const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec& offset);
arma::mat weibull_gr(const arma::vec& theta, const arma::mat& X, Nullable<arma::mat>& Z,
                     const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
                     const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec& offset);
arma::mat weibull_hess(const arma::vec& theta, const arma::mat& X, Nullable<arma::mat>& Z,
                       const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
                       const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec& offset);


//void weibull_inplaceT(double &llval, arma::mat &grObsval, arma::vec &grval, arma::mat &hessval,
//                     const arma::vec& theta, const arma::mat& X, const Nullable<arma::mat>& Z,
//                     const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
//                     const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec& offset);
//typedef void (*funcPtr)(double &llval, arma::mat &grObsval, arma::vec &grval, arma::mat &hessval,
//                        const arma::vec& theta, const arma::mat& X, const Nullable<arma::mat>& Z,
//                        const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
//                        const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec& offset);
//void exec_weibull_inplaceT(SEXP func_, double &llval, arma::mat &grObsval, arma::vec &grval, arma::mat &hessval,
//                           const arma::vec& theta, const arma::mat& X, const Nullable<arma::mat> Z,
//                           const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
//                           const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec& offset);


//template <typename T>
//arma::vec ll(T func, const arma::vec& theta, const arma::mat& X, Nullable<arma::mat> Z,
//             const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
//             const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec & offset) {
//  return (func(theta, X, Z, tt0, tt, d, pfixed, w, offset));
//};
//
//template <typename T>
//arma::mat gr(T func, const arma::vec& theta, const arma::mat& X, Nullable<arma::mat> Z,
//                 const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
//                 const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec & offset) {
//  return (func(theta, X, Z, tt0, tt, d, pfixed, w, offset));
//};
//
//template <typename T>
//arma::mat hess(T func, const arma::vec& theta, const arma::mat& X, Nullable<arma::mat> Z,
//                 const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
//                 const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec & offset) {
//  return (func(theta, X, Z, tt0, tt, d, pfixed, w, offset));
//};
//
//// Define the type for the likelihood etc function pointers
//typedef SEXP (*LikelihoodFunc)(const arma::vec& theta, const arma::mat& X, Nullable<arma::mat> Z,
//              const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
//              const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec & offset);
//typedef SEXP (*GradFunc)(const arma::vec& theta, const arma::mat& X, Nullable<arma::mat> Z,
//              const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
//              const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec & offset);
//typedef SEXP (*HessFunc)(const arma::vec& theta, const arma::mat& X, Nullable<arma::mat> Z,
//              const arma::vec& tt0, const arma::vec& tt, const arma::vec& d,
//              const Nullable<arma::vec>& pfixed, const arma::vec& w, const arma::vec & offset);
//
#endif

