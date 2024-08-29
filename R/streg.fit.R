#' @title Fit a Parametric Survival Model
#' @description Raw interface for fitting parametric survival regression models. See
#' [rstreg::streg()] for the formula and data based interface.
#' @param x \code{model.matrix} for the regression parameters.
#' @param z \code{model.matrix} for the auxiliary parameters
#' @param y \code{Survival outcome of class [survival::Surv()]}
#' @param weights observation weights
#' @param offset offset for the linear predictor of the regression equation
#' @param init vector of initial values
#' @param pfixed value to fix auxiliary parameter p at.
#' @param control control list passed to [roptim]
#' @param dist survival distribution
#' @param max.method maximisation method passed to [roptim]
#' @param ... for future methods
#'
#' @export
streg.fit <- function(x, z, y, weights, offset, init, pfixed, control, dist, max.method, ...) {
  if (!is.matrix(x))
    stop("Invalid x matrix ")
  if (is.null(pfixed) && !is.matrix(z))
    stop("Invalid z matrix")
  n <- nrow(x)
  nvar <- ncol(x)
  nz <- ncol(z)
  ny <- ncol(y)
  if (is.null(offset))
    offset <- rep(0, n)
  if (missing(weights) || is.null(weights))
    weights <- rep(1, n)
  else if (any(weights <= 0))
    stop("Invalid weights, must be >0")
  if (!is.null(pfixed))
    pfixed <- rep(pfixed, n)
  nvar2 <- nvar + ifelse(!is.null(pfixed), 0L, nz)
  if (is.numeric(init)) {
    if (length(init) == nvar && (nvar2 > nvar)) {
      init <- c(init, rep(0,nz))
    }
    if (length(init) != nvar2)
      stop("Wrong length for initial parameters")
  }
  else init <- rep(0, nvar2)
  if (dist=="weibull" && !is.null(z)) {
    colnames(z) <- paste0("(log_p)", colnames(z))
  }
  names(init) <- c(colnames(x),colnames(z))
  if (ny==3) {
    status <- y[, "status"]
    exit <- y[, "stop"]
    enter <- y[, "start"]
  }
  else {
    status <- y[, "status"]
    exit <- y[, "time"]
    enter <- rep(0, length(exit))
  }
#  fit <- maxLik(weibull_ll, grad=weibull_gr, hess=weibull_hess, start=init,
#                X=x, Z=z, tt=exit, tt0=enter, d=status, pfixed=pfixed, w=weights, offset=offset, method=max.method, control=control)
  fit <- weibull(X=x, Z=z, tt=exit, tt0=enter, d=status,
                 pfixed=pfixed, w=weights, offset=offset,
                 method=max.method, control=control, print=0, init_theta=init)
  fit$par <- as.vector(fit$par)
  names(fit$par) <- names(init)
  fit$hessian_numerical <- fit$hessian
  fit$linear.predictors <- c(x %*% fit$par[colnames(x)] + offset)
  fit$gradientObs <- weibull_gr(fit$par,x,z,enter,exit,status,pfixed,weights,offset)
  dimnames(fit$gradientObs) <- list(rownames(x), names(init))
  fit$hessian <- weibull_hess(fit$par,x,z,enter,exit,status,pfixed,weights,offset)
  dimnames(fit$hessian) <- list(names(init), names(init))
  return(fit)
}
