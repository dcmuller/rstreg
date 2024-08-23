#' @title Regression for a Parametric Survival Model
#' @description
#' Fit a parametric survival regression model. Currently this function only implements
#' the Weibull proportional hazards model with parameters \eqn{\lambda} and \eqn{\text{p}}.
#' This parameterisation is the same as used in Stata's \code{streg} command, but
#' differs from the parameterisation in [survival::survreg()] which uses the location-scale
#' framework fits models with the accelerated failure time metric. Unlike [survival::survreg()],
#' this function also allows delayed entry/counting process survival data. In future
#' this function will offer further parametric survival functions in both the proportional
#' hazards and accelerated failure time metrics.
#'
#'
#' @param formula a formula expression. The response should be a survival object as returned by the
#' [survival::Surv()] function.
#' @param data a data frame in which to interpret the variables named in the \code{formula},
#' \code{weights}, or \code{subset} arguments.
#' @param weights optional vector of observation weights
#' @param subset subset of the observations to be used in the fit
#' @param na.action a missing-data filter function, applied to the model.frame, after any subset
#' argument has been used.
#' @param dist Parametric survival distribution. Currently only \code{"weibull"} is available.
#' @param pfixed optional real value to fix the Weibull parameter p. e.g. set \code{pfixed=1} to
#' fit an exponential survival model.
#' @param init optional vector of initial values for the parameters.
#' @param max.method optimisation method passed to [maxLik::maxLik()]
#' @param control list of control values passed to [maxLik::maxLik()]
#' @param model if TRUE returns the model frame
#' @param x if TRUE returns the x matrix
#' @param z if TRUE returns the z matrix (predictor matrix for the Weibull parameter p)
#' @param y if TRUE returns the response
#' @param robust Use robust sandwich variance-covariance instead of the asymptotic
#' variance-covariance. If there is a cluster argument this defaults to TRUE, and the returned
#' robust variance-covariance matrix will be adjusted for the number of clusters rather than
#' the number of observations.
#' @param cluster optional variable that identifies groups of observations that is used in the
#' calculation of the robust variance-covariance
#' @param ... other arguments passed to [maxLik::maxLik()]
#'
#' @details
#' Currently fits proportional hazards model of the form \eqn{h(\mathbf{t})=h_0(\mathbf{t})\exp(g(\mathbf{X}))}.
#' The Weibull proportional hazards model has Survivor function \eqn{\exp(-\lambda t^p)}, and
#' is parameterised by \eqn{\lambda = \exp{\mathbf{X\beta}}}. When \code{formula} contains a
#' [survival::strata()] term, the Weibull parameter \eqn{p} is parameterised by
#' \eqn{p=\exp(\mathbf{Z\gamma})}, thus fitting separate baseline hazards for each level of
#' the \code{strata} term.
#'
#'
#' @return an object of class \code{streg} is returned.
#'
#' @author David C Muller
#'
#' @import survival
#'
#' @export
#'
streg <- function(formula, data, weights, subset, na.action, dist = "weibull", pfixed=NULL,
                  init = NULL, max.method="NR", control=NULL, model=TRUE,
                  x=TRUE, z=TRUE, y = TRUE, robust = FALSE, cluster,
                  ...)
{
  Call <- match.call()
  if (missing(formula))
    stop("a formula argument is required")
  ss <- c("cluster", "offset")
  Terms <- if (missing(data))
    terms(formula, specials = ss)
  else terms(formula, specials = ss, data = data)
  tcl <- attr(Terms, "specials")$cluster
  if (length(tcl) > 1)
    stop("a formula cannot have multiple cluster terms")
  if (length(tcl) > 0) {
    factors <- attr(Terms, "factors")
    if (any(factors[tcl, ] > 1))
      stop("cluster() cannot be in an interaction")
    if (attr(Terms, "response") == 0)
      stop("formula must have a Surv response")
    if (is.null(Call$cluster))
      Call$cluster <- attr(Terms, "variables")[[1 + tcl]][[2]]
    else warning("cluster appears both in a formula and as an argument, formula term ignored")
    Terms <- survival:::drop.special(Terms, tcl)
    formula <- Call$formula <- formula(Terms)
  }
  indx <- match(c("formula", "data", "weights", "subset",
                  "na.action", "cluster"), names(Call), nomatch = 0)
  if (indx[1] == 0)
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  special <- c("strata")
  temp$formula <- if (missing(data))
    terms(formula, special)
  else terms(formula, special, data = data)
  m <- eval(temp, parent.frame())
  Terms <- attr(m, "terms")
  weights <- model.extract(m, "weights")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")
  type <- attr(Y, "type")
  cluster <- model.extract(m, "cluster")
  if (length(cluster)) {
    if (missing(robust))
      robust <- TRUE
    cluster <- as.numeric(as.factor(cluster))
  }
  else if (robust)
    cluster <- 1:nrow(Y)
  strats <- attr(Terms, "specials")$strata
  keepz <- NULL
  if (length(strats)) {
    temp <- untangle.specials(Terms, "strata", 1)
    keepz <- temp$terms
    if (length(temp$vars) == 1)
      strata.keep <- m[[temp$vars]]
    else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
    strata <- as.numeric(strata.keep)
    nstrata <- max(strata)
  }
  else {
    nstrata <- 1
    strata <- 0
  }

  X <- model.matrix(Terms, m)
  xlevels <- .getXlevels(Terms, m)
  contr.save <- attr(X, "contrasts")
  if (!all(is.finite(X)))
    stop("data contains an infinite predictor")
  n <- nrow(X)
  nvar <- ncol(X)
  if (nstrata>1) {
    zTerms <- Terms[keepz]
    attr(zTerms, "intercept") <- attr(Terms, "intercept")
  }
  else {
    zTerms <- Terms[0]
    attr(zTerms, "intercept") <- attr(Terms, "intercept")
  }
  if (nstrata>1 && !is.null(pfixed))
    stop("cannot provide a fixed Weibull paramater pfixed in conjunction with strata()")
  if (is.null(pfixed)) {
    Z <- model.matrix(zTerms, m)
    zlevels <- .getXlevels(zTerms, m)
    contr.save.Z <- attr(Z, "contrasts")
    if (!all(is.finite(Z)))
      stop("data contains an infinite predictor")
    nvarZ <- ncol(Z)
  }
  else {
    Z <- NULL
    nvarZ <- NULL
    zlevels <- NULL
    contr.save.Z <- NULL
  }

  offset <- model.offset(m)
  if (length(offset) == 0 || all(offset == 0))
    offset <- rep(0, n)
  if (is.character(dist)) {
    dist <- match.arg(dist, c("exponential", "weibull", "gompertz", "lognormal", "loglogistic"))
  }
  logcorrect <- 0
  exactsurv <- Y[, ncol(Y)] == 1
  if (any(exactsurv)) {
    if (is.null(weights))
      logcorrect <- sum(log(Y[exactsurv, (ncol(Y)-1)]))
    else logcorrect <- sum(weights[exactsurv] * log(Y[exactsurv, (ncol(Y)-1)]))
  }
  Ysave <- Y

  fit <- streg.fit(X, Z, Y, weights, offset, init = init, pfixed,
                     max.method=max.method, control=control, dist = dist,
                     ...)
  if (is.character(fit))
    fit <- list(fail = fit)
  else {
    fit$maximum <- fit$maximum + logcorrect
    nvar <- ncol(X)
    fit$coefficients <- fit$estimate[1:nvar]
    fit$df <- length(fit$estimate)
    fit$df.residual <- n - fit$df
    fit$terms <- Terms
    fit$contrasts <- contr.save
    if (length(xlevels))
      fit$xlevels <- xlevels
    fit$means <- apply(X, 2, mean)
    fit$zterms <- zTerms
    if (length(contr.save.Z))
      fit$contrasts.Z <- contr.save.Z
    if (length(zlevels))
      fit$zlevels <- zlevels
    if (!is.null(weights))
      fit$weights <- weights
    if (!is.null(pfixed))
      fit$pfixed <- pfixed
    fit$call <- Call
    fit$dist <- dist
    fit$metric <- "PH"
    if (model)
      fit$model <- m
    if (x)
      fit$x <- X
    if (z && !is.null(Z))
      colnames(Z) <- paste0("(log_p)", colnames(Z))
      fit$z <- Z
    if (y)
      fit$y <- Ysave
    ## return variance matrices
    fit$var <- solve(-fit$hessian)
    if (robust) {
      fit$naive.var <- fit$var
      if (!model)
        fit$model <- m
      if (length(cluster))
        fit$var <- fit$naive.var%*%((max(cluster)/(max(cluster)-1))*crossprod(rowsum(fit$gradientObs, cluster)))%*%fit$naive.var
      else fit$var <- fit$naive.var%*%((n/(n-1))*crossprod(rowsum(fit$gradientObs)))%*%fit$naive.var
      if (!model)
        fit$model <- NULL
    }
    singular <- (diag(fit$var) == 0)[1:length(fit$coefficients)]
    if (any(singular))
      fit$coefficients[singular] <- NA
    na.action <- attr(m, "na.action")
    if (length(na.action))
      fit$na.action <- na.action
    class(fit) <- c("streg", "survreg", class(fit))
  }
  fit
}
