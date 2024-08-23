#' @title Predicted Values for a 'streg' Object
#'
#' @description
#' Predicted values for a \code{streg} object
#'
#' @param object result of a model fit using the \code{streg} function
#' @param newdata data for prediction. If absent predictions are for observations
#' used in the original fit.
#' @param type the type of predicted value. \code{"lp"}, \code{"linear"}, and
#' \code{"xb"} provide the linear predictor, and \code{"hr"} provides the antilog of
#' the linear predictor. \code{"p"} and \code{"log.p"} provide the Weibull shape
#' parameter and it's logarithm, respectively (constant for models that do not include a
#' [survival::strata()] term). \code{"log.hazard"} and \code{"hazard"} provide the
#' hazard on the log and natural scale respectively, as do \code{"log.cum.hazard"} and
#' \code{"cum.hazard"} for the cumulative hazard. \code{"survival"} provides the survival
#' probability. \code{"uquantile"} provides the predicted quantile on the linear
#' predictor scale, and \code{"quantile"} on the data scale. \code{terms} provides the
#' matrix of centered terms for the linear predictor.
#' @param nocons if TRUE removes the intercept from the linear predictor. Only valid for
#' types \code{"lp"}, \code{"linear"}, \code{"xb"}, and \code{"hr"}.
#' @param nooffset if TRUE removes the offset from the linear predictor. Only valid for
#' types \code{"lp"}, \code{"linear"}, \code{"xb"}, and \code{"hr"}.
#' @param se.fit if TRUE, includes the standard errors of the prediction in the result.
#' The delta method is used for predictions that are non-linear transformations of the
#' linear predictors.
#' @param ci.fit if TRUE, includes a matrix with lower and upper limits for the confidence
#' interval. For type \code{"hr"}, intervals are calculated on the linear predictor scale
#' and transformed. Similarly for types \code{"hazard"} and \code{"cum.hazard"}, intervals
#' are calculated on the log scale and then transformed. Intervals for type
#' \code{"survival"} are calculated for the log cumulative hazard and transformed.
#' @param ci.level level for calculated confidence intervals
#' @param terms subset of terms. The default for residual type "terms" is a matrix
#' with one column for every term, excluding the intercept, in the model.
#' @param perc vector of percentiles. This is used only for quantile predictions.
#' @param na.action applies only when the newdata argument is present, and defines the
#' missing value action for the new data. The default is to include all observations.
#' @param ... for future methods
#'
#' @return A vector or matrix of predicted values. If \code{se.fit} or \code{ci.fit}
#' are TRUE, the return type is a list containing, where relevant, \code{fit},
#' \code{se.fit}, and \code{ci.fit} objects.
#'
#' @seealso [rstreg::streg()]
#' @author David C Muller
#'
#' @import survival
#' @import maxLik
#'
#' @export

predict.streg <- function (object, newdata, type = c("lp", "linear", "xb", "hr",
                                                     "log.p", "p",
                                                     "log.hazard", "hazard",
                                                     "log.cum.hazard", "cum.hazard",
                                                     "survival", "quantile", "uquantile", "terms"),
                           nocons=FALSE, nooffset=FALSE,
                           se.fit = FALSE, ci.fit=FALSE, ci.level=0.95,
                           terms = NULL, perc = c(0.1, 0.5, 0.9), na.action = na.pass, ...)
{
  if (ci.fit) {
    if (type=="quantile" || type=="uquantile")
      stop("ci.fit=TRUE is not valid for quantile predictions")
    se.fit <- TRUE
    ci.crit <- c(1,-1)*qnorm((1- ci.level)/2)
  }
  type <- match.arg(type)
  if (type == "linear" || type=="xb")
    type <- "lp"
  if (!(type == "lp" || type=="hr") && nocons)
    stop(paste0("nocons=TRUE is not a valid option for type = '", type,"'"))
  n <- nrow(object$gradientObs)
  Terms <- object$terms
  zTerms <- object$zterms
  if (!inherits(Terms, "terms"))
    stop("invalid terms component of  object")
  if (!inherits(zTerms, "terms") && !is.null(object$z))
    stop("invalid Z terms component of  object")
  strata <- attr(Terms, "specials")$strata
  Terms <- delete.response(Terms)
  zTerms <- delete.response(zTerms)
  coef <- object$coefficients
  intercept <- attr(Terms, "intercept")
  nvar <- length(object$coefficients)
  fixedscale <- (nvar == ncol(object$var))
  vvb <- object$var[1:nvar, 1:nvar]
  vvfull <- object$var
  coeffull <- object$estimate
  if (!fixedscale) {
    coefz <- coeffull[(length(coef)+1):length(coeffull)]
    vvz <- vvfull[(length(coef)+1):length(coeffull),(length(coef)+1):length(coeffull)]
  }
  else {
    coefz <- log(object$pfixed)
    vvz <- as.matrix(0)
  }
  if (missing(newdata) && (type != "lp" || se.fit))
    need.x <- need.z <- TRUE
  else need.x <- need.z <- FALSE
  if (type=="log.hazard" || type=="hazard" ||
      type=="log.cum.hazard" || type=="cum.hazard" || type=="survival")
    need.y <- TRUE
  else need.y <- FALSE
  if (!missing(newdata)) {
    newframeX <- stats::model.frame(Terms, data = newdata,
                                   na.action = na.action, xlev = object$xlevels)
    na.action.used <- attr(newframeX, "na.action")
    newframeZ <- stats::model.frame(zTerms, data=newdata,
                                    na.action = na.action, xlev = object$zlevels)
    if (need.y) {
      if (!any(lapply(newdata, class)=="Surv"))
        stop(paste0("prediction type '", type, "' requires one variable in newdata to be a valid Surv object"))
      if (sum(as.numeric(lapply(newdata, class)=="Surv"))>1)
        stop("multiple survival objects found in newdata")
      newY <- newdata[,(lapply(newdata, class)=="Surv")]
    }
  }
  else na.action.used <- object$na.action
  if (length(strata) && !fixedscale) {
    mf <- stats::model.frame(object)
    temp <- untangle.specials(Terms, "strata", 1)
    if (length(temp$vars) == 1)
      strata.keep <- mf[[temp$vars]]
    else strata.keep <- strata(mf[, temp$vars], shortlabel = TRUE)
    strata <- as.numeric(strata.keep)
    nstrata <- max(strata)
    if (missing(newdata)) {
      if (need.x) {
        x <- object[["x"]]
        if (is.null(x))
          x <- model.matrix(object, mf)
      }
      if (need.z) {
        z <- object[["z"]]
        if (is.null(z))
          z <- model.matrix(zTerms, mf)
      }
      if (need.y) {
        y <- object[["y"]]
        if (is.null(y))
          stop("this call to predict needs an streg fit with option y=TRUE")
      }
      offset <- model.offset(mf)
      if (is.null(offset))
        offset <- 0
    }
    else if (!missing(newdata)) {
      if (length(temp$vars) == 1)
        newstrat <- newframeX[[temp$vars]]
      else newstrat <- strata(newframeX[, temp$vars], shortlabel = TRUE)
      strata <- match(newstrat, levels(strata.keep))
      x <- model.matrix(object, newframeX)
      z <- model.matrix(zTerms, newframeZ)
      offset <- model.offset(newframeX)
      if (is.null(offset))
        offset <- 0
      if (need.y)
        y <- newY
    }
  }
  else {
    nstrata <- 1
    if (missing(newdata)) {
      mf <- stats::model.frame(object)
      strata <- rep(1L, n)
      if (need.x)
        x <- model.matrix(object)
      if (need.z)
        z <- model.matrix(zTerms, mf)
      if (need.y)
        y <- object$y
      offset <- model.offset(mf)
      if (is.null(offset))
        offset <- 0
    }
    else {
      x <- model.matrix(object, newframeX)
      z <- model.matrix(zTerms, newframeZ)
      if (need.y)
        y <- newY
      strata <- rep(1L, nrow(x))
      offset <- model.offset(newframeX)
      if (is.null(offset))
        offset <- 0
    }
  }
  if (type == "terms" && intercept)
    x <- sweep(x, 2, object$means)
  if (is.character(object$dist))
    dd <- streg.distributions[[object$dist]]
  else dd <- object$dist

  if (!is.null(dd$dist))
    dd <- streg.distributions[[dd$dist]]
  if ((type == "lp" || type == "hr")) {
    if (nocons)
      x[, "(Intercept)"] <- 0
    pred <- drop(x %*% coef) + as.numeric(nooffset)*offset
    if (se.fit)
      se <- sqrt(diag(x %*% vvb %*% t(x)))
    if (ci.fit) {
      ci <- pred + outer(se, ci.crit)
      colnames(ci) <- c(paste0(c("l","u"),100*ci.level))
    }
    if (type == "hr") {
      pred <- exp(pred)
      if (se.fit)
        se <- se*pred
    if (ci.fit)
      ci <- exp(ci)
    }
  }
  else if (type == "log.p" || type == "p") {
    pred <- drop(z %*% coefz)
    if (se.fit)
      se <- sqrt(diag(z %*% vvz %*% t(z)))
    if (ci.fit) {
      ci <- pred + outer(se, ci.crit)
      colnames(ci) <- c(paste0(c("l","u"),100*ci.level))
    }
    if (type == "p") {
      pred <- exp(pred)
      if (se.fit)
        se <- se*pred
      if (ci.fit)
        ci <- exp(ci)
    }
  }
  else if (type=="log.hazard" || type=="hazard") {
    if (ncol(y)==2)
      tt <- y[, 1]
    else tt <- y[, 2]
    xb <- drop(x %*% coef)
    lp <- drop(z %*% coefz)
    pred <- dd$lh(xb,lp,tt)
    if (se.fit) {
      gr <- dd$dlh(x,z,lp,tt)
      se <- sqrt(diag(gr %*% vvfull %*% t(gr)))
    }
    if (ci.fit) {
      ci <- pred + outer(se, ci.crit)
      colnames(ci) <- c(paste0(c("l","u"),100*ci.level))
    }
    if (type=="hazard") {
      pred <- exp(pred)
      if (se.fit)
        se <- se*pred
      if (ci.fit) {
        ci <- exp(ci)
        colnames(ci) <- c(paste0(c("l","u"),100*ci.level))
      }
    }
  }
  else if (type=="log.cum.hazard" || type=="cum.hazard" || type=="survival") {
    if (ncol(y)==2) {
      tt <- y[, 1]
      tt0 <- rep(0, length(tt))
    }
    else {
      tt <- y[, 2]
      tt0 <- y[, 1]
    }
    xb <- drop(x %*% coef)
    lp <- drop(z %*% coefz)
    pred <- dd$lH(xb,lp,tt,tt0)
    if (se.fit) {
      gr <- dd$dlH(x,z,lp,tt,tt0)
      se <- sqrt(diag(gr %*% vvfull %*% t(gr)))
    }
    if (ci.fit) {
      ci <- pred + outer(se, ci.crit)
      colnames(ci) <- c(paste0(c("l","u"),100*ci.level))
    }
    if (type=="cum.hazard" || type=="survival") {
      pred <- exp(pred)
      if (se.fit)
        se <- se*pred
      if (ci.fit) {
        ci <- exp(ci)
        colnames(ci) <- c(paste0(c("l","u"),100*ci.level))
      }
    }
    if (type=="survival") {
      pred <- exp(-pred)
      if (se.fit)
        se <- se*pred
      if (ci.fit) {
        ci <- exp(-ci)
        colnames(ci) <- c(paste0(c("u","l"), 100*ci.level))
        ci <- ci[,c(2,1)]
      }
    }
  }
  else if (type == "quantile" || type == "uquantile") {
    xb <- drop(x %*% coef)
    lp <- drop(z %*% coefz)
    pred <- lapply(perc, dd$qf, xb=xb, lp=lp)
    if (se.fit) {
      gr <- lapply(perc, dd$dqf, x, z, xb, lp)
      se <- lapply(gr, function(x) sqrt(diag(x %*% vvfull %*% t(x))))
    }
    if (type=="uquantile") {
      pred <- lapply(pred,log)
      if(se.fit) {
        for (i in 1L:length(pred))
          se[[i]] <- se[[i]]/exp(pred[[i]])
      }
    }
    pred <- do.call(cbind, pred)
    se <- do.call(cbind, se)
  }
  else {
    asgn <- attrassign(x, Terms)
    hasintercept <- attr(Terms, "intercept") > 0
    if (hasintercept)
      asgn$"(Intercept)" <- NULL
    nterms <- length(asgn)
    pred <- matrix(ncol = nterms, nrow = NROW(x))
    dimnames(pred) <- list(rownames(x), names(asgn))
    if (se.fit) {
      se <- matrix(ncol = nterms, nrow = NROW(x))
      dimnames(se) <- list(rownames(x), names(asgn))
      R <- object$var
    }
    for (i in 1:nterms) {
      ii <- asgn[[i]]
      pred[, i] <- x[, ii, drop = FALSE] %*% (coef[ii])
      if (se.fit) {
        for (j in (1:NROW(x))) {
          xi <- x[j, ii, drop = FALSE] * (coef[ii])
          vci <- R[ii, ii]
          se[j, i] <- sqrt(sum(xi %*% vci %*% t(xi)))
        }
      }
    }
    if (!is.null(terms)) {
      pred <- pred[, terms, drop = FALSE]
      if (se.fit)
        se <- se[, terms, drop = FALSE]
    }
  }
  if (!is.null(na.action.used)) {
    pred <- naresid(na.action.used, pred)
    if (se.fit)
      se <- napredict(na.action.used, se)
    if (ci.fit)
      ci <- napredict(na.action.used, ci)
  }
  if (ci.fit)
    list(fit = pred, se.fit = se, ci.fit=ci)
  else if (se.fit)
    list(fit=pred, se.fit=se)
  else pred
}


streg.distributions <- vector(mode="list")
streg.distributions[["weibull"]] <- list(
  name="Weibull",
  lh = function(xb, lp, tt) {xb + exp(lp)*log(tt) - log(tt)},
  dlh = function(x,z,lp,tt) {
    cbind(x,
          z*exp(lp)*log(tt))
  },
  lH = function(xb,lp,tt,tt0) {
    xb + log(exp(exp(lp)*log(tt)) - ifelse(tt0>0, exp(exp(lp)*log(tt0)), 0))
  },
  dlH = function(x,z,lp,tt,tt0) {
    cbind(x,
          z * ((exp(lp) *
                  (log(tt)*exp(exp(lp)*log(tt)) - ifelse(tt0>0, log(tt0)*exp(exp(lp)*log(tt0)),0)))/
                 ((exp(exp(lp)*log(tt)) - ifelse(tt0>0, exp(exp(lp)*log(tt0)),0)))) )
  },
  qf = function(qq,xb,lp) {exp(1/exp(lp)*log(-log((1-qq))/exp(xb)))},
  dqf = function(qq,x,z,xb,lp) {
    cbind(
      -(x * (exp(log(-(log(1 - qq)/exp(xb)))/exp(lp))/exp(lp))),
      -(z * (exp(log(-(log(1 - qq)/exp(xb)))/exp(lp)) * (log(-(log(1 - qq)/exp(xb)))/exp(lp))))
    )
  }
)
