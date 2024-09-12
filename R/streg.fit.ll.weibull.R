streg.fit.ll.weibull <- function(theta,X,Z,tt0=NULL,tt,d,pfixed,w=NULL,offset=offset,...) {
  if (is.null(tt0)) {
    tt0 <- rep(0, length(d))
  }
  if (is.null(w)) {
    w <- rep(1, length(d))
  }
  beta <- theta[colnames(X)]
  log_lambda <- drop(X%*%beta) + offset
  lambda <- exp(log_lambda)
  gamma <- theta[colnames(Z)]
  if (!is.null(pfixed))
    log_p <- log(pfixed)
  else
    log_p <- drop(Z%*%gamma)
  p <- exp(log_p)

  log_tt <- log(tt)
  log_tt0 <- ifelse(tt0>0,log(tt0),0)

  ll <- w*(d*(log_lambda + log_p + (p-1)*log(tt)) - lambda*exp(p*log(tt)) + ifelse(tt0>0, lambda*exp(p*log(tt0)), 0))

  ## matrix of observation-level gradient vectors (scores)
  gr <- matrix(NA, nrow=length(d), ncol=length(theta))
  colnames(gr) <- names(theta)
  gr[, colnames(X)] <- w*(X*(d+lambda*(ifelse(tt0>0,tt0^p,0) - tt^p)))
  if (!is.null(Z))
    gr[, colnames(Z)] <- w*(Z*(d*(1 + p*log(tt)) + lambda*p*(ifelse(tt0>0, tt0^p*log(tt0),0) - exp(log(tt)*p)*log(tt))))
  colnames(gr) <- names(theta)

  ## build hessian matrix of second derivatives
  H <- matrix(NA, nrow=length(theta), ncol=length(theta))
  colnames(H) <- names(theta)
  rownames(H) <- names(theta)
  Hbb <- t(w*X*as.vector(lambda) * as.vector(ifelse(tt0>0, exp(p*log(tt0)),0) - exp(p*log(tt)))) %*% X
  H[names(beta),names(beta)] <- Hbb
  if (!is.null(Z)) {
    Hpb <- t(w*X*as.vector(lambda)*as.vector(p)*
               as.vector(ifelse(tt0>0,exp(p*log(tt0)) * log(tt0),0) - exp(p*log(tt)) * log(tt))) %*% Z
    Hpp <- t(w * Z *
               as.vector(d*p*log_tt - lambda*exp(p*log_tt)*p*log_tt*(1+p*log_tt) +
                           lambda*exp(p*log_tt0)*p*log_tt0*(1+p*log_tt0))) %*% Z
    H[rownames(Hpb), colnames(Hpb)] <- Hpb
    H[colnames(Hpb), rownames(Hpb)] <- t(Hpb)
    H[rownames(Hpp), colnames(Hpp)] <- Hpp
  }

  attr(ll, "gradient") <- gr
  attr(ll, "Hessian") <- H
  return(ll)
}
