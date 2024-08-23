#' @export
streg.fit <- function(x, z, y, weights, offset, init, pfixed, control, dist, max.method, ...) {
  if (!is.matrix(x))
    stop("Invalid X matrix ")
  if (is.null(pfixed) && !is.matrix(z))
    stop("Invalid Z matrix")
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
  fit <- maxLik(streg.fit.ll.weibull, start=init,
                X=x, Z=z, tt=exit, tt0=enter, d=status, pfixed=pfixed, w=weights, method=max.method, control=control)
  fit$linear.predictors <- c(x %*% fit$estimate[colnames(x)] + offset)
  return(fit)
}
