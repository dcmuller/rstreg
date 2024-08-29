#' @export
print.streg <- function (x, ...)
{
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  if (!is.null(x$fail)) {
    cat(" streg failed.", x$fail, "\n")
    return(invisible(x))
  }
  coef <- x$par
  if (any(nas <- is.na(coef))) {
    if (is.null(names(coef)))
      names(coef) <- paste("b", 1:length(coef), sep = "")
    cat("\nCoefficients: (", sum(nas), " not defined because of singularities)\n",
        sep = "")
  }
  else cat("\nCoefficients:\n")
  print(coef, ...)
  if (!is.null(x$pfixed))
    cat("\np fixed at", format(x$pfixed), "\n")

  nobs <- length(x$linear.predictors)
  cat("\nLog Likelihood (model) = ", format(round(x$value[1], 3)))
  cat("\n")
  omit <- x$na.action
  if (length(omit))
    cat("n=", nobs, " (", naprint(omit), ")\n", sep = "")
  else cat("n=", nobs, "\n")
  invisible(x)
}
