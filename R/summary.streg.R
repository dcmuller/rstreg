#' @export
summary.streg <- function(object, ...)
{
  if (!is.null(object$fail)) {
    warning(" streg failed.", object$fail, "   No summary provided\n")
    return(invisible(object))
  }

  coef <- object$par
  cname <- names(object$par)

  n <- length(object$linear.predictors)
  p <- sum(!is.na(coef))
  if(!p) {
    warning("This model has zero rank --- no summary is provided")
    return(invisible(object))
  }

  if (is.null(object$naive.var)){
    table <- matrix(rep(coef, 6), ncol = 6)
    dimnames(table) <- list(cname, c("Value", "Std. Error", "z", "p", "lci.95", "uci.95"))
    stds <- sqrt(diag(object$var))
    table[, 2] <- stds
    table[, 3] <- table[, 1]/stds
    table[, 4] <- 2*pnorm(-abs(table[,3]))
    table[, 5] <- table[, 1] + qnorm(0.025)*stds
    table[, 6] <- table[, 1] + qnorm(0.975)*stds
  } else {
    table <- matrix(rep(coef, 7), ncol = 7)
    dimnames(table) <- list(cname, c("Value", "Std. Err","(Naive SE)", "z", "p", "lci.95", "uci.95"))
    stds <- sqrt(diag(object$var))
    table[, 2] <- stds
    table[, 3] <- sqrt(diag(object$naive.var))
    table[, 4] <- table[, 1]/stds
    table[, 5] <- 2*pnorm(-abs(table[,4]))
    table[, 6] <- table[, 1] + qnorm(0.025)*stds
    table[, 7] <- table[, 1] + qnorm(0.975)*stds
  }
  if (object$metric=="PH"){
    eform.table <- exp(table[, c("Value", "lci.95", "uci.95")])
    colnames(eform.table)[1] <- "Haz.Ratio"
  }

  x <- object[match(c('call', 'df', 'value', 'iterations', 'na.action',
                      'coefficients', 'estimates', 'var', 'pfixed'),
                    names(object), nomatch=0)]
  x <- c(x, list(table=table, eform.table=eform.table, n=n),
         robust=!is.null(object$naive.var))

  class(x) <- 'summary.streg'
  x
}
