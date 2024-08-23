#' @export
print.summary.streg <- function(x, digits = max(options()$digits - 4, 3),
                                  signif.stars=FALSE, ...) {

  if(is.null(digits))
    digits <- options()$digits
  cat("\nCall:\n")
  dput(x$call)

  print(x$table, digits = digits)
  if (!is.null(x$eform.table))
    print(x$eform.table, digits=digits)
  if (!is.null(x$pfixed))
    cat("\nWeibull parameter fixed at",format(x$pfixed, digits=digits),"\n")

  df  <- x$df
  cat("Log Likelihood (model)=", format(round(x$maximum[1],3)))
  cat("\n")
  if (x$robust) cat("(Log Likelihood assumes independent observations)\n")
  cat("Number of Iterations:", format(trunc(x$iterations)),
      "\n")
  omit <- x$na.action
  if (length(omit))
    cat("n=", x$n, " (", naprint(omit), ")\n", sep="")
  else cat("n=", x$n, "\n")

  invisible(NULL)
}
