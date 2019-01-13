
print.summary.mra <- function(x, ...){

  cat('Call:\n')
  print(x$call)

  if(!is.null(x$TAB.b)){
    cat('\n')
    cat('Test for Presence of Confounders:\n')

    printCoefmat(x$TAB.b, P.values = TRUE, has.Pvalue = TRUE)
  }

  cat('\n')
  cat('Case-Control Data:\n')

  printCoefmat(x$TAB.o, P.values = TRUE, has.Pvalue = TRUE)

  cat('\n')
  cat('Exposure Data:\n')

  printCoefmat(x$TAB.e, P.values = TRUE, has.Pvalue = TRUE)


}
