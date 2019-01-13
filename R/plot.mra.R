
plot.mra <- function(x, ...){

  resid <- resid(x)
  fitted <- fitted(x)

  plot(fitted, resid, type = 'p',
       xlab = 'Fitted values', ylab = 'Residuals',
       main = 'Residuals vs Fitted')
  abline(h = 0, lty = 'dotted', col = 'grey')
  lines(lowess(fitted, resid), col = 'red')

}
