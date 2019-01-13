
mra.default <- function(oformula, odata, eformula, edata){

  dat <- init(oformula, odata,
              eformula, edata,
              scale = TRUE)

  tsr <- TSR(dat)
  wald <- Wald.test(dat)

  i <- dat$map$bet
  lm <- LM.test(dat, wald$par[i], sqrt(wald$vcov[i, i]))

  ct <- c.test(dat, wald$par, wald$vcov)

  fit <- reorganize(dat, wald, lm, ct, tsr)
  fit$call <- match.call()
  class(fit) <- 'mra'

  fit

}
