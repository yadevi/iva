
diagno <- function(par, dat){

  map <- dat$map
  use <- dat$use

  eg <- dat$edata$eg
  ed <- dat$edata$ed
  z <- dat$edata$z

  alp0 <- par[map$alp0]
  alp <- par[map$alp]

  fit <- as.vector(eg %*% alp)
  if(use$phi){
    ex <- dat$edata$ex
    phi <- par[map$phi]
    fit <- fit + as.vector(ex %*% phi)
  }
  res <- as.vector(z - alp0 - fit)

  fit <- subset(fit, ed == 0)
  res <- subset(res, ed == 0)

  list(residuals = res, fitted.values = fit)

}
