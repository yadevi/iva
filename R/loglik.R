

loglik <- function(par, dat){

  map <- dat$map
  use <- dat$use

  od <- dat$odata$od
  ox <- dat$odata$ox
  og <- dat$odata$og

  ed <- dat$edata$ed
  z <- dat$edata$z
  ex <- dat$edata$ex
  eg <- dat$edata$eg

  bet <- par[map$bet]
  a <- par[map$a]
  c0 <- par[map$c0]
  alp0 <- par[map$alp0]
  alp <- par[map$alp]

  n0 <- sum(1-od)
  n1 <- sum(od)
  eta.og <- as.vector(og %*% alp)
  tmp <- a + bet * eta.og
  if(use$gam){
    gam <- par[map$gam]
    eta.ox <- as.vector(ox %*% gam)
    tmp <- tmp + eta.ox
  }
  delta <- exp(tmp)
  ell.o <- sum(od * tmp) - sum(log(n0 + n1 * delta))

  ############

  eta.eg <- as.vector(eg %*% alp)
  r0 <- z - alp0 - eta.eg

  if(use$phi){
    phi <- par[map$phi]
    eta.ex <- as.vector(ex %*% phi)
    r0 <- r0 - eta.ex
  }

  ell.e <- - .5 * sum((1-ed) * (log(c0) + r0^2/c0))

  #############################

  if(!use$c1){
    return(ell.o + ell.e)
  }

  # has case in exposure data
  alp1 <- par[map$alp1]
  r1 <- r0 + alp0 - alp1
  c1 <- par[map$c1]
  ell.e <- ell.e - .5 * sum(ed * (log(c1) + r1^2/c1))

  ell.o + ell.e

}

