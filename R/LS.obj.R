

LS.obj <- function(par, dat){

  map <- dat$map.e
  use <- dat$use

  ed <- dat$edata$ed
  z <- dat$edata$z
  ex <- dat$edata$ex
  eg <- dat$edata$eg

  c0 <- par[map$c0]
  alp0 <- par[map$alp0]
  alp <- par[map$alp]

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
    return(-ell.e)
  }

  # has case in exposure data
  alp1 <- par[map$alp1]
  r1 <- r0 + alp0 - alp1
  c1 <- par[map$c1]
  ell.e <- ell.e - .5 * sum(ed * (log(c1) + r1^2/c1))

  -ell.e

}