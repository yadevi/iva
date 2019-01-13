
LS.deriv <- function(par, dat){

  np <- length(par)
  ret <- rep(NA, np)

  #########################
  map <- dat$map.e
  use <- dat$use

  ed <- dat$edata$ed
  z <- dat$edata$z
  ex <- dat$edata$ex
  eg <- dat$edata$eg

  c0 <- par[map$c0]
  alp0 <- par[map$alp0]
  alp <- par[map$alp]

  ##################

  eta.eg <- as.vector(eg %*% alp)
  r0 <- z - alp0 - eta.eg

  if(use$phi){
    phi <- par[map$phi]
    eta.ex <- as.vector(ex %*% phi)
    r0 <- r0 - eta.ex
  }

  ##################

  ret[map$c0] <- -.5 * sum((1-ed) * (1/c0 - 1/c0^2 * r0^2))
  ret[map$alp0] <- 1/c0 * sum((1-ed) * r0)

  ret[map$alp] <- 1/c0 * as.vector(t(eg) %*% ((1-ed) * r0))
  if(use$phi){
    ret[map$phi] <- 1/c0 * as.vector(t(ex) %*% ((1-ed) * r0))
  }

  ##################
  if(!use$c1){
    return(-ret)
  }

  # has case in exposure data
  alp1 <- par[map$alp1]
  r1 <- r0 + alp0 - alp1
  c1 <- par[map$c1]

  ret[map$c1] <- -.5 * sum(ed * (1/c1 - 1/c1^2 * r1^2))
  ret[map$alp1] <- 1/c1 * sum(ed * r1)

  ret[map$alp] <- ret[map$alp] + 1/c1 * as.vector(t(eg) %*% (ed * r1))
  if(use$phi){
    ret[map$phi] <- ret[map$phi] + 1/c1 * as.vector(t(ex) %*% (ed * r1))
  }

  -ret


}


