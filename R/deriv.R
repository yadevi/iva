
deriv <- function(par, dat){

  np <- length(par)
  ret <- rep(NA, np)

  #########################
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
  pr <- 1/(n0 + n1 * delta)
  Delta <- -n1 * pr * delta

  ##################

  eta.eg <- as.vector(eg %*% alp)
  r0 <- z - alp0 - eta.eg

  if(use$phi){
    phi <- par[map$phi]
    eta.ex <- as.vector(ex %*% phi)
    r0 <- r0 - eta.ex
  }

  ##################

  ret[map$bet] <- sum((od + Delta) * eta.og)
  ret[map$a] <- sum(od + Delta)
  ret[map$c0] <- -.5 * sum((1-ed) * (1/c0 - 1/c0^2 * r0^2))
  ret[map$alp0] <- 1/c0 * sum((1-ed) * r0)
  if(use$gam){
    ret[map$gam] <- as.vector(t(ox) %*% (od + Delta))
  }

  ret[map$alp] <- bet * as.vector(t(og) %*% (od + Delta)) + 1/c0 * as.vector(t(eg) %*% ((1-ed) * r0))
  if(use$phi){
    ret[map$phi] <- 1/c0 * as.vector(t(ex) %*% ((1-ed) * r0))
  }

  ##################
  if(!use$c1){
    return(ret)
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

  ret


}


