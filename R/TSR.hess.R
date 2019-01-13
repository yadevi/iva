
TSR.hess <- function(par, dat){

  np <- length(par)
  ret <- matrix(0, np, np)

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
  xi <- Delta * (1 + Delta)

  ##################

  eta.eg <- as.vector(eg %*% alp)
  r0 <- z - alp0 - eta.eg

  if(use$phi){
    phi <- par[map$phi]
    eta.ex <- as.vector(ex %*% phi)
    r0 <- r0 - eta.ex
  }

  ##################

  if(use$c1){# has case in exposure data
    alp1 <- par[map$alp1]
    r1 <- r0 + alp0 - alp1
    c1 <- par[map$c1]
  }

  ##################

  ret[map$c0, map$c0] <- -.5 /c0^2 * sum(1 - ed)

  ##################

  if(use$c1){
    ret[map$c1, map$c1] <- -.5/c1^2 * sum(ed)
  }

  ##################

  ret[map$alp0, map$alp0] <- -1/c0 * sum(1 - ed)
  ret[map$alp, map$alp0] <- -1/c0 * (t(eg) %*% (1 - ed))
  ret[map$alp0, map$alp] <- t(ret[map$alp, map$alp0])
  if(use$phi){
    ret[map$phi, map$alp0] <- -1/c0 * (t(ex) %*% (1 - ed))
    ret[map$alp0, map$phi] <- t(ret[map$phi, map$alp0])
  }

  ##################

  if(use$alp1){
    ret[map$alp1, map$alp1] <- -1/c1 * sum(ed)
    ret[map$alp, map$alp1] <- -1/c1 * (t(eg) %*% ed)
    ret[map$alp1, map$alp] <- t(ret[map$alp, map$alp1])
    if(use$phi){
      ret[map$phi, map$alp1] <- -1/c1 * (t(ex) %*% ed)
      ret[map$alp1, map$phi] <- t(ret[map$phi, map$alp1])
    }
  }

  ##################

  ret[map$alp, map$alp] <- - 1/c0 * (t(eg) %*% ((1 - ed)* eg))
  if(use$c1){
    ret[map$alp, map$alp] <- ret[map$alp, map$alp] - 1/c1 * t(eg) %*% (ed * eg)
  }
  ret[map$a, map$alp] <- bet * n0 * t(t(og) %*% (Delta * pr))

  ret[map$bet, map$alp] <- bet * n0 * t(t(og) %*% (Delta * eta.og * pr))

  if(use$phi){
    ret[map$alp, map$phi] <- -1/c0 * (t(eg) %*% ((1 - ed) * ex))
    if(use$c1){
      ret[map$alp, map$phi] <- ret[map$alp, map$phi] - 1/c1 * (t(eg) %*% (ed * ex))
    }
    ret[map$phi, map$alp] <- t(ret[map$alp, map$phi])
  }
  if(use$gam){
    ret[map$gam, map$alp] <- bet * n0 * (t(ox) %*% ((Delta * pr) * og))
  }

  ##################

  ret[map$a, map$a] <- n0 * sum(Delta * pr)
  ret[map$a, map$bet] <- n0 * sum(Delta * eta.og * pr)
  ret[map$bet, map$a] <- ret[map$a, map$bet]
  if(use$gam){
    ret[map$gam, map$a] <- n0 * (t(ox) %*% (Delta * pr))
    ret[map$a, map$gam] <- t(ret[map$gam, map$a])
  }

  ##################

  ret[map$bet, map$bet] <- n0 * sum(Delta * eta.og^2 * pr)
  if(use$gam){
    ret[map$gam, map$bet] <- n0 * (t(ox) %*% (Delta * eta.og * pr))
    ret[map$bet, map$gam] <- t(ret[map$gam, map$bet])
  }

  ##################

  if(use$phi){
    ret[map$phi, map$phi] <- -1/c0 * (t(ex) %*% ((1 - ed) * ex))
    if(use$c1){
      ret[map$phi, map$phi] <- ret[map$phi, map$phi] - 1/c1 * (t(ex) %*% (ed * ex))
    }
  }

  ##################

  if(use$gam){
    ret[map$gam, map$gam] <- n0 * (t(ox) %*% ((Delta * pr) * ox))
  }

  rownames(ret) <- names(par)
  colnames(ret) <- names(par)

  ret


}


