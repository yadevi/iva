
Lemma4 <- function(par, dat){

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

  ret[map$c0, map$c0] <- -.5/c0^2 * sum(1 - ed)

  if(use$c1){
    ret[map$c1, map$c1] <- -.5/c1^2 * sum(ed)
  }

  ret[map$alp0, map$alp0] <- -1/c0 * sum(1 - ed)
  ret[map$alp0, map$alp] <- -1/c0 * t(t(eg) %*% (1 - ed))
  ret[map$alp, map$alp0] <- t(ret[map$alp0, map$alp])
  if(use$phi){
    ret[map$alp0, map$phi] <- -1/c0 * t(t(ex) %*% (1 - ed))
    ret[map$phi, map$alp0] <- t(ret[map$alp0, map$phi])
  }

  if(use$alp1){
    ret[map$alp1, map$alp1] <- -1/c1 * sum(ed)
    ret[map$alp1, map$alp] <- -1/c1 * t(t(eg) %*% ed)
    ret[map$alp, map$alp1] <- t(ret[map$alp1, map$alp])
    if(use$phi){
      ret[map$alp1, map$phi] <- -1/c1 * t(t(ex) %*% ed)
      ret[map$phi, map$alp1] <- t(ret[map$alp1, map$phi])
    }
  }

  ret[map$alp, map$alp] <- bet^2 * n0 * (t(og) %*% ((Delta * pr) * og)) -
    1/c0 * (t(eg) %*% ((1 - ed)* eg))
  if(use$c1){
    ret[map$alp, map$alp] <- ret[map$alp, map$alp] - 1/c1 * t(eg) %*% (ed * eg)
  }
  ret[map$alp, map$a] <- bet * n0 * (t(og) %*% (Delta * pr))
  ret[map$a, map$alp] <- t(ret[map$alp, map$a])
  ret[map$alp, map$bet] <- bet * n0 * (t(og) %*% (Delta * eta.og * pr))
  ret[map$bet, map$alp] <- t(ret[map$alp, map$bet])
  if(use$phi){
    ret[map$alp, map$phi] <- -1/c0 * (t(eg) %*% ((1 - ed) * ex))
    if(use$c1){
      ret[map$alp, map$phi] <- ret[map$alp, map$phi] - 1/c1 * (t(eg) %*% (ed * ex))
    }
    ret[map$phi, map$alp] <- t(ret[map$alp, map$phi])
  }
  if(use$gam){
    ret[map$alp, map$gam] <- bet * n0 * t(og) %*% (Delta * pr * ox)
    ret[map$gam, map$alp] <- t(ret[map$alp, map$gam])
  }

  ret[map$a, map$a] <- n0 * sum(pr * Delta)
  ret[map$a, map$bet] <- n0 * sum(pr * Delta * eta.og)
  ret[map$bet, map$a] <- ret[map$a, map$bet]
  if(use$gam){
    ret[map$a, map$gam] <- n0 * t(t(ox) %*% (Delta * pr))
    ret[map$gam, map$a] <- t(ret[map$a, map$gam])
  }

  ret[map$bet, map$bet] <- n0 * sum(pr * Delta * eta.og^2)
  if(use$gam){
    ret[map$bet, map$gam] <- n0 * t(t(ox) %*% (Delta * eta.og * pr))
    ret[map$gam, map$bet] <- t(ret[map$bet, map$gam])
  }

  if(use$phi){
    ret[map$phi, map$phi] <- -1/c0 * (t(ex) %*% ((1 - ed) * ex))
    if(use$c1){
      ret[map$phi, map$phi] <- ret[map$phi, map$phi] - 1/c1 * (t(ex) %*% (ed * ex))
    }
  }

  if(use$gam){
    ret[map$gam, map$gam] <- n0 * (t(ox) %*% (Delta * pr * ox))
  }

  Jn <- -ret
  Jn

}

