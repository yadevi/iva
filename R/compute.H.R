
## H in Lemma 5
compute.H <- function(par, dat){
  
  np <- length(par)
  ret <- matrix(0, np, np)
  rownames(ret) <- names(par)
  colnames(ret) <- names(par)
  
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
  
  ret[map$c0, map$c0] <- 1/4/c0^4 * sum((1 - ed) * r0^4) - 3/4/c0^2 * sum(1 - ed)
  ret[map$c0, map$alp0] <- 1/2/c0^3 * sum((1 - ed) * r0^3)
  ret[map$alp0, map$c0] <- t(ret[map$c0, map$alp0])
  ret[map$c0, map$alp] <- 1/2/c0^3 * t(t(eg) %*% ((1 - ed) * r0^3))
  ret[map$alp, map$c0] <- t(ret[map$c0, map$alp])
  if(use$c1){
    ret[map$c1, map$c1] <- 1/4/c1^4 * sum(ed * r1^4) - 3/4/c1^2 * sum(ed)
    ret[map$c1, map$alp1] <- 1/2/c1^3 * sum(ed * r1^3)
    ret[map$alp1, map$c1] <- t(ret[map$c1, map$alp1])
    ret[map$c1, map$alp] <- 1/2/c1^3 * t(t(eg) %*% (ed * r1^3))
    ret[map$alp, map$c1] <- t(ret[map$c1, map$alp])
  }
  
  if(use$phi){
    ret[map$c0, map$phi] <- 1/2/c0^3 * t(t(ex) %*% ((1 - ed) * r0^3))
    ret[map$phi, map$c0] <- t(ret[map$c0, map$phi])
    if(use$c1){
      ret[map$c1, map$phi] <- 1/2/c1^3 * t(t(ex) %*% (ed * r0^3))
      ret[map$phi, map$c1] <- t(ret[map$c1, map$phi])
    }
  }
  
  ret
  
}


