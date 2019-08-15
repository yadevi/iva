
# individual-level score

ind.score <- function(par, dat){
  
  np <- length(par)
  
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
  
  #############
  
  oid <- names(od)
  eid <- names(ed)
  id <- c(oid, setdiff(eid, oid))
  sample <- matrix(0, length(id), np)
  rownames(sample) <- id
  colnames(sample) <- c(names(par))
  
  #############
  
  sample[eid, map$c0] <- -1/2 * (1 - ed) * (1/c0 - 1/c0^2 * r0^2)
  if(use$c1){
    sample[eid, map$c1] <- -1/2 * ed * (1/c1 - 1/c1^2 * r1^2)
  }
  
  sample[eid, map$alp0] <- 1/c0 * (1 - ed) * r0
  if(use$alp1){
    sample[eid, map$alp1] <- 1/c1 * ed * r1
  }
  
  sample[oid, map$alp] <- bet * (od + Delta) * og
  sample[eid, map$alp] <- sample[eid, map$alp] + 1/c0 * (1 - ed) * r0 * eg
  if(use$c1){
    sample[eid, map$alp] <- sample[eid, map$alp] + 1/c1 * ed * r1 * eg
  }
  
  sample[oid, map$a] <- (od + Delta)
  sample[oid, map$bet] <- (od + Delta) * eta.og
  
  if(use$phi){
    sample[eid, map$phi] <- 1/c0 * (1 - ed) * r0 * ex
    if(use$c1){
      sample[eid, map$phi] <- sample[eid, map$phi] + 1/c1 * ed * r1 * ex
    }
  }
  
  if(use$gam){
    sample[oid, map$gam] <- (od + Delta) * ox
  }
  
  #############
  
  sample
  
}
