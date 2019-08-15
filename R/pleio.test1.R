

pleio.test1 <- function(dat, par){
  
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
  sample <- matrix(0, length(id), np + 1)
  rownames(sample) <- id
  colnames(sample) <- c(names(par), 'tau')
  
  #############
  
  U <- rowSums(og^2)
  sample[oid, 'tau'] <- (od + Delta) * U
  
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
  
  vcov2 <- matrix(0, np + 1, np + 1)
  rownames(vcov2) <- colnames(sample)
  colnames(vcov2) <- colnames(sample)
  
  omega1 <- setdiff(oid, eid)
  if(length(omega1) > 0){
    od1 <- od[omega1]
    if(any(od1 == 0)){
      omega10 <- names(od1)[od1 == 0]
      n10 <- length(omega10)
      vcov2 <- vcov2 + (n10 - 1) * cov(sample[omega10, ])
    }
    
    if(any(od1 == 1)){
      omega11 <- names(od1)[od1 == 1]
      n11 <- length(omega11)
      vcov2 <- vcov2 + (n11 - 1) * cov(sample[omega11, ])
    }
  }
  
  omega2 <- setdiff(eid, oid)
  if(length(omega2) > 0){
    ed2 <- ed[omega2]
    if(any(ed2 == 0)){
      omega20 <- names(ed2)[ed2 == 0]
      n20 <- length(omega20)
      vcov2 <- vcov2 + (n20 - 1) * cov(sample[omega20, ])
    }
    
    if(any(ed2 == 1)){
      omega21 <- names(ed2)[ed2 == 1]
      n21 <- length(omega21)
      vcov2 <- vcov2 + (n21 - 1) * cov(sample[omega21, ])
    }
  }
  
  omega3 <- intersect(oid, eid)
  if(length(omega3) > 0){
    od3 <- od[omega3]
    if(any(od3 == 0)){
      omega30 <- names(od3)[od3 == 0]
      n30 <- length(omega30)
      vcov2 <- vcov2 + (n30 - 1) * cov(sample[omega30, ])
    }
    
    if(any(od3 == 1)){
      omega31 <- names(od3)[od3 == 1]
      n31 <- length(omega31)
      vcov2 <- vcov2 + (n31 - 1) * cov(sample[omega31, ])
    }
  }
  
  #############
  #############
  ## S = (D + Delta) * U
  ## S is the score statistic for the test
  
  S <- t(U) %*% (od + Delta)
  
  ##########################################################
  ## compute cov(S) at parameters estimated under the null #
  ##########################################################
  
  ## zeta = (a, gam, bet, alp)
  ## zeta is the parameters used in defining S
  ## zeta is a sub-vector of par
  ## gam refers to mu in the note, the coefficients of covariates in outcome model
  
  zeta.id <- c(map$a, map$gam, map$bet, map$alp)
  
  ## compute dS = partial S/partial zeta
  if(use$gam){
    V <- as.matrix(data.frame(1, ox, og %*% alp, bet * og))
  }else{
    V <- as.matrix(data.frame(1, og %*% alp, bet * og))
  }
  colnames(V) <- NULL
  
  dS <- t(U) %*% (xi * V)
  
  ############
  
  H <- solve(-deriv2(par, dat))[zeta.id, ]
  
  mat <- cbind(dS %*% H, 1)
  vcov <- mat %*% vcov2 %*% t(mat)
  stat <- as.vector(t(S) %*% solve(vcov) %*% S)
  
  pval <- pchisq(stat, df = 1, lower.tail = FALSE)
  
  list(pval = pval, S = S, vcov = vcov)
  
}
