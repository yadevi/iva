

pleio.test <- function(dat, par, vcov){
  
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
  #############
  ## S = (D + Delta) * U
  ## S is the score statistic for the test
  U1 <- scale(rowSums(og))
  U2 <- scale(rowSums(og^2))
  U <- matrix(c(U1, U2), ncol = 2, byrow = FALSE)
  
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
  #dS <- n0 * (t(U) %*% (pr * Delta * V))
  
  ############
  
  id.c <- map$c0
  if(use$c1){
    id.c <- c(id.c, map$c1)
  }
  
  Jn <- -deriv2(par, dat)
  zeta.name <- colnames(Jn)[zeta.id]
  H <- solve(Jn[-id.c, -id.c])[zeta.name, ]
  
  mat1 <- cbind(dS %*% H, diag(2))
  mat2 <- -deriv3(par, dat)
  mat2 <- mat2 - (n0 + n1)/n0/n1 * mat2[, map$a, drop = FALSE] %*% mat2[map$a, , drop =FALSE]
  mat2 <- mat2[-id.c, -id.c]
  
  vcov2 <- mat1 %*% mat2 %*% t(mat1)
  
  #############
  
  stat2 <- as.vector(t(S) %*% solve(vcov2) %*% S)
  pval2 <- pchisq(stat2, df = 2, lower.tail = FALSE)
  
  list(p=pval2, S=S, vcov=vcov2)
  
}
