

observed.information <- function(par, dat){
  
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
  
  sample <- ind.score(par, dat)
  oid <- names(od)
  eid <- names(ed)
  id <- c(oid, setdiff(eid, oid))
  
  omega1 <- setdiff(oid, eid)
  omega2 <- setdiff(eid, oid)
  omega3 <- intersect(oid, eid)
  
  omega10 <- omega11 <- NULL
  omega20 <- omega21 <- NULL
  omega30 <- omega31 <- NULL
  
  has.omega1 <- ifelse(length(omega1) > 1, TRUE, FALSE)
  has.omega2 <- ifelse(length(omega2) > 1, TRUE, FALSE)
  has.omega3 <- ifelse(length(omega3) > 1, TRUE, FALSE)
  
  has.omega10 <- has.omega11 <- FALSE
  has.omega20 <- has.omega21 <- FALSE
  has.omega30 <- has.omega31 <- FALSE
  
  if(has.omega1){
    od1 <- od[omega1]
    if(any(od1 == 0)){
      omega10 <- names(od1)[od1 == 0]
      has.omega10 <- TRUE
    }
    
    if(any(od1 == 1)){
      omega11 <- names(od1)[od1 == 1]
      has.omega11 <- TRUE
    }
  }
  
  
  if(has.omega2){
    ed2 <- ed[omega2]
    if(any(ed2 == 0)){
      omega20 <- names(ed2)[ed2 == 0]
      has.omega20 <- TRUE
    }
    
    if(any(ed2 == 1)){
      omega21 <- names(ed2)[ed2 == 1]
      has.omega21 <- TRUE
    }
  }
  
  
  if(has.omega3){
    od3 <- od[omega3]
    if(any(od3 == 0)){
      omega30 <- names(od3)[od3 == 0]
      has.omega30 <- TRUE
    }
    
    if(any(od3 == 1)){
      omega31 <- names(od3)[od3 == 1]
      has.omega31 <- TRUE
    }
  }
  
  #################
  
  vid <- c(map$c0, map$alp0)
  if(use$phi){
    vid <- c(vid, map$phi)
  }
  sid <- c(omega20, omega30)
  n <- length(sid)
  ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(sample[sid, vid])
  
  if(use$c1){
    vid <- c(map$c1, map$alp1)
    if(use$phi){
      vid <- c(vid, map$phi)
    }
    sid <- c(omega21, omega31)
    n <- length(sid)
    ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(sample[sid, vid])
  }
  
  vid <- c(map$a, map$bet)
  if(use$gam){
    vid <- c(vid, map$gam)
  }
  sid <- c(omega10, omega30)
  n <- length(sid)
  ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(sample[sid, vid])
  
  sid <- c(omega11, omega31)
  n <- length(sid)
  ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(sample[sid, vid])
  
  ###############
  
  if(has.omega30){
    sid <- omega30
    n <- length(sid)
    vid1 <- map$c0
    vid2 <- map$a
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    
    vid2 <- map$bet
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    
    if(use$gam){
      vid2 <- map$gam
      ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
      ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    }
  }
  
  if(has.omega30){
    sid <- omega30
    n <- length(sid)
    vid1 <- map$alp0
    vid2 <- map$a
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    
    vid2 <- map$bet
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    
    if(use$gam){
      vid2 <- map$gam
      ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
      ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    }
  }
  
  if(has.omega30){
    sid <- omega30
    n <- length(sid)
    if(use$phi){
      vid1 <- map$phi
      vid2 <- map$a
      ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
      ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
      
      vid2 <- map$bet
      ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
      ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
      
      if(use$gam){
        vid2 <- map$gam
        ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
        ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
      }
    }
  }
  
  ###############
  
  if(has.omega31){
    sid <- omega31
    n <- length(sid)
    vid1 <- map$c1
    vid2 <- map$a
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    
    vid2 <- map$bet
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    
    if(use$gam){
      vid2 <- map$gam
      ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
      ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    }
  }
  
  if(has.omega31){
    sid <- omega31
    n <- length(sid)
    vid1 <- map$alp1
    vid2 <- map$a
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    
    vid2 <- map$bet
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    
    if(use$gam){
      vid2 <- map$gam
      ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
      ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
    }
  }
  
  if(has.omega31){
    sid <- omega31
    n <- length(sid)
    if(use$phi){
      vid1 <- map$phi
      vid2 <- map$a
      ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
      ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
      
      vid2 <- map$bet
      ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
      ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
      
      if(use$gam){
        vid2 <- map$gam
        ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], sample[sid, vid2])
        ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(sample[sid, vid2], sample[sid, vid1])
      }
    }
  }
  
  ###############
  ## alp v.s. alp
  
  sid <- c(omega10, omega30)
  n <- length(sid)
  tmp1 <- bet * (od + Delta) * og
  rownames(tmp1) <- names(od)
  vid <- map$alp
  ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(tmp1[sid, ])
  
  sid <- c(omega11, omega31)
  n <- length(sid)
  tmp1 <- bet * (od + Delta) * og
  rownames(tmp1) <- names(od)
  vid <- map$alp
  ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(tmp1[sid, ])
  
  sid <- c(omega20, omega30)
  n <- length(sid)
  tmp2 <- 1/c0 * (1 - ed) * r0 * eg
  rownames(tmp2) <- names(ed)
  vid <- map$alp
  ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(tmp2[sid, ])
  
  if(has.omega30){
    sid <- omega30
    n <- length(sid)
    ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(tmp1[sid, ], tmp2[sid, ])
    ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(tmp2[sid, ], tmp1[sid, ])
  }
  
  if(use$c1){
    sid <- c(omega21, omega31)
    n <- length(sid)
    tmp3 <- 1/c1 * ed * r1 * eg
    rownames(tmp3) <- names(ed)
    vid <- map$alp
    ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(tmp3[sid, ])
    
    if(has.omega31){
      sid <- omega31
      n <- length(sid)
      ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(tmp1[sid, ], tmp3[sid, ])
      ret[vid, vid] <- ret[vid, vid] + (n - 1) * cov(tmp3[sid, ], tmp1[sid, ])
    }
  }
  
  ii=c(map$c0,map$c1,map$alp0,map$alp1,map$alp)
  
  ###############
  ## c0 v.s. alp
  
  if(has.omega30){
    sid <- omega30
    n <- length(sid)
    vid1 <- map$c0
    vid2 <- map$alp
    tmp <- bet * (od + Delta) * og
    rownames(tmp) <- names(od)
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], tmp[sid, ])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(tmp[sid, ], sample[sid, vid1])
  }
  
  sid <- c(omega20, omega30)
  n <- length(sid)
  vid1 <- map$c0
  vid2 <- map$alp
  tmp <- 1/c0 * (1 - ed) * r0 * eg
  rownames(tmp) <- names(ed)
  ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], tmp[sid, ])
  ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(tmp[sid, ], sample[sid, vid1])
  
  ###############
  ## alp0 v.s. alp
  
  if(has.omega30){
    sid <- omega30
    n <- length(sid)
    vid1 <- map$alp0
    vid2 <- map$alp
    tmp <- bet * (od + Delta) * og
    rownames(tmp) <- names(od)
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], tmp[sid, ])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(tmp[sid, ], sample[sid, vid1])
  }
  
  sid <- c(omega20, omega30)
  n <- length(sid)
  vid1 <- map$alp0
  vid2 <- map$alp
  tmp <- 1/c0 * (1 - ed) * r0 * eg
  rownames(tmp) <- names(ed)
  ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], tmp[sid, ])
  ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(tmp[sid, ], sample[sid, vid1])
  
  print(eigen(ret)$values)
  
  ###############
  ## alp1 v.s. alp
  
  if(use$alp1){
    if(has.omega31){
      sid <- omega31
      n <- length(sid)
      vid1 <- map$alp1
      vid2 <- map$alp
      tmp <- bet * (od + Delta) * og
      rownames(tmp) <- names(od)
      ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], tmp[sid, ])
      ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(tmp[sid, ], sample[sid, vid1])
    }
    
    sid <- c(omega21, omega31)
    n <- length(sid)
    vid1 <- map$alp1
    vid2 <- map$alp
    tmp <- 1/c1 * ed * r1 * eg
    rownames(tmp) <- names(ed)
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], tmp[sid, ])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(tmp[sid, ], sample[sid, vid1])
  }
  
  print(eigen(ret)$values)
  ###############
  ## c1 v.s. alp
  
  if(use$c1){
    if(has.omega31){
      sid <- omega31
      n <- length(sid)
      vid1 <- map$c1
      vid2 <- map$alp
      tmp <- bet * (od + Delta) * og
      rownames(tmp) <- names(od)
      ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], tmp[sid, ])
      ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(tmp[sid, ], sample[sid, vid1])
    }
    
    sid <- c(omega21, omega31)
    n <- length(sid)
    vid1 <- map$c1
    vid2 <- map$alp
    tmp <- 1/c1 * ed * r1 * eg
    rownames(tmp) <- names(ed)
    ret[vid1, vid2] <- ret[vid1, vid2] + (n - 1) * cov(sample[sid, vid1], tmp[sid, ])
    ret[vid2, vid1] <- ret[vid2, vid1] + (n - 1) * cov(tmp[sid, ], sample[sid, vid1])
  }
  
  print(eigen(ret)$values)
  NULL
  
}
