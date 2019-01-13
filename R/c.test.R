
# test presence of confounder
c.test <- function(dat, par, vcov){

  use <- dat$use
  map <- dat$map

  if(!use$alp1){
    return(NULL)
  }

  alp1 <- par[map$alp1]
  alp0 <- par[map$alp0]
  bet <- par[map$bet]
  c0 <- par[map$c0]
  alp.D <- alp1 - alp0
  b <- alp.D - bet * c0
  names(b) <- 'Confounder'
  id <- c(map$alp1, map$alp0, map$bet, map$c0)

  tmp <- c(1, -1, -c0, -bet)
  se <- sqrt(t(tmp) %*% vcov[id, id] %*% tmp)
  se <- as.vector(se)
  ncp <- b^2/se^2
  power <- pchisq(qchisq(.95, df = 1),
                  df = 1, ncp = ncp, lower.tail = FALSE)
  pval <- pchisq(ncp, df = 1, lower.tail = FALSE)

  list(b = b, se = se, pval = pval)

}
