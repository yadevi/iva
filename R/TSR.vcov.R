
TSR.vcov <- function(par, dat){

  J <- -TSR.hess(par, dat)/dat$n

  map <- dat$map
  id <- c(map$a, map$bet, map$gam)
  N <- J[-id, -id]
  C <- J[id, id]
  B <- J[id, -id]

  rho0 <- dat$n0/dat$n
  rho1 <- dat$n1/dat$n

  I <- J
  I[] <- 0
  I[-id, -id] <- N
  C1 <- C[, 'a', drop = FALSE]
  I[id, id] <- C - (rho0 + rho1)/rho0/rho1 * C1 %*% t(C1)

  vcov <- solve(J) %*% I %*% t(solve(J))/dat$n

  if(any(eigen(vcov)$values <= 1e-6)){
    warning('TSR may be misleading as vcov is not a positive-definite matrix')
  }

  vcov

}
