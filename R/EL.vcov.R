# Observed information is used for Wald's test
# https://www.jstor.org/stable/pdf/2335893.pdf
#  Bradley Efron and David V. Hinkley, Biometrika, 65:3, 457-482

EL.vcov <- function(par, dat){

  J <- -deriv2(par, dat)/dat$n
  rho0 <- dat$n0/dat$n
  rho1 <- dat$n1/dat$n

  map <- dat$map
  vcov <- solve(J)
  vcov[map$a, map$a] <- vcov[map$a, map$a] - (rho0 + rho1)/rho0/rho1
  vcov <- vcov/dat$n
  colnames(vcov) <- names(par)
  rownames(vcov) <- names(par)

  if(any(eigen(vcov)$values <= 1e-6)){
    warning('EL may be misleading as vcov is not a positive-definite matrix')
  }

  vcov

}

