

Lemma5 <- function(par, dat){

  map <- dat$map
  n0 <- dat$n0
  n1 <- dat$n1

  Jn <- Lemma4(par, dat)
  Jna <- Jn[, map$a, drop = FALSE]
  Hn <- compute.H(par, dat)
  In <- Jn - (n0 + n1)/n0/n1 * Jna %*% t(Jna)
  # In <- In + Hn

  In

}
