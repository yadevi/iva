
TSR <- function(dat){

  dat <- TSR.init(dat)
  par <- gmm(dat)
  suppressWarnings(vcov <- TSR.vcov(par, dat))

  list(par = par, vcov = vcov)

}
