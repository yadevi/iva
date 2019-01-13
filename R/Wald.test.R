

Wald.test <- function(dat){

  sol <- EL.mle(dat)
  suppressWarnings(vcov <- EL.vcov(sol$par, dat))

  pval.wald <- EL.wald(sol$par, vcov)

  diag <- diagno(sol$par, dat)

  par <- sol$par
  sol$par <- NULL

  i <- dat$map$bet
  ncp <- par[i]^2/vcov[i, i]
  power <- pchisq(qchisq(.95, df = 1),
                  df = 1, ncp = ncp, lower.tail = FALSE)

  wald <- list(par = par,
               vcov = vcov,
               pval = pval.wald,
               sol = sol,
               diag = diag)

}