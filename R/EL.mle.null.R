

EL.mle.null <- function(dat, bet.null){

  map <- dat$map
  ## I assume map$bet == 1
  if(map$bet != 1){
    stop('debug EL.mle.null')
  }
  par.null <- dat$par[-dat$map$bet]

  fn0 <- -nloglik.null(par.null, dat, bet.null)

  suppressWarnings(
  sol <- ucminf(par.null, nloglik.null, gr = nderiv.null,
                control = list(maxeval = 5000),
                dat = dat, bet.null = bet.null))

  par <- c(bet = bet.null, sol$par)
  J <- -deriv2(par, dat)[-map$bet, -map$bet]/dat$n
  min.ev <- min(eigen(J)$values)

  if(sol$convergence > 0 && -sol$value > fn0 && min.ev > 1e-6){
    value <- nloglik.null(sol$par, dat, bet.null)
    return(list(par = par,
                value = sol$value,
                msg = sol$message,
                info = sol$info))
  }

  stop('ucminf fails')

}
