

EL.mle <- function(dat, max.iter = 1){

  par <- dat$par
  fn0 <- loglik(par, dat)
  iter <- 0
  while (TRUE) {
    if(iter > max.iter){
      stop('ucminf fails')
    }

    sol <- ucminf(par, nloglik, gr = nderiv,
                  control = list(maxeval = 5000), dat = dat)

    if(sol$convergence %in% c(1,2,4) && -sol$value > fn0 && sol$info['maxgradient'] < 1e-5){
      J <- -deriv2(sol$par, dat)/dat$n
      min.ev <- min(eigen(J)$values)
      if(min.ev > 1e-6){
        return(list(par = sol$par,
                    value = sol$value,
                    msg = sol$message,
                    info = sol$info))
      }
    }

    par[dat$map$bet] <- runif(1, -.2, .2)
    iter <- iter + 1
  }

  stop('fail to find an estimate')

}
