
nderiv.null <- function(par.null, dat, bet.null){

  ## I assume map$bet == 1
  if(dat$map$bet != 1){
    stop('debug nderiv.null')
  }
  par <- c(bet = bet.null, par.null)

  nderiv(par, dat)[-1]

}


