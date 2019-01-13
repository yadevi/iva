

nloglik.null <- function(par.null, dat, bet.null){

  ## I assume map$bet == 1
  if(dat$map$bet != 1){
    stop('debug nloglik.null')
  }
  par <- c(bet = bet.null, par.null)

  nloglik(par, dat)

}


