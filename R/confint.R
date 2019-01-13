

confint <- function(dat, bet, se){

  left <- seq(bet, bet - 6 * se, length.out = 20)[-1]
  right <- seq(bet, bet + 6 * se, length.out = 20)[-1]
  offset <- qchisq(.95, df = 1)

  x <- bet
  fn <- LM(bet, dat, FALSE, offset)
  lci <- -Inf
  for(b in left){
    stat <- LM(b, dat, FALSE, offset)
    fn <- c(stat, fn)
    x <- c(b, x)
    if(stat > 0){
      lci <- b
      break
    }
  }


  uci <- Inf
  for(b in right){
    stat <- LM(b, dat, FALSE, offset)
    fn <- c(fn, stat)
    x <- c(x, b)
    if(stat > 0){
      uci <- b
      break
    }
  }

  # plot(x, fn, type='b', ylim = c(min(fn), max(c(fn, 0))))
  # abline(v = bet, col = 'red')
  # abline(h = 0, col = 'green')

  if(!is.infinite(lci)){
    zero <- uniroot(LM, c(lci, bet),
                    dat = dat, return.p = FALSE, offset = offset)
    lci <- zero$root
  }

  if(!is.infinite(uci)){
    zero <- uniroot(LM, c(bet, uci),
                    dat = dat, return.p = FALSE, offset = offset)
    uci <- zero$root
  }

  ci <- c(lci, uci)
  names(ci) <- c('95% LCI', '95% UCI')

  ci

}

