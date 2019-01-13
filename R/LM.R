

LM <- function(bet.null, dat, return.p = TRUE, offset = 0, par.alt = NULL){

  sol <- EL.mle.null(dat, bet.null)
  par <- sol$par

  i <- dat$map$bet
  s <- deriv(par, dat)[i]

  if(!is.null(par.alt)){
    par <- par.alt
  }

  Jn <- Lemma4(par, dat)
  In <- Lemma5(par, dat)


  V <- In[i, i] - Jn[i, -i] %*% solve(Jn[-i, -i]) %*% In[-i, i] -
    In[i, -i] %*% solve(Jn[-i, -i]) %*% Jn[-i, i] +
    Jn[i, -i] %*% solve(Jn[-i, -i]) %*% In[-i, -i] %*% solve(Jn[-i, -i]) %*% Jn[-i, i]

  stat <- s^2/V
  pval <- pchisq(stat, df = 1, lower.tail = FALSE)

  ret <- ifelse(return.p, pval, stat - offset)

  ret

}

