

LM.test <- function(dat, bet, se){

  pval <- LM(0, dat, TRUE)

  ci <- confint(dat, bet, se)

  list(pval = pval, ci = ci)

}

