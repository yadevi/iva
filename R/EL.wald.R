
EL.wald <- function(par, vcov){

  pval.wald <- pchisq(par['bet']^2/vcov['bet', 'bet'], df = 1, lower.tail = FALSE)
  pval.wald

}

