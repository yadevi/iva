
summary.mra <- function(object, ...){

  map <- object$wald$map
  var <- object$wald$var
  use <- object$wald$use

  par <- coef(object)
  vcov <- vcov(object)

  id.o <- c(map$a, map$bet, map$gam)
  id.e <- c(map$alp0, map$alp1, map$alp, map$phi)
  id.sigma <- c(map$c0, map$c1)
  par.o <- par[id.o]
  par.e <- par[id.e]
  sigma2 <- par[id.sigma]

  names(par.o) <- c('(Intercept)',
                    var$var.z,
                    var$var.ox)
  if(use$alp1){
    names(par.e) <- c('(Intercept) Ctrl',
                      '(Intercept) Case',
                      var$var.g,
                      var$var.ex)
  }else{
    names(par.e) <- c('(Intercept) Ctrl',
                      var$var.g,
                      var$var.ex)
  }

  vcov.o <- vcov[id.o, id.o]
  vcov.e <- vcov[id.e, id.e]

  ##########################

  names(par.o)[2] <- paste0(names(par.o)[2], ' (Causal)')
  se.o <- sqrt(diag(vcov.o))
  tval.o <- par.o / se.o

  TAB.o <- cbind(Estimate = par.o,
               StdErr = se.o,
               t.value = tval.o,
               p.value = pchisq(tval.o^2, df = 1, lower.tail = FALSE))
  colnames(TAB.o) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')

  se.e <- sqrt(diag(vcov.e))
  tval.e <- par.e / se.e

  TAB.e <- cbind(Estimate = par.e,
                 StdErr = se.e,
                 t.value = tval.e,
                 p.value = pchisq(tval.e^2, df = 1, lower.tail = FALSE))
  colnames(TAB.e) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')

  TAB.b <- NULL
  if(!is.null(object$ct)){
    b <- object$ct$b
    se.b <- object$ct$se
    tval.b <- b / se.b
    TAB.b <- cbind(Estimate = b,
                   StdErr = se.b,
                   t.value = tval.b,
                   p.value = pchisq(tval.b^2, df = 1, lower.tail = FALSE))
    colnames(TAB.b) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
  }

  ret <- list(call = object$call,
              TAB.o = TAB.o,
              TAB.e = TAB.e,
              TAB.b = TAB.b)

  class(ret) <- 'summary.mra'
  ret

}