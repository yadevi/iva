

gmm <- function(dat){

  suppressWarnings(
  sol <- ucminf(dat$par.e, LS.obj, gr = LS.deriv,
                control = list(maxeval = 5000), dat = dat))

  par.e <- sol$par

  odata <- dat$odata$dat
  odata[, dat$var$var.z] <- dat$odata$og %*% sol$par[dat$map.e$alp]
  odata[, dat$var$var.g] <- NULL
  ofit <- glm(paste0(dat$var$var.d, '~.'), data = odata,
              family = 'binomial')

  par.o <- dat$par.o
  par.o[dat$map.o$a] <- coef(ofit)['(Intercept)']
  par.o[dat$map.o$bet] <- coef(ofit)[dat$var$var.z]
  if(dat$use$gam){
    par.o[dat$map.o$gam] <- coef(ofit)[dat$var$var.ox]
  }


  par <- dat$par
  par[] <- NA
  map <- dat$map
  map.e <- dat$map.e
  map.o <- dat$map.o

  var <- dat$var
  par[map$bet] <- coef(ofit)[var$var.z]
  par[map$a] <- coef(ofit)['(Intercept)'] - log(dat$n0/dat$n1)
  par[map$c0] <- par.e[map.e$c0]
  par[map$alp0] <- par.e[map.e$alp0]
  par[map$alp] <- par.e[map.e$alp]
  par[map$c1] <- par.e[map.e$c1]
  par[map$alp1] <- par.e[map.e$alp1]
  par[map$phi] <- par.e[map.e$phi]
  par[map$gam] <- coef(ofit)[var$var.ox]

  par

}