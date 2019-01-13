
TSR.init <- function(dat){

  use <- dat$use
  var <- dat$var
  var.d <- var$var.d
  var.z <- var$var.z
  var.g <- var$var.g
  var.ex <- var$var.ex
  var.ox <- var$var.ox

  odata <- dat$odata
  edata <- dat$edata

  ofit <- glm(paste0(var.d, '~.'),
              data = odata$dat,
              family = 'binomial')

  efit0 <- lm(paste0(var.z, '~.'),
              data = subset(edata$dat0))

  if(use$alp1){
    efit1 <- lm(paste0(var.z, '~.'),
                data = subset(edata$dat1))
  }

  ################################

  k <- length(var.g) # number of IV
  np <- 2 + k + length(var.ex) + 2 * (use$alp1)
  par.e <- rep(NA, np)

  ######################

  map.e <- list(c0 = 1,
                alp0 = 2,
                alp = 3:(2+k))

  par.e[map.e$c0] <- summary(efit0)$sigma^2
  par.e[map.e$alp0] <- coef(efit0)['(Intercept)']
  par.e[map.e$alp] <- coef(efit0)[var.g]

  names(par.e)[map.e$c0] <- 'c0'
  names(par.e)[map.e$alp0] <- 'alp0'
  names(par.e)[map.e$alp] <- paste0('alp.', var.g)

  last.loc <- 2+k

  if(use$alp1){
    map.e$c1 <- last.loc + 1
    map.e$alp1 <- last.loc + 2
    last.loc <- last.loc + 2

    par.e[map.e$c1] <- summary(efit1)$sigma^2
    par.e[map.e$alp1] <- coef(efit1)['(Intercept)']

    names(par.e)[map.e$c1] <- 'c1'
    names(par.e)[map.e$alp1] <- 'alp1'
  }

  ek <- length(var.ex)
  if(ek > 0){
    map.e$phi <- (last.loc + 1):(last.loc + ek)
    use$phi <- TRUE
    last.loc <- last.loc + ek

    par.e[map.e$phi] <- coef(efit0)[var.ex]
    names(par.e)[map.e$phi] <- paste0('phi.', var.ex)
  }

  map.o <- list(a = 1,
                bet = 2)
  np <- 2 + length(var.ox)
  par.o <- rep(NA, np)
  names(par.o)[1:2] <- c('a', 'bet')
  if(use$gam){
    map.o$gam <- 3:np
    names(par.o)[map.o$gam] <- paste0('gam.', var.ox)
  }

  n <- dat$n
  n0 <- dat$n0
  n1 <- dat$n1

  list(n = n, n0 = n0, n1 = n1,
       odata = odata, edata = edata,
       map.o = map.o, map.e = map.e, map = dat$map,
       use = use, var = var,
       par.o = par.o, par.e = par.e, par = dat$par)

}

