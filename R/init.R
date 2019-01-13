

init <- function(oformula, odata,
                 eformula, edata,
                 scale = TRUE){

  n <- length(unique(c(odata$id, edata$id)))

  po <- parse.odata(oformula, odata)
  pe <- parse.edata(eformula, edata)

  if(po$var.d != pe$var.d){
    stop('Different outcome detected in two datasets')
  }

  if(!setequal(po$var.g, pe$var.g)){
    stop('Different instruments detected in two datasets')
  }

  var.d <- po$var.d
  var.z <- pe$var.z
  var.g <- po$var.g
  var.ox <- po$var.ox
  var.ex <- pe$var.ex

  odata <- po$odata
  edata <- pe$edata

  od <- odata[, var.d]
  if(!is.null(var.ox)){
    ox <- as.matrix(odata[, var.ox, drop = FALSE])
  }else{
    ox <- NULL
  }
  og <- as.matrix(odata[, var.g, drop = FALSE])

  ed <- edata[, var.d]

  # scale the exposure for numerical stability
  if(scale){
    # message('mra standarizes the exposure to have zero mean and unit variance for numerical stability')
    mu <- mean(edata[, var.z])
    sigma <- sd(edata[, var.z])
    edata[, var.z] <- scale(edata[, var.z])
  }else{
    mu <- 0
    sigma <- 1
  }

  z <- edata[, var.z]
  if(!is.null(var.ex)){
    ex <- as.matrix(edata[, var.ex, drop = FALSE])
  }else{
    ex <- NULL
  }
  eg <- as.matrix(edata[, var.g, drop = FALSE])

  #########################

  ofit <- glm(paste0(var.d, '~.'),
              data = odata[, c(var.d, var.ox, var.g)],
              family = 'binomial')
  par.og <- coef(ofit)[var.g]

  efit0 <- lm(paste0(var.z, '~.'),
              data = edata[ed == 0, c(var.z, var.ex, var.g)])

  if(any(ed == 1)){
    efit1 <- lm(paste0(var.z, '~.'),
                data = edata[ed == 1, c(var.z, var.ex, var.g)])
  }

  ################################

  k <- length(var.g) # number of IV
  np <- 4 + k + length(var.ox) + length(var.ex) + 2 * any(ed == 1)
  par <- rep(NA, np)
  n0 <- sum(1 - od)
  n1 <- sum(od)

  ######################

  map <- list(bet = 1,
              a = 2,
              c0 = 3,
              alp0 = 4,
              alp = 5:(4+k))

  par[map$bet] <- 0
  par[map$a] <- coef(ofit)['(Intercept)'] - log(n1/n0)
  par[map$c0] <- summary(efit0)$sigma^2
  par[map$alp0] <- coef(efit0)['(Intercept)']
  par[map$alp] <- coef(efit0)[var.g]

  names(par)[map$bet] <- 'bet'
  names(par)[map$a] <- 'a'
  names(par)[map$c0] <- 'c0'
  names(par)[map$alp0] <- 'alp0'
  names(par)[map$alp] <- paste0('alp.', var.g)

  use <- list(bet = TRUE,
              a = TRUE,
              c0 = TRUE,
              alp0 = TRUE,
              alp = TRUE,
              c1 = FALSE,
              alp1 = FALSE,
              phi = FALSE,
              gam = FALSE,
              shared.x = FALSE)

  last.loc <- 4+k

  if(any(ed == 1)){
    map$c1 <- last.loc + 1
    map$alp1 <- last.loc + 2
    use$c1 <- TRUE
    use$alp1 <- TRUE
    last.loc <- last.loc + 2

    par[map$c1] <- summary(efit1)$sigma^2
    par[map$alp1] <- coef(efit1)['(Intercept)']

    names(par)[map$c1] <- 'c1'
    names(par)[map$alp1] <- 'alp1'
  }

  ek <- length(var.ex)
  if(ek > 0){
    map$phi <- (last.loc + 1):(last.loc + ek)
    use$phi <- TRUE
    last.loc <- last.loc + ek

    par[map$phi] <- coef(efit0)[var.ex]
    names(par)[map$phi] <- paste0('phi.', var.ex)
  }

  ok <- length(var.ox)
  if(ok > 0){
    map$gam <- (last.loc + 1):(last.loc + ok)
    use$gam <- TRUE
    last.loc <- last.loc + ok

    par[map$gam] <- coef(ofit)[var.ox]
    names(par)[map$gam] <- paste0('gam.', var.ox)
  }

  if(ek > 0 && ok > 0){
    var.x <- intersect(var.ox, var.ex)
    tmp <- 1:length(par)
    names(tmp) <- names(par)

    map$gam.x <- tmp[paste0('gam.', var.x)]
    map$phi.x <- tmp[paste0('phi.', var.x)]
    use$shared.x <- TRUE
  }

  ################################

  n0 <- sum(1 - od)
  n1 <- sum(od)

  odata <- list(od = od,
                ox = ox,
                og = og,
                dat = odata[, c(var.d, var.ox, var.g)])

  edata <- list(ed = ed,
                z = z,
                ex = ex,
                eg = eg,
                dat0 = edata[ed == 0, c(var.z, var.ex, var.g)],
                dat1 = edata[ed == 1, c(var.z, var.ex, var.g)])

  var <- list(var.d = var.d,
              var.z = var.z,
              var.g = var.g,
              var.ex = var.ex,
              var.ox = var.ox)

  list(n = n, n0 = n0, n1 = n1,
       odata = odata, edata = edata,
       map = map, use = use,
       par = par, var = var,
       mu = mu, sigma = sigma)

}

