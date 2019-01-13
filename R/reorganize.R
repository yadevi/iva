
reorganize <- function(dat, wald, lm, ct, tsr){

  var <- dat$var
  map <- dat$map
  use <- dat$use

  if(use$shared.x){
    np <- length(wald$par)
    tmp <- diag(np)
    rownames(tmp) <- rownames(wald$vcov)
    colnames(tmp) <- colnames(wald$vcov)
    tmp[map$gam.x, map$bet] <- -wald$par[map$phi.x]
    tmp[map$gam.x, map$phi.x] <- wald$par[map$bet]

    wald$vcov <- tmp %*% wald$vcov %*% t(tmp)
    wald$par[map$gam.x] <- wald$par[map$gam.x] - wald$par[map$bet] * wald$par[map$phi.x]

    tmp[map$gam.x, map$bet] <- -tsr$par[map$phi.x]
    tmp[map$gam.x, map$phi.x] <- tsr$par[map$bet]

    tsr$vcov <- tmp %*% tsr$vcov %*% t(tmp)
    tsr$par[map$gam.x] <- tsr$par[map$gam.x] - tsr$par[map$bet] * tsr$par[map$phi.x]
  }

  wald$map <- map
  wald$var <- var
  wald$use <- use

  par <- wald$par
  vcov <- wald$vcov
  residuals <- wald$diag$residuals
  fitted.values <- wald$diag$fitted.values

  mu <- dat$mu
  sigma <- dat$sigma
  residuals <- residuals * sigma
  fitted.values <- fitted.values * sigma

  wald$par <- NULL
  wald$vcov <- NULL
  wald$diag <- NULL


  ###############

  tmp <- rescale(par, vcov, map, mu, sigma)
  par <- tmp$par
  vcov <- tmp$vcov

  sigma2 <- par[c(map$c0, map$c1)]

  lm$ci <- lm$ci / sigma

  tmp <- rescale(tsr$par, tsr$vcov, map, mu, sigma)
  tsr$par <- tmp$par
  tsr$vcov <- tmp$vcov

  # par and vcov returned are different from the estimated parameters
  # if shared covariates are used in both models.
  # returned par and its vcov are subject to parameterization used in underlying models
  # see Assumptions 1 and 2 for underlying models and their parameters
  fit <- list(sigma2 = sigma2,
              residuals = residuals, fitted.values = fitted.values,
              coefficients = par, vcov = vcov,
              wald = wald, lm = lm, ct = ct, tsr = tsr)

}