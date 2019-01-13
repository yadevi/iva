
rescale <- function(par, vcov, map, mu, sigma){

  mat <- diag(length(par))
  rownames(mat) <- names(par)
  colnames(mat) <- names(par)
  mat[map$c0, map$c0] <- sigma^2
  mat[map$c1, map$c1] <- sigma^2
  mat[map$alp0, map$alp0] <- sigma
  mat[map$alp1, map$alp1] <- sigma
  diag(mat[map$phi, map$phi]) <- sigma
  diag(mat[map$alp, map$alp]) <- sigma
  mat[map$bet, map$bet] <- 1/sigma

  par <- as.vector(mat %*% par)
  names(par) <- colnames(vcov)
  par[map$alp0] <- par[map$alp0] + mu
  par[map$alp1] <- par[map$alp1] + mu
  vcov <- mat %*% vcov %*% t(mat)

  list(par = par, vcov = vcov)

}
