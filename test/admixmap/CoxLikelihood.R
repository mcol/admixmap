logl.OneIndiv <- function(i, j, x, beta, n, r, alpha, d){
  xbeta = x%*%beta  
  l <- x[i,j]*beta[j]* sum(n[i,]*r[i,]) - exp(xbeta[i])*sum(n[i,]*alpha*d)
  return(l)
}

logl <- function(x, beta, j, n, r, alpha, d){
  logl <- 0
  for(i in 1:nrow(x)){
    logl <- logl + logl.OneIndiv(i, j, x, beta, n, r, alpha, d)
  }
  return (logl)
}
logl(covariate, c(0,0), 2, n, r, alpha, d)
