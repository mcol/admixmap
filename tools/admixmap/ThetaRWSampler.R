softmax <- function(a){
  mu <- exp(a) / sum(exp(a))
  return(mu)
}
inv.softmax <- function(mu){
  a <- log(mu) - sum(log(mu))/length(mu)
  return(a)
}
##dirichletLogDensity <- function(x, a){
  ##K <- length(a)
  ##theta <- c(x[1:(K-1)], 1-sum(x[1:(K-1)]))
  ##f <- lgamma(sum(a))
  ##for(i in 1:K){
    ##if(a[i]>0)
      ##f <- f + ( a[i] - 1 ) * log( theta[i] ) - lgamma( a[i] )
  ##}
  ##f              
##}
dirichletLogDensity <- function(x, alpha){
  K <- length(a)
  theta <- c(x[1:(K-1)], 1-sum(x[1:(K-1)]))
  
  f <- lgamma(sum(alpha))
  for(i in 1:K){
    if(alpha[i]>0)
      f <- f + ( alpha[i] ) * log( theta[i] ) - lgamma( alpha[i] )
  }
  f              
}


stepsize <- 1.0
theta <- c(0.5,0.5)
sum.theta <- theta
N.iters <- 10000
plot(1, xlim=c(0,N.iters), ylim = c(0, 1), type='n')
avec <- 0
thetavec <- 0.5
points(0, theta[1])
##points(0, 0)

n <- 0.1 ## proportion of allele1
alpha <- c(1+n, 24-n)
n.accepted <- 0
for( it in 1:N.iters){
old <- theta[1]
  a <- inv.softmax(theta)
  ##old <- a[1]
  ##new <- old
  a[1] <- rnorm(1, a[1], stepsize)
  a[2] <- 0-a[1]
  theta.prop <- softmax(a)
##u <- ##rnorm(1, theta[1], 0.1)
##runif(1)
  ##rbeta(1,1,1)
##theta.prop <- c(u, 1-u)


  logp <- dirichletLogDensity(theta.prop, alpha) - dirichletLogDensity(theta, alpha)
  if(logp >=0 ){accept <- TRUE
  }else{
    if(log(runif(1)) < logp){accept <- TRUE;
    }else accept <- FALSE
  }
  if(accept){
    n.accepted <- n.accepted+1;
    theta <- theta.prop;
    ##new <- a[1];
  }
  sum.theta <- sum.theta+theta
  ##points(it, theta[1])
  lines(x=c(it-1, it), y=c(old, theta[1]))
  avec <- c(avec,new)
thetavec <- c(thetavec, theta)

}
##plot(c(1:(N.iters+1)), avec, type='l')
print("Acceptance rate = ")
print(n.accepted / N.iters)
print("Ergodic averages: ")
print(sum.theta / N.iters)
      
