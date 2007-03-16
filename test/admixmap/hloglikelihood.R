logL <- function(h, beta, d, x){
   return(h*sum(d * (log(beta) + log(x))) - sum(lgamma(h*d)))
}

N <- 100

d<-rgamma(N, shape=10, rate=10)
h<-500000
alpha<-h*d
beta<-10

x <-rgamma(N, shape=alpha, rate=10)
hvalues <-seq(1, 700000, by=10000)

y <- numeric(hvalues)
for(i in 1:length(hvalues)) {
  y[i] <- logL(hvalues[i], beta, d, x)
}

plot(hvalues, y, xlab="h", ylab="logl")
  
