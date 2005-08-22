## script to simulate test data for admixmap
simulateHaploidAlleles <- function(M,rho,x,L) {   
  gameteAncestry <- numeric(L)  
  randAnc <- runif(L)  
  gameteAncestry[1] <- ifelse(randAnc[1] > M, 2, 1) # M is prhob pop 1.
  for(locus in 2:L) {    
    T <- t(matrix(data=c(M + (1-M)*exp(-rho*x[locus]),  1-M - (1-M)*exp(-rho*x[locus]),
                    M - M*exp(-rho*x[locus]), 1-M + M*exp(-rho*x[locus]) ),
                  nrow=2, ncol=2))
    gameteAncestry[locus] <- ifelse(randAnc[locus] > T[gameteAncestry[locus-1], 1], 2, 1)
  }
  simulateHaploidAlleles <-
    ifelse(runif(L) > alleleFreqs[1,gameteAncestry], 2, 1)
}

simulateGenotypes <- function(M1,M2, rho,x,L) {  
  paternalGamete <- simulateHaploidAlleles(M1,rho,x,L)  
  maternalGamete <- simulateHaploidAlleles(M2,rho,x,L)  
  simulateGenotypes <- paste(paternalGamete, ",", maternalGamete, sep="")
}

##########################################################################
## Start of script
N <- 200
rho <- 6 ## sum-of-intensities
L <- 40 #default # number of loci
##L <- 200 
##marker spacings#
# x <- 0.1 #evenly spaced#
#  x <- runif(L,0.01,0.05) #default, spacings between 1 and 5 cM
x <- rep(100,L) #unlinked markers
# x <- rep(0.01,L) #denser markers

## alleleFreqs <- matrix(c(1, 0, 0, 1), nrow=2, ncol=2)
alleleFreqs <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow=2, ncol=2)
genotypes <- character(L)

## write("outcome\n", file="regression/data/genotypes.txt", append=FALSE)
outcome <- numeric(N)
avM <- numeric(N)
for(individual in 1:N) {
  M1 <- rbeta(1, 4, 1)
  # M2 <- rbeta(1, 4, 1)# random mating
  M2 <- M1 #assortative mating
  avM[individual] <- 0.5*(M1 + M2)
  obs <- simulateGenotypes(M1, M2, rho, x, L)  #make some genotypes missing
  for(locus in 1:L) if(runif(n=1) < 0.1)
    obs[locus]<-""
  genotypes <- rbind(genotypes, obs)
  ## simulate outcome
  outcome[individual] <- rnorm(1, mean=(50 + 25*avM[individual]), sd=3) #linear regression
  ofam <- gaussian
  ##outcome[individual] <- rbinom(1, 1, 1 / (1+exp(-2*avM[individual])))  # binary outcome
  ##ofam <- binomial
}
m2 <- 1 - avM
reg.true <-summary.glm(glm(outcome ~ m2, family = ofam))

outcome.table <- data.frame(outcome, row.names=NULL) # write outcome variable to file
write.table(outcome.table, file="sim/outcome.txt", row.names=FALSE, col.names=TRUE)
Mvector.table <- data.frame(avM, row.names=NULL)
write.table(outcome, file="sim/Mvalues.txt", row.names=FALSE,
            col.names=TRUE)
##write genotypes to file
genotypes <- genotypes[-1,]
id = as.character(seq(1:N))
genotypes <- data.frame(id, genotypes, row.names=NULL)
write.table(genotypes, file="sim/genotypes.txt",
            row.names=FALSE)
## write locus file
loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-1]),  rep(2,L),  x, row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "Distance")
write.table(loci, file="sim/loci.txt", row.names=FALSE)
## write allelefreqs files
allelefreqs <- data.frame(as.vector(loci[,1]), rep(alleleFreqs[1,1],L), rep(alleleFreqs[1,2],L))
dimnames(allelefreqs)[[2]] <- c("Locus", "Pop1", "Pop2")

## old format allelefreqfile
write.table(allelefreqs, file="sim/allelefreqfile.txt", row.names=FALSE)
priorallelefreqs <- numeric(2)
locusnames <- character(0)
for(locus in 1:L) {
  priorallelefreqs <- rbind(priorallelefreqs, 500*alleleFreqs + 0.5)
  locusnames <- c(locusnames, rep(as.vector(loci[locus,1]), 2))
}
## priorallelefreqfile
priorallelefreqs <- data.frame(locusnames, priorallelefreqs[-1,], row.names=NULL)
write.table(priorallelefreqs, file="sim/priorallelefreqs.txt",
            row.names=FALSE)                               

# run analysis with no outcome vars  
system("perl sim.pl")

noreg.quantiles <- read.table(file="sim/PosteriorQuantiles.txt", header=T)
print(noreg.quantiles)
m2.estimates <- read.table("sim/IndividualVarPosteriorMeans.txt", header=T)[, 2]
reg.estimates <- summary.glm(glm(outcome ~ m2.estimates, family=ofam))

plot(m2, m2.estimates, xlim=c(0,1), ylim=c(0,1))
lines(c(0,1),c(0,1), type="l")

system("perl simRegression.pl")

reg.quantiles <- read.table(file="sim/PosteriorQuantiles.txt", header=T)
print(reg.quantiles)
m2.estimates.r <- read.table("sim/IndividualVarPosteriorMeans.txt", header=T)[, 2]
reg.estimates.r <- summary.glm(glm(outcome ~ m2.estimates.r))

plot(m2.estimates, m2.estimates.r, xlim=c(0,1), ylim=c(0,1))
lines(c(0,1),c(0,1), type="l")
