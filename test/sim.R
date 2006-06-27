## script to simulate test data for admixmap
simulateHaploidAlleles <- function(M,rho,x,L) {   
  gameteAncestry <- numeric(L)  
  randAnc <- runif(L)
  f <- numeric(L)
  gameteAncestry[1] <- ifelse(randAnc[1] > M, 2, 1) # M is prob pop 1.
  for(locus in 2:L) {
    if(!is.na(x[locus])) {
      f[locus] <- exp(-rho*x[locus])
    } else {
      f[locus] <- 0
    }
    T <- t(matrix(data=c(M + (1-M)*f[locus],  1-M - (1-M)*f[locus],
                    M - M*f[locus], 1-M + M*f[locus] ),
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
K <- 2
N <- 500
rho <- 6 ## sum-of-intensities
L <- 50 #default # number of loci
beta <- 2 # regression slope

##marker spacings#
# x <- 0.1 #evenly spaced#
#x <- runif(L,0.01,0.05) #default, spacings between 1 and 5 cM
x <- rep(NA,L) #unlinked markers
# x <- rep(0.01,L) #denser markers

## simulate allele freqs
alleleFreqs <- matrix(data=NA, nrow=2*L, ncol=K)
for(locus in 1:L) {
  alleleFreqs[2*locus - 1, 1] <- rbeta(1, 20, 80) # freqs allele 1 in subpop 1
  alleleFreqs[2*locus - 1, 2] <- rbeta(1, 80, 20) # freqs allele 1 in subpop 2
  alleleFreqs[2*locus, ] <- 1 - alleleFreqs[2*locus - 1, ] # freqs allele 2
}

genotypes <- character(L)
outcome <- numeric(N)
avM <- numeric(N)
popadmixparams <- c(3, 1) # population admixture params for pop1, pop2
popM <- popadmixparams[2] / sum(popadmixparams) # mean admixture proportions
for(individual in 1:N) {
  M1 <- rbeta(1, popadmixparams[1], popadmixparams[2])
  # M2 <- rbeta(1, popadmixparams[1], popadmixparams[2])# random mating
  M2 <- M1 #assortative mating
  avM[individual] <- 1 - 0.5*(M1 + M2)
  obs <- simulateGenotypes(M1, M2, rho, x, L)  #make some genotypes missing
  for(locus in 1:L) if(runif(n=1) < 0.1)
    obs[locus]<-""
  genotypes <- rbind(genotypes, obs)
  ## simulate outcome
  alpha <- -beta*popM
  #outcome[individual] <- rnorm(1, mean=(alpha + beta*avM[individual]), sd=1) #linear regression
  #ofam <- gaussian
  outcome[individual] <- rbinom(1, 1, 1 / (1+exp(-alpha - beta*avM[individual])))  # binary outcome
  ofam <- binomial
}
reg.true <-summary.glm(glm(outcome ~ avM, family = ofam))

outcome.table <- data.frame(outcome, row.names=NULL) # write outcome variable to file
write.table(outcome.table, file="simdata/outcome.txt", row.names=FALSE, sep="\t")
Mvector.table <- data.frame(avM, row.names=NULL)
write.table(outcome, file="simdata/Mvalues.txt", row.names=FALSE, sep="\t")
##write genotypes file
genotypes <- genotypes[-1, ]
id = as.character(seq(1:N))
genotypes <- data.frame(id, genotypes, row.names=NULL)
write.table(genotypes, file="simdata/genotypes.txt", row.names=FALSE, sep="\t")

## write locus file
loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-1]),  rep(2,L),  x, row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "Distance")
write.table(loci, file="simdata/loci.txt", row.names=FALSE, sep="\t")

## write priorallelefreqs file
priorallelefreqs <- t(matrix(c(20,80,80,20), ncol=2*L, nrow=2))
locusnames <- rep(loci[, 1], each=2)
priorallelefreqs <- data.frame(locusnames, priorallelefreqs)
dimnames(priorallelefreqs)[[2]] <- c("Locus", "Pop1", "Pop2")
write.table(priorallelefreqs, file="simdata/priorallelefreqs.txt", sep="\t", row.names=FALSE)
