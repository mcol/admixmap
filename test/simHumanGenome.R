## script to simulate test data for admixmap
#######################################################################
simulateHaploidAlleles <- function(M, rho, x, L, alleleFreqs) {
  ## M is proportionate admixture from pop 1 
  gameteAncestry <- numeric(L)  
  simAlleles <- numeric(L)

  randAnc <- runif(L)  
  gameteAncestry[1] <- ifelse(randAnc[1] < M, 1, 2) 
  ## Tmatrix[i,j] is prob ancestry at t+1 = j given ancestry at t = i
  for(locus in 2:L) {
    if(is.na(x[locus])) {
      Tmatrix <- t(matrix(data=c(M, 1-M, M, 1-M), nrow=2, ncol=2))
    } else {
      Tmatrix <- t(matrix(data=c(M + (1-M)*exp(-rho*x[locus]),  1-M - (1-M)*exp(-rho*x[locus]),
                      M - M*exp(-rho*x[locus]), 1-M + M*exp(-rho*x[locus]) ),
                    nrow=2, ncol=2))
    }
    ## simulate ancestry at locus
    gameteAncestry[locus] <- ifelse(randAnc[locus] < Tmatrix[gameteAncestry[locus-1], 1], 1, 2)
  }

  for(locus in 1:L) {
    ## simulate allele conditional on ancestry at locus
simAlleles[locus] <- ifelse(runif(1) < alleleFreqs[2*locus - 1, gameteAncestry[locus]], 1, 2)
  }
  return(simAlleles)
}

simulateGenotypes <- function(M1,M2, rho,x,L, alleleFreqs) {  
  paternalGamete <- simulateHaploidAlleles(M1,rho,x,L, alleleFreqs)  
  maternalGamete <- simulateHaploidAlleles(M2,rho,x,L, alleleFreqs)  
  simulateGenotypes <- paste(paternalGamete, ",", maternalGamete, sep="")
}

distanceFromLast <- function(v.Chr, v.Position) {
## given a table of loci with chr number and position in cM or megabases,
## this function returns vector of distances from last locus
## distance from last coded as NA for first locus on each chromosome 
  OrderedPosition <- 10000*v.Chr + v.Position
  v.DistanceFromLast <- c(10000,
                          (OrderedPosition[-1]) - OrderedPosition[-length(OrderedPosition)])
  v.DistanceFromLast[v.DistanceFromLast > 1000] <- NA
  return(v.DistanceFromLast)
}

##########################################################################
# undebug(simulateHaploidAlleles)
# debug(simulateHaploidAlleles)
## Start of script
N <- 200
K <- 2 # num subpopulations
rho <- 6 # sum-of-intensities
spacing <- 40 # 40 cM spacing gives 99 loci
eta <- 5 # allele freq dispersion parameter
beta <- 2 # regression slope
logistic <- F # logistic or linear

x <- numeric(0)
chr <- numeric(0)
numChr <- 22
## chromosome lengths in cM
chr.L <- c(292,272,233,212,197,201,184,166,166,181,156,169,117,128,110,130,128,123,109,96,59,58)
length <- sum(chr.L)
for(chromosome in 1:22) {
  positions <- seq(0, chr.L[chromosome], spacing)
  x <- c( x, positions) 
  chr <- c(chr, rep(chromosome, length(positions)))
}
x <- 0.01*distanceFromLast(chr, x)
L <- length(x) # number of loci

## simulate correlated allele freqs
mu <- numeric(L) # ancestral freqs allele 1
alleleFreqs <- matrix(data=NA, nrow=2*L, ncol=K)
for(locus in 1:L) {
  mu[locus] <- rbeta(1, 2, 2)
  alleleFreqs[2*locus - 1, ] <- rbeta(2, mu[locus]*eta, (1 - mu[locus])*eta)
                                        # freqs allele 1 in each of K subpops
  alleleFreqs[2*locus, ] <- 1 - alleleFreqs[2*locus - 1, ] # freqs allele 2
}
#alleleFreqs[,1] <- 1
#alleleFreqs[,2] <- 0

## simulate genotypes and outcome variable
genotypes <- character(L)
outcome <- numeric(N)
avM <- numeric(N)
popadmixparams <- c(1, 4) # population admixture params for pop1, pop2
popM <- popadmixparams[1] / sum(popadmixparams)
for(individual in 1:N) {
  M1 <- 1 - rbeta(1, popadmixparams[1], popadmixparams[2]) ## M1 is prob pop 1
  ## M2 <- rbeta(1, 4, 1)# random mating
  M2 <- M1 #assortative mating
  avM[individual] <- 1 - 0.5*(M1 + M2)
  obs <- simulateGenotypes(M1, M2, rho, x, L, alleleFreqs)
  ##make some genotypes missing
  ##for(locus in 1:L) if(runif(n=1) < 0.1)
  ##    obs[locus]<-""
  genotypes <- rbind(genotypes, obs)

  ## simulate outcome
  alpha <- -beta*popM
  if(logistic) { # logistic regression with approx equal numbers of cases and controls
    outcome[individual] <- rbinom(1, 1, 1 / (1+exp(-(alpha + beta*avM[individual]))))  
    ofam <- binomial
  } else { # linear regression
    outcome[individual] <- rnorm(1, mean=(alpha + beta*avM[individual]), sd=1) 
    ofam <- gaussian
  }
}


reg.true <-summary.glm(glm(outcome ~ avM, family = ofam))

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
x[is.na(x)] <- 100
loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-1]),  rep(2,L),  x, row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "Distance")
write.table(loci, file="sim/loci.txt", row.names=FALSE)

  ## write allelefreqs files
trueallelefreqs <- data.frame(rep(loci[,1], each=2), alleleFreqs)
dimnames(trueallelefreqs)[[2]] <- c("Locus", "Pop1", "Pop2")
write.table(trueallelefreqs, file="sim/trueallelefreqs.txt",
            row.names=FALSE)                               
## write reference prior allele freqs
priorallelefreqs <- data.frame(rep(loci[,1], each=2), matrix(data=0.5, nrow=2*L, ncol=K))
dimnames(priorallelefreqs)[[2]] <- c("Locus", "Pop1", "Pop2")
write.table(priorallelefreqs, file="sim/priorallelefreqs.txt",
            row.names=FALSE)                               

# run analysis with outcome var
if(logistic) {
  system("perl simLogistic.pl")
} else system("perl simLinear.pl")

reg.quantiles <- read.table(file="sim/PosteriorQuantiles.txt", header=T)
print(reg.quantiles)
## plot posterior means of individual admixture against true values 
avM.estimates <- read.table("sim/IndividualVarPosteriorMeans.txt", header=T)[, 2]
reg.estimates <- summary.glm(glm(outcome ~ avM.estimates, family=ofam))

plot(avM, avM.estimates, xlim=c(0,1), ylim=c(0,1))
lines(c(0,1),c(0,1), type="l")

