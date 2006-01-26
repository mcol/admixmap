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

simulateGenotypes <- function(M1,M2, rho1, rho2, x,L, alleleFreqs) {  
  paternalGamete <- simulateHaploidAlleles(M1,rho1,x,L, alleleFreqs)  
  maternalGamete <- simulateHaploidAlleles(M2,rho2,x,L, alleleFreqs)  
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
## Start of script
N <- 1
K <- 2 # num subpopulations
rho1 <- 8 # sum-of-intensities
rho2 <- 2
spacing <- 5 # 40 cM spacing gives 99 loci

x <- numeric(0) # distance from last in morgans
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

## simulate allele freqs
alleleFreqs <- matrix(data=NA, nrow=2*L, ncol=K)
pbar <- 0.8
for(locus in 1:L) {
  alleleFreqs[2*locus - 1, 1] <- 0.95 # rbeta(1, 10*pbar, 10*(1-pbar))
  alleleFreqs[2*locus - 1, 2] <- 0.05 # rbeta(1, 10*(1-pbar), 10*pbar)
                                        # freqs allele 1 
  alleleFreqs[2*locus, ] <- 1 - alleleFreqs[2*locus - 1, ] # freqs allele 2
}
p1 <- alleleFreqs[seq(2, 2*L, by=2), 1]
p2 <- alleleFreqs[seq(2, 2*L, by=2), 2]
f <- (p1 - p2)^2 / ((p1+p2)*(2 - p1 - p2))
print(mean(f))

## simulate genotypes
for(individual in 1:N) {
  M1 <- 0.05
  M2 <- 0.95 
  obs <- simulateGenotypes(M1, M2, rho1, rho2, x, L, alleleFreqs)
}

##write genotypes to file
genotypes <- data.frame(matrix(obs, nrow=1))
id = as.character(seq(1:N))
genotypes <- data.frame(id, genotypes, row.names=NULL)
write.table(genotypes, file="data/genotypes.txt", sep="\t",
            row.names=FALSE)

## write locus file
x[is.na(x)] <- 100
loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-1]),  rep(2,L),  x, row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "Distance")
write.table(loci, file="data/loci.txt", sep="\t", row.names=FALSE)

## write allelefreqs files
trueallelefreqs <- data.frame(rep(loci[,1], each=2), alleleFreqs)
dimnames(trueallelefreqs)[[2]] <- c("Locus", "Pop1", "Pop2")
write.table(trueallelefreqs, file="data/trueallelefreqs.txt", sep="\t",
            row.names=FALSE)                               


