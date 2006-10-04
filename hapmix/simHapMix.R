## script to simulate test data for admixmap
#######################################################################
library("combinat")

rDiscrete <- function(p) {
  s <- 0
  z <- 0.0
  j <- 1
  x <- runif(1)
  while(s==0) {
    z <- z + p[j]
    if(x < z) {
      s <- j
    }
    j <- j + 1
  }
  return(s)
}

simulateAncestry <- function(mu, f, L) {
  ## mu vector of mixture proportions of length K,
  ## rho scalar arrival rate, d vector of distances, L num loci
  Anc <- integer(L)
  for(locus in 1:L) {
    ## xi is indicator variable for > 0 arrivals
    xi <- runif(1) >  f[locus]
    if(xi) {
      Anc[locus] <- rDiscrete(mu)
    } else {
      Anc[locus] <- Anc[locus - 1]
    }
  }
  return(Anc)  
}

simulateHaploidAlleles <- function(mu, f, L, Freqs) {
  Anc <- simulateAncestry(mu, f, L)
  ## Freqs is array S x K x L     
  Alleles <- integer(L)
  for(locus in 1:L) {
    Alleles[locus] <- rDiscrete(Freqs[, Anc[locus], locus])
  }
  return(Alleles)
}

simulateGenotypes <- function(mu1,mu2, f, L, alleleFreqs, allele1.counts) {  
  paternalGamete <- simulateHaploidAlleles(mu1, f, L, alleleFreqs)  
  maternalGamete <- simulateHaploidAlleles(mu2, f, L, alleleFreqs)
  allele1.counts <- (paternalGamete==1) + (maternalGamete==1)
  simulateGenotypes <- paste(paternalGamete, ",", maternalGamete, sep="")
  return(list(genotypes=simulateGenotypes, counts=allele1.counts))
}

distanceFromLast <- function(v.Chr, v.Position) {
  ## given a table of loci with chr number and position in cM or megabases,
  ## this function returns vector of distances from last locus
  ## distance from last coded as NA for first locus on each chromosome 
  OrderedPosition <- 10000*v.Chr + v.Position
  L <- length(OrderedPosition)
  v.DistanceFromLast <- c(10000, (OrderedPosition[-1]) - OrderedPosition[-L])
  v.DistanceFromLast[v.DistanceFromLast > 1000] <- NA
  return(v.DistanceFromLast)
}

##########################################################################
## Start of script
## chromosome lengths in cM
#chr.L <- c(292,272,233,212,197,201,184,166,166,181,156,169,117,128,110,130,128,123,109,96,59,58)
chr.L <- c(20, 20) ## trial runs with 2 chr
#chr.L <- 40##one longer chromosome
numChr <- length(chr.L)

N <- 100##number of individuals
K <- 4##number of block states
rhoalpha <-  40 #    |
rhobeta0 <- 5 ## prior on rho
rhobeta1 <- 5##      |
rhobeta <- rgamma(1, rhobeta0, rhobeta1)
## mean r = rhoalpha*rhobeta1 / (rhobeta0 - 1)
# var =  r (r + 1) / (rhobeta0 - 2) 

spacing <- 0.01 # spacing in cM

## assign map distances
x <- numeric(0)
chr <- integer(0)
for(chromosome in 1:numChr) {
  positions <- seq(0, chr.L[chromosome], spacing)
  x <- c( x, positions) 
  chr <- c(chr, rep(chromosome, length(positions)))
}
distances <- distanceFromLast(chr, x)
L <- length(x) # number of loci

f <- numeric(L)
for(locus in 1:L) {
  if(is.na(distances[locus])) {
    f[locus] <- 0.0
  } else {
    ##rhobeta <- rgamma(1, shape=rhobeta0, rate=rhobeta1)
    rho <- rgamma(1, shape=rhoalpha, rate=rhobeta)
    f[locus] <- exp(-rho*distances[locus])
  }
}

## elements of Dirichlet param vector for prior on allele freqs
mu <- rep(1/K, K)
alleleFreqs <- array(data=NA, dim=c(2, K, L))

freqs.shape <- 3
freqs.rate <- 10
for(locus in 1:L) {
freqs.alpha  <- rgamma(1, shape=freqs.shape, rate=freqs.rate) / K
##  for(state in 1:K) {
##p <- rep(0, K)
##while( (min(p)<(1e-9)) || (max(p)>=(1-(1e-9)))){
   p  <- rbeta(K, freqs.alpha, freqs.alpha) # freqs allele 1
##  }
alleleFreqs[1, , locus] <- p
alleleFreqs[2, , locus] <- 1 - p # freqs allele 2
##  }
}
allele1.counts <- rep(0, L)
genotypes <- matrix(data="0,0", nrow=N, ncol=L)
for(individual in 1:N) {
  g.list <- simulateGenotypes(mu, mu, f, L, alleleFreqs, allele1.counts)
  genotypes[individual, ] <- g.list$genotypes
  allele1.counts <- allele1.counts + g.list$counts
##  genotypes[individual, ] <-simulateHaploidAlleles(mu, f, L, alleleFreqs)  
}

## write genotypes file
id = as.character(seq(1:N))
sex <- rep(1, N)##for all males, irrelevant if no X-chromosome
genotypes <- data.frame(id, sex, genotypes, row.names=NULL)
write.table(genotypes, file="data/genotypes.txt", sep="\t", row.names=FALSE)
## write locus file
distances[is.na(distances)] <- 100

##next 2 lines for X-only data
##loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-c(1:2)]),  rep(2,L),  distances, rep("X", L), row.names=NULL)
##dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "Distanceincm", "Chrm")

##next 2 lines for autosomal-only data
loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-c(1:2)]),  rep(2,L),  distances, row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "Distanceincm")

write.table(loci, file="data/loci.txt", row.names=FALSE)

## write allelefreqsfile
freqstable <- numeric(K)
locusnames <- character(0)
for(locus in 1:L) {
  locusnames <- c(locusnames, rep(paste("X", locus, sep=""), 2))
  freqstable <- rbind(freqstable, alleleFreqs[1:2, , locus])
}
statelabels <- paste("state", seq(1:K), sep="")
freqstable <- data.frame(locusnames, freqstable[-1, ])
dimnames(freqstable)[[2]] <- c("locus", statelabels)
write.table(freqstable, file="data/allelefreqs.txt", sep="\t", row.names=F)

  
