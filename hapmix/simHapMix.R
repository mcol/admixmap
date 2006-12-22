## script to simulate test data for admixmap
#######################################################################
#library("combinat")

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
##cat(Anc, " ")
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
#chr.L <-
##c(292,272,233,212,197,201,184,166,166,181,156,169,117,128,110,130,128,123,109,96,59,58)
chr.L <- c(1, 1) ## trial runs with 2 chr
#chr.L <- 20
numChr <- length(chr.L)

N <- 100##number of individuals
K <- 6##number of block states
DiploidData = F

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

##use this to read distances from file
#L<-5000
#K<-4
#N<-60
#distances<-read.table("/ichec/work/ndlif006b/genepi/hapmap/Eur/chr22data/loci5000.txt", header=T, comment.char="", na.strings="#")[,3]


##generate numbers of arrivals and locus correlations
h <- 400 ##average number arrivals per Mb
lambda.rate <- 10##rgamma(1, shape=rhobeta0, rate=rhobeta1)
f <- numeric(L)
for(locus in 1:L) {
  if(is.na(distances[locus])) {
    f[locus] <- 0.0
  } else {
    lambda <- rgamma(1, h * distances[locus], lambda.rate)
    f[locus] <- exp( -lambda )
  }
}

##set block mixture proportions
mu <- rep(1/K, K)

## generate allele freqs
alleleFreqs <- array(data=NA, dim=c(2, K, L))
##use this to read freqs from file
#freqs.alpha<-read.table("/ichec/work/ndlif006b/genepi/hapmap/Eur/Results1States/AlleleFreqPosteriorMeans.txt",
#header=T)[,2]
#freqs.alpha[freqs.alpha==0]<-0.001
#freqs.alpha[freqs.alpha==1]<-0.999

alpha.shape <- 2
alpha.rate <- 10
for(locus in 1:L) {
freqs.alpha <- rgamma(1, shape=alpha.shape, rate=alpha.rate)/K##Gamma with mean 0.1
##p <- rep(0, K)
##while( (min(p)<(1e-9)) || (max(p)>=(1-(1e-9)))){
   #p <- rbeta(K, freqs.alpha[locus*2-1], freqs.alpha[locus*2]) # freqs  allele 1
   p <- rbeta(K, freqs.alpha, freqs.alpha) # freqs allele 1
    #p <- c(0,1)
##  }
alleleFreqs[1, , locus] <- p     # freqs allele 1
alleleFreqs[2, , locus] <- 1 - p # freqs allele 2
##  }
}
allele1.counts <- rep(0, L)
genotypes.diploid <- matrix(data="0,0", nrow=N, ncol=L)
genotypes.haploid <- matrix(data="0", nrow=2*N, ncol=L)
for(individual in 1:N) {

##  g.list <- simulateGenotypes(mu, mu, f, L, alleleFreqs, allele1.counts)
## genotypes.diploid[individual, ] <- g.list$genotypes
## allele1.counts <- allele1.counts + g.list$counts

  paternalGamete <- simulateHaploidAlleles(mu, f, L, alleleFreqs)
  maternalGamete <- simulateHaploidAlleles(mu, f, L, alleleFreqs)
  genotypes.haploid[2*individual-1, ] <- paternalGamete
  genotypes.haploid[2*individual, ] <- maternalGamete
  genotypes.diploid[individual,] <- paste(paternalGamete, ",", maternalGamete, sep="")

}

##write diploid genotypes
id = as.character(seq(1:N))
sex <- rep(1, N)##for all males, irrelevant if no X-chromosome
##genotypes <- data.frame(id, sex, genotypes, row.names=NULL)
genotypes <- data.frame(id, genotypes.diploid, row.names=NULL)
write.table(genotypes, file="data/genotypes.txt", sep="\t", row.names=FALSE)

##write haploid genotypes
id = as.character(seq(1:(2*N)))
sex <- rep(1, 2*N)##for all males, irrelevant if no X-chromosome
##genotypes <- data.frame(id, sex, genotypes, row.names=NULL)
genotypes <- data.frame(id, genotypes.haploid, row.names=NULL)
write.table(genotypes, file="data/genotypes_haploid.txt", sep="\t", row.names=FALSE)

## write locus file
distances[is.na(distances)] <- 100
##chr.names<-rep("X", L)##all X-chromosome

##next 2 lines for X-only data
##loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-c(1:2)]),  rep(2,L),  distances, rep("X", L), row.names=NULL)
##dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "DistanceinMb", "Chrm")

##next 2 lines for autosomal-only data
loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-1]),  rep(2,L),  distances, row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "DistanceinMb")

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

##priorallelefreqs <- freqstable[seq(1, (dim(freqstable)[[1]]), by=2)]
##write.table(freqstable+0.5, file="data/priorallelefreqs.txt", sep="\t", row.names=F)

  
