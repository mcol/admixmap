## script to simulate test data for admixmap
#######################################################################
#library("combinat")

############################################
## Describe Model here                     #
############################################

##where to write data files
datadir <- "hapmixsimdata"

N <- 60 ##number of haploid individuals
K <- 8 ##number of block states
NumCases <- 100
NumControls<-200

## chromosome lengths in Mb

chr.L <- 0.5        ##single chromosome
#chr.L <- c(10, 10) ## trial runs with 2 chr

numChr <- length(chr.L)##number of chromosomes

spacing <- 0.002 # spacing in Mb

h <- 400         ##average number arrivals per unit distance
lambda.rate <- 10##rgamma(1, shape=rhobeta0, rate=rhobeta1)

extreme.allele.freqs <- F##set to True for allelefreqs of 0 and 1

freq.dispersion.prior.shape <- 1
freq.dispersion.prior.rate <- 10

mixture.proportions.Dirichlet.Prior.params <- rep(1, K)
fixed.mixture.proportions <- F##set to True for mixture proportions fixed at 1/K
                              ##set to F to sample mixture proportions for each locus  
########################################################################
####################################
## functions required for script   #
####################################

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
  ## mu matrix of mixture proportions of dimension L*K,
  ## rho scalar arrival rate, d vector of distances, L num loci
  Anc <- integer(L)
  for(locus in 1:L) {
    ## xi is indicator variable for > 0 arrivals
    xi <- runif(1) >  f[locus]
    if(xi) {
      Anc[locus] <- rDiscrete(mu[locus,])
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

rdirichlet <- function(params) {
  S <- length(params)
  p <- numeric(S)
  if(S==2) {
    p[1] <- rbeta(1, params[1], params[2])
    p[2] <- 1 - p[1]
  } else {
    for(i in 1:S) {
      p[i] <- rgamma(1, params[i], 1)
    }
    sum.p <- sum(p)
    p <- p / sum.p
  }
  return(p)
}
##############################
## Start of script          ##
##############################

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

##generate numbers of arrivals and locus correlations
f <- numeric(L)
for(locus in 1:L) {
  if(is.na(distances[locus])) {
    f[locus] <- 0.0
  } else {
    lambda <- rgamma(1, h * distances[locus], lambda.rate)
    f[locus] <- exp( -lambda )
  }
}

##set default block mixture proportions
mixture.proportions <- matrix(1/K, nrow=L, ncol=K)

##create matrix of allele freqs
alleleFreqs <- array(data=NA, dim=c(2, K, L))

for(locus in 1:L) {
## generate allele freqs
  if(extreme.allele.freqs){
     freq.allele1 <- c(0,1)
  }else{
   freqs.dispersion <- rgamma(1, shape=freq.dispersion.prior.shape, rate=freq.dispersion.prior.rate)
   freq.proportion.allele1 <- runif(1)
   freq.proportion.allele2 <- 1-freq.proportion.allele1
   freqs.beta.parameter1 <- freqs.dispersion * freq.proportion.allele1
   freqs.beta.parameter2 <- freqs.dispersion * freq.proportion.allele2

    ##while( (min(p)<(1e-9)) || (max(p)>=(1-(1e-9)))){

     freq.allele1 <- rbeta(K, freqs.beta.parameter.1, freqs.beta.parameter2) # freqs allele 1
     freq.allele1 <- c(0,1)
    ##  }
  }
  alleleFreqs[1, , locus] <- freq.allele1     # freqs allele 1
  alleleFreqs[2, , locus] <- 1 - freq.allele1 # freqs allele 2

  ##simulate mixture proportions
  if(!fixed.mixture.proportions)
    mixture.proportions[locus,] <- rdirichlet(mixture.proportions.Dirichlet.Prior.params)

}##end of locus loop

allele1.counts <- rep(0, L)
genotypes.diploid <- matrix(data="0,0", nrow=N, ncol=L)
genotypes.haploid <- matrix(data="0", nrow=2*N, ncol=L)
for(individual in 1:N) {

##  g.list <- simulateGenotypes(mu, mu, f, L, alleleFreqs, allele1.counts)
## genotypes.diploid[individual, ] <- g.list$genotypes
## allele1.counts <- allele1.counts + g.list$counts

  paternalGamete <- simulateHaploidAlleles(mixture.proportions, f, L, alleleFreqs)
  maternalGamete <- simulateHaploidAlleles(mixture.proportions, f, L, alleleFreqs)
  genotypes.haploid[2*individual-1, ] <- paternalGamete
  genotypes.haploid[2*individual, ] <- maternalGamete
  genotypes.diploid[individual,] <- paste(paternalGamete, ",", maternalGamete, sep="")

}

##write diploid genotypes
id = as.character(seq(1:N))
sex <- rep(1, N)##for all males, irrelevant if no X-chromosome
#genotypes <- data.frame(id, sex, genotypes, row.names=NULL)
genotypes <- data.frame(id, genotypes.diploid)
dimnames(genotypes)[[2]] <- c("ID", paste("X", 1:L, sep=""))
write.table(genotypes, file=paste(datadir,"genotypes_diploid.txt", sep="/"), sep="\t", row.names=F, col.names=T)

##write haploid genotypes
id = as.character(seq(1:(2*N)))
sex <- rep(1, 2*N)##for all males, irrelevant if no X-chromosome
##genotypes <- data.frame(id, sex, genotypes, row.names=NULL)
genotypes <- data.frame(id, genotypes.haploid)
write.table(genotypes, file=paste(datadir,"genotypes_haploid.txt", sep="/"), sep="\t", row.names=F, col.names=T)

## write locus file
distances[is.na(distances)] <- 100
##chr.names<-rep("X", L)##all X-chromosome

##next 2 lines for X-only data
##loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-c(1:2)]),  rep(2,L),  distances, rep("X", L), row.names=NULL)
##dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "DistanceinMb", "Chrm")

##next 2 lines for autosomal-only data
loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-1]),  rep(2,L),  distances, row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "DistanceinMb")

write.table(loci, file=paste(datadir,"loci.txt", sep="/"), row.names=FALSE)

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
write.table(freqstable, file=paste(datadir,"allelefreqs.txt", sep="/"), sep="\t", row.names=F)

##priorallelefreqs <- freqstable[seq(1, (dim(freqstable)[[1]]), by=2)]
##write.table(freqstable+0.5, file="data/priorallelefreqs.txt", sep="\t", row.names=F)

##write mixture proportions to file
mixture.props <- data.frame(mu)
mixture.props <- rbind(mixture.props, colMeans(mixture.props))
dimnames(mixture.props)[[1]][L+1] <- "average"
write.table(round(mixture.props, 4), file=paste(datadir, "TrueMixtureProps.txt", sep="/"), row.names=T, col.names=T)


##generate a case-control genotypes file  
##done exactly as for diploid data but with fewer individuals and outputting only some of the loci
NN <- NumCases + NumControls

CCLoci <- seq(from=2, to=L, by=2)#even-numbered loci
genotypes.cc <- matrix(data="0,0", nrow=NN, ncol=L)
for(individual in 1:NN) {

##  g.list <- simulateGenotypes(mu, mu, f, L, alleleFreqs, allele1.counts)
## genotypes.diploid[individual, ] <- g.list$genotypes
## allele1.counts <- allele1.counts + g.list$counts

  paternalGamete <- simulateHaploidAlleles(mu, f, L, alleleFreqs)
  maternalGamete <- simulateHaploidAlleles(mu, f, L, alleleFreqs)
  genotypes.cc[individual,] <- paste(paternalGamete, ",", maternalGamete, sep="")

}
id = as.character(seq(1:NN))

##write even-numbered genotypes
genotypes <- data.frame(id, genotypes.cc[,seq(from=2, to=L, by=2)])
##write even-numbered genotypes to file
dimnames(genotypes)[[2]] <- c("ID", paste("X", CCLoci, sep=""))
write.table(genotypes, file=paste(datadir,"genotypes_casectrl.txt", sep="/"), sep="\t", row.names=F, col.names=T)

##write an outcome variable file
##outcome <- data.frame(Outcome = sample(size=NN, x=c(0,1), replace=T))
outcome <- data.frame(Outcome = c(rep(0, NumControls), rep(1, NumCases)))
write.table(outcome, file=paste(datadir,"outcome.txt", sep="/"), row.names=F, col.names=T)
