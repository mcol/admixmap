library(R.utils)

## script to simulate test data for admixmap
#######################################################################

simulateMeiosis <- function(x, T) {
  ## returns vector of segregation indicators at loci 1 to T in a single meiosis
  seg <- integer(T) # takes values 0 or 1
  randmei <- runif(T)  
  f <- exp(-x) # prob 0 arrivals in preceding interval
  seg[1] <- ifelse(randmei[1] < 0.5, 0, 1)
  for(locus in 2:T) { # calculate transition matrix
    if(is.na(x[locus])) { ## distance from last missing
      Tmatrix <- t(matrix(data=c(0.5, 0.5, 0.5, 0.5), nrow=2, ncol=2))
    } else {
      Tmatrix <- 0.5 * t(matrix(data=c(1 + f[locus],  1 - f[locus],
                                  1 - f[locus], 1 + f[locus] ),
                                nrow=2, ncol=2))
    }
    ## simulate crossovers
    seg[locus] <- ifelse(randmei[locus] < Tmatrix[1+seg[locus-1], 1], 0, 1)
  }
  return(seg)
}

simulateGameteFromParent <- function(parent.genotypes, x, T) {
  ## returns vector of simulated alleles on one gamete
  t.alleles <- integer(T)
  parent.alleles <- matrix(unlist(strsplit(parent.genotypes, ",")), nrow=2)
  seg <- simulateMeiosis(x, T)
  for(locus in 1:T) {
    t.alleles[locus] <- parent.alleles[1+seg[locus], locus]
  }
  return(t.alleles)
}

simulateOffspringFromParents <- function(sex, parent1.genotypes, parent2.genotypes, x, L, Xchr.L) {
  ## returns vector of simulated offspring genotypes given parental genotypes 
  offspring.gamete1 <- simulateGameteFromParent(parent1.genotypes, x, L+Xchr.L) #gamete from father
  offspring.gamete2 <- simulateGameteFromParent(parent2.genotypes, x, L+Xchr.L) #gamete from mother
  if(Xchr.L > 0 & sex==1) { #male - code X chr genotypes as homozygous for maternal allele
    offspring.gamete1[(L+1):(L+Xchr.L)] <-  offspring.gamete1[(L+1):(L+Xchr.L)] # 
  }
  offspring.genotypes <- paste(offspring.gamete1, offspring.gamete2, sep=",")
  return(offspring.genotypes)
}

simulateHaploidAlleles <- function(M, rho, x, L, alleleFreqs) {
  ## returns vector of simulated alleles on one admixed gamete
  ## M is proportionate admixture from pop 1 
  gameteAncestry <- integer(L) # takes values 1 or 2
  simAlleles <- integer(L)
  randAnc <- runif(L)  
  gameteAncestry[1] <- ifelse(randAnc[1] < M, 1, 2)
  f <- exp(-rho * 0.01 * x) # prob 0 arrivals in preceding interval
  ## Tmatrix[i,j] is prob ancestry at t+1 = j given ancestry at t = i
  for(locus in 2:L) {
    if(is.na(x[locus])) { ## distance from last missing
      Tmatrix <- t(matrix(data=c(M, 1-M, M, 1-M), nrow=2, ncol=2))
    } else {
      Tmatrix <- t(matrix(data=c(M + (1-M)*f[locus],  1-M - (1-M)*f[locus],
                            M - M*f[locus], 1-M + M*f[locus] ),
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

simulateGenotypes <- function(sex, M1, M2, rho, x, L, Xchr.L, alleleFreqs) {
  # returns a vector of genotypes for a single diploid individual
  maternalGamete <- simulateHaploidAlleles(M2,rho,x,L+Xchr.L, alleleFreqs)  
  paternalGamete <- simulateHaploidAlleles(M1,rho, x, L+Xchr.L, alleleFreqs)
  diploidAlleles <- paste(paternalGamete, maternalGamete, sep=",")
  if(Xchr.L > 0 & sex==1) { ## haploid at X chr loci - code as homozygous for maternal allele
    diploidAlleles[(L+1):(L+Xchr.L)] <- gsub("([12]),([12])", "\\2,\\2", diploidAlleles[(L+1):(L+Xchr.L)])
    ## or code as haploid with maternal gamete only
    ## diploidAlleles[(L+1):(L+Xchr.L)] <- gsub("([12]),([12])", "\\2", diploidAlleles[(L+1):(L+Xchr.L)])
  }
  return(diploidAlleles)
}

simulateSibPair <- function(sex2, popadmixparams, rho, dist, L, Xchr.L, alleleFreqs) {
  ## simulate sib-pair
  ## simulate founder gametes (parental genotypes)
  ind <- simulateIndividual(1, popadmixparams, rho, dist, L, Xchr.L, alleleFreqs)
  father.genotypes  <- ind$genotypes
  father.avM <- ind$avM
  ind <- simulateIndividual(2, popadmixparams, rho, dist, L, Xchr.L, alleleFreqs)
  mother.genotypes  <- ind$genotypes
  mother.avM <- ind$avM
  ## simulate offspring gametes from founder gametes
  sib1.genotypes <- simulateOffspringFromParents(sex2[1], father.genotypes, mother.genotypes, x, L, Xchr.L)
  sib2.genotypes <- simulateOffspringFromParents(sex2[2], father.genotypes, mother.genotypes, x, L, Xchr.L)
  return(list(sib1=sib1.genotypes, sib2=sib2.genotypes))
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

simulateIndividual <- function(sex, popadmixparams, rho, dist, L, Xchr.L, alleleFreqs) {
  M1 <- 1 - rbeta(1, popadmixparams[1], popadmixparams[2]) ## M1 is prob pop 1
  M2 <- M1 #assortative mating
  avM <- 1 - 0.5*(M1 + M2)
  genotypes <- simulateGenotypes(sex, M1, M2, rho, dist, L, Xchr.L, alleleFreqs)
  ##make some genotypes missing
  ##for(locus in 1:L) if(runif(n=1) < 0.1)
  ##    obs[locus]<-"0,0"
  return(list(genotypes=genotypes, avM=avM))  
}

##########################################################################
## Start of script
## specify genome and marker panel
numChr <- 22
## chromosome lengths in cM
chr.L <- c(292,272,233,212,197,201,184,166,166,181,156,169,117,128,110,130,128,
           123,109,96,59,58,120)  # last value is length of X chr
## assign map distances
x <- numeric(0)
chr <- numeric(0)
length <- sum(chr.L)
chr.labels <- c(as.character(1:22), "X")
spacing <- 40 # 40 cM spacing gives 99 autosomal loci
for(chromosome in 1:23) {
  positions <- seq(0, chr.L[chromosome], spacing)
  x <- c( x, positions) 
  chr <- c(chr, rep(chr.labels[chromosome], length(positions)))
}
chrnum <- chr
chrnum[chr=="X"] <- "23"
dist <- distanceFromLast(as.integer(chrnum), x)
Xchr.L <- length(positions) # number of X chr loci
L <- length(x) - Xchr.L # number of autosomal loci

rho <- 6 # sum-of-intensities
eta <- 5 # allele freq dispersion parameter #10 is upper limit with 200 obs and admixmparams Di(1,2)
NumSubPops <- 2 # num subpopulations
popadmixparams <- c(2, 6) # population admixture params for pop1, pop2

## simulate allelefreqs
mu <- numeric(L+Xchr.L) # ancestral freqs allele 1
alleleFreqs <- matrix(data=NA, nrow=2*(L+Xchr.L), ncol=NumSubPops)
for(locus in 1:(L+Xchr.L)) {
  mu[locus] <- rbeta(1, 2, 2)
  alleleFreqs[2*locus - 1, ] <- rbeta(2, mu[locus]*eta, (1 - mu[locus])*eta)
                                        # freqs allele 1 in each of NumSubPops subpops
  alleleFreqs[2*locus, ] <- 1 - alleleFreqs[2*locus - 1, ] # freqs allele 2
}

## set allele freqs as required for testing purposes
## use 0, 1, and 1, 0 to make all markers informative
## use 0.5, 0.5 and 0.5, 0.5 to make all markers uninformative

alleleFreqs[,1] <- c(1, 0)
alleleFreqs[,2] <- c(0, 1)

##############################################################
popM <- popadmixparams[1] / sum(popadmixparams) # mean admixture proportions
beta <- 2 # regression slope for effect of admixture
alpha <- -beta*popM 
logistic <- TRUE # logistic or linear

N.ind <- 100
N.sibpairs <- 0

####################################################################
N <- N.ind + 4*N.sibpairs
## first 6 cols of pedfile: 
ped6.ind <- matrix(data=0, ncol=6, nrow=N.ind)
ped6.ind[, 1] <- seq(1:N.ind)
ped6.ind[, 2] <- seq(1:N.ind)

ped6.sibpair <-  matrix(data=0, ncol=6, nrow=4*N.sibpairs)
ped6.sibpair[, 1] <- N.ind + rep(1:N.sibpairs, each=4) # 4 members of each pedigree
ped6.sibpair[, 2] <- seq((N.ind+1):N)

## simulate unrelated individuals
## all female
sex.ind <- rep(2, N.ind) # 2 - rbinom(N.ind, 1,  0.5)
## simulate genotypes 
genotypes.ind <- matrix(data="0,0", nrow=N.ind, ncol=L+Xchr.L)
outcome.ind <- integer(N.ind)
avM <- numeric(N.ind)
for(i in 1:N.ind) {
  ind <- simulateIndividual(sex.ind[i], popadmixparams, rho, dist, L, Xchr.L, alleleFreqs)
  genotypes.ind[i, ]  <- ind$genotypes
  avM[i] <- ind$avM
  
  ## simulate outcome
  if(logistic) { # logistic regression with approx equal numbers of cases and controls
    outcome.ind[i] <- rbinom(1, 1, 1 / (1+exp(-(alpha + beta*avM[i] ))))  
  } else { # linear regression
    outcome.ind[i] <- rnorm(1, mean=(alpha + beta*avM[i]), sd=1) 
  }
}
ped6.ind[1:N.ind, 5] <- sex.ind
ped6.ind[1:N.ind, 6] <- 1 + outcome.ind

if(N.sibpairs > 0) {
  ## simulate sibpairs
  genotypes.sibpair <- matrix(data="0,0", nrow=4*N.sibpairs, ncol=L+Xchr.L)
  for(i in 1:N.sibpairs) {
    sex2 <- 1 + rbinom(2, 1, 0.5)
    sibpair <- simulateSibPair(sex2, popadmixparams, rho, dist, L, Xchr.L, alleleFreqs)
    genotypes.sibpair[4*i - 1, ] <- sibpair$sib1
    genotypes.sibpair[4*i, ] <- sibpair$sib2
    ped6.sibpair[4*i-1, 3] <- ped6.sibpair[4*i-3, 2] 
    ped6.sibpair[4*i-1, 4] <- ped6.sibpair[4*i-2, 2]
    ped6.sibpair[4*i, 3] <- ped6.sibpair[4*i-3, 2] 
    ped6.sibpair[4*i, 4] <- ped6.sibpair[4*i-2, 2]
    ped6.sibpair[(4*i-3):(4*i), 5] <- c(1, 2, sex2)  
  }
  
  genotypes <- data.frame(rbind(genotypes.ind, genotypes.sibpair))
} else {
  genotypes <- data.frame(genotypes.ind)
}

locusnames <- paste("c", chr, seq(1:(L+Xchr.L)), sep=".")
dimnames(genotypes)[[2]] <- locusnames
         
ped6.sibpair[, 6] <- 2 ## affected sib-pairs

ped6 <- rbind(ped6.ind, ped6.sibpair)
colnames(ped6) <- c("famid", "individ", "patid", "matid", "sex", "outcome")

ped <- data.frame(ped6, genotypes)
  
###############################################################################
# write simulated data to files 
mkdirs("data") # returns FALSE if directory already exists

## write pedfile format: 
write.table(ped, file="data/genotypes.ped", sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)
## write in standard ADMIXMAP format for unrelated individuals
write.table(ped[, -c(1, 3:4, 6)], file="data/genotypes.txt", quote=FALSE,
            row.names=FALSE, col.names=TRUE)
            
## write outcome variable to file
outcome.table <- data.frame(ped6[, 6] - 1, row.names=NULL) 
write.table(outcome.table, file="data/outcome.txt", row.names=FALSE, col.names="outcome")
## write true admixture proportions to file
Mvector.table <- data.frame(avM, row.names=NULL)
write.table(Mvector.table, file="data/Mvalues.txt", row.names=FALSE,
            col.names=TRUE)

## write locus file
x[is.na(x)] <- NA
loci <- data.frame(locusnames, rep(2, L + Xchr.L),  dist, chr, row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "cM", "chr")
write.table(loci, file="data/loci.txt", row.names=FALSE, quote=FALSE)

## write allelefreqs files
trueallelefreqs <- data.frame(rep(loci[,1], each=2), alleleFreqs)
dimnames(trueallelefreqs)[[2]] <- c("Locus", "Pop1", "Pop2")
write.table(trueallelefreqs, file="data/trueallelefreqs.txt",
            row.names=FALSE, quote=FALSE)                               

## draw allele counts given true values
alleleCounts <- array(data=NA, dim=c(L+Xchr.L, 2, 2)) # 2nd dim indexes subpops, 3rd dim indexes alleles
for(locus in 1:(L+Xchr.L)) {
  alleleCounts[locus,1,1] <- rbinom(1, 2*N, alleleFreqs[2*locus-1, 1])
  alleleCounts[locus,2,1] <- rbinom(1, 2*N, alleleFreqs[2*locus-1, 2])
  alleleCounts[locus,,2] <- 2*N - alleleCounts[locus,,1]
}

priorallelefreqs <- as.data.frame(0.5+array(as.vector(aperm(alleleCounts, c(3,1,2))),
                                            dim=c((L+Xchr.L)*2,2)))
priorallelefreqs <- data.frame(locusnames, priorallelefreqs)
write.table(priorallelefreqs, file="data/priorallelefreqs.txt", row.names=FALSE, sep="\t",
            quote=FALSE)


  
