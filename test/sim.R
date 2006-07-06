## script to simulate test data for admixmap
simulateHaploidAlleles <- function(M,rho,x,L) {   
  gameteAncestry <- numeric(L)  
  randAnc <- runif(L)
  f <- numeric(L)
  gameteAncestry[1] <- ifelse(randAnc[1] > M, 2, 1) # M is prob pop 1.
  if (L > 1) {
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
  }
  simulateHaploidAlleles <- ifelse(runif(L) > alleleFreqs[1,gameteAncestry], 2, 1)
}

simulateAutosomalGenotypes <- function(M1,M2, rho,x,L) {
  maternalGamete <- simulateHaploidAlleles(M2,rho,x,L)  
  paternalGamete <- simulateHaploidAlleles(M1,rho,x,L)
  g <- paste(paternalGamete, ",", maternalGamete, sep="")
  return(g)
}

simulateXGenotypes <- function(M1,M2, rho,x,L, male) {
  maternalGamete <- simulateHaploidAlleles(M2,rho,x,L)  
  if(!male) {
    paternalGamete <- simulateHaploidAlleles(M1,rho,x,L)
    g <- paste(paternalGamete, ",", maternalGamete, sep="")
  } else {
    g <- paste(maternalGamete, ",", maternalGamete, sep="") # 2nd element should be 0
  }
  return(g)
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
K <- 2
N <- 200
rho <- 6 ## sum-of-intensities
rhoX <- 0.5*rho
spacing <- 40 # 40 cM spacing gives 99 loci
L <- 50 #default # number of loci
beta <- 2 # regression slope
popadmixparams <- c(3, 1) # population admixture params for pop1, pop2
afreqparams <- c(10, 90)

##marker spacings#
# x <- 0.1 #evenly spaced#
#x <- runif(L,0.01,0.05) #default, spacings between 1 and 5 cM
#x <- rep(NA,L) #unlinked markers
# x <- rep(0.01,L) #denser markers

x <- numeric(0)
chr <- numeric(0)
numChr <- 22 # not including X chr
## chromosome lengths in cM
chr.L <- c(292,272,233,212,197,201,184,166,166,181,156,169,117,128,110,130,128,123,109,96,59,58)
length <- sum(chr.L)
for(chromosome in 1:22) {
  positions <- seq(0, chr.L[chromosome], spacing)
  x <- c( x, positions) 
  chr <- c(chr, rep(chromosome, length(positions)))
}
x <- 0.01*distanceFromLast(chr, x)
L <- length(x) # number of autosomal loci

positionsX <- as.vector(seq(0, 188, 200))#spacing))
chrX <- as.vector(rep(1, length(positionsX)))
xX <- as.vector(0.01*distanceFromLast(chrX, positionsX))
LX <- length(xX) # number of X loci

## simulate allele freqs
alleleFreqs <- matrix(data=NA, nrow=2*(L+LX), ncol=K)
for(locus in 1:(L + LX)) {
  alleleFreqs[2*locus - 1, 1] <- rbeta(1, afreqparams[1], afreqparams[2]) # freqs allele 1 in subpop 1
  alleleFreqs[2*locus - 1, 2] <- rbeta(1, afreqparams[2], afreqparams[1]) # freqs allele 1 in subpop 2
  alleleFreqs[2*locus, ] <- 1 - alleleFreqs[2*locus - 1, ] # freqs allele 2
}

genotypes <- character(L+LX)
outcome <- numeric(N)
avM <- numeric(N)
male <- seq(0, N) #rbinom(N, 1, 0.5)
popM <- popadmixparams[2] / sum(popadmixparams) # mean admixture proportions
for(individual in 1:N) {
  M1 <- rbeta(1, popadmixparams[1], popadmixparams[2])
  # M2 <- rbeta(1, popadmixparams[1], popadmixparams[2])# random mating
  M2 <- M1 #assortative mating
  avM[individual] <- 1 - 0.5*(M1 + M2)
  obs <- simulateAutosomalGenotypes(M1, M2, rho, x, L)  #make some genotypes missing
  for(locus in 1:L) if(runif(n=1) < 0.1)
    obs[locus] <- "0,0"
  obs <- c(obs, simulateXGenotypes(M1, M2, rhoX, xX, LX, as.logical(male[individual])))
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
genotypes <- data.frame(id, 2-male, genotypes, row.names=NULL)
write.table(genotypes, file="simdata/genotypes.txt", row.names=FALSE, sep="\t")

## write locus file
chrlabels <- c(as.character(chr), rep("X", LX))
loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-(1:2)]),  rep(2,L+LX), c(x, xX), chrlabels, 
                   row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "Distance", "Chr")
write.table(loci, file="simdata/loci.txt", row.names=FALSE, sep="\t")

## write priorallelefreqs file
priorallelefreqs <- t(matrix(c(afreqparams[1], afreqparams[2],
                               afreqparams[2], afreqparams[1]), ncol=2*(L+LX), nrow=2))
locusnames <- rep(loci[, 1], each=2)
priorallelefreqs <- data.frame(locusnames, priorallelefreqs)
dimnames(priorallelefreqs)[[2]] <- c("Locus", "Pop1", "Pop2")
write.table(priorallelefreqs, file="simdata/priorallelefreqs.txt", sep="\t", row.names=FALSE)
