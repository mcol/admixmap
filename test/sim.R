## script to simulate test data for admixmap

simulateCox <- function(N, TT, beta, x) {
  r <- matrix(data = 0, nrow = N, ncol = TT)
  n <- matrix(data = TRUE, nrow = N, ncol = TT) # all individuals enter at time 0
  entryt <- rep(0, N)
  tfail <- numeric(N)
  fail <- rep(1, N)
  d <- numeric(TT)
  exp.xbeta <- numeric(N)
  for(i in 1:N) {
    exp.xbeta[i] <- exp(x[i]*beta)
  }
  ## sample vector of TT baseline hazard rates from gamma distribution
  alpha <- rgamma(TT, 1, 1)
  ## sample failure times in t th interval from exponential distribution
  sumd <- 0
  for(t in 1:TT) {
    tmin <- 10000 
    for(i in 1:N) {
      if(n[i, t]) {
        lambda <- alpha[t] * exp.xbeta[i]
        f <- rexp(1, lambda)
        if(f < tmin) {
          tmin <- f
          j <- i  # indexes individual with min fail time
        }
      }
    }
    d[t] <- tmin
    sumd <- sumd + d[t]
    n[j, t:TT] <- FALSE
    tfail[j] <- sumd
  }
  outcome <- data.frame(entryt, tfail, fail)
  return(outcome)
}  
  
rdiscrete <- function(probs) {
  v <- as.vector(rmultinom(1, 1, probs))
  return(match(1, v))
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

  
simulateHaploidAlleles <- function(M,rho,x,L, freqs, S) {   
  gameteAncestry <- integer(L)  
  f <- numeric(L)
  alleles <- integer(L)
  randAnc <- runif(L)
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
  for(locus in 1:L) {
    locus.rows <- seq(from=1+S*(locus-1), to=S*locus)
    probs <- freqs[locus.rows, gameteAncestry[locus]]
    alleles[locus] <- rdiscrete(probs)
  }
  return(alleles)
}

simulateAutosomalGenotypes <- function(M1,M2, rho,x,L, freqs, S) {
  maternalGamete <- simulateHaploidAlleles(M2,rho,x,L, freqs, S)  
  paternalGamete <- simulateHaploidAlleles(M1,rho,x,L, freqs, S)
  g <- paste(paternalGamete, ",", maternalGamete, sep="")
  return(g)
}

simulateXGenotypes <- function(M1,M2, rho,x,LX, freqsX, S, male) {
  maternalGamete <- simulateHaploidAlleles(M2,rho,x,LX, freqsX, S)  
  if(!male) {
    paternalGamete <- simulateHaploidAlleles(M1,rho,x,LX, freqsX, S)
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
spacingX <- 30 #30
L <- 50  # default number of autosomal loci
beta <- 2 # regression slope
model <- 3 # 1 linear, 2 logistic, 3 cox
popadmixparams <- c(3, 1) # population admixture params for pop1, pop2
S <- 2 # number of alleles
afreqparams <- matrix(c(10, 90, 90, 10), nrow=2, ncol=2)
#afreqparams <- matrix(c(10, 10, 80, 80, 10, 10), nrow=3, ncol=2)

##marker spacings#
# x <- 0.1 #evenly spaced#
#x <- runif(L,0.01,0.05) #default, spacings between 1 and 5 cM
#x <- rep(NA,L) #unlinked markers
# x <- rep(0.01,L) #denser markers

x <- numeric(0)
chr <- integer(0)
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

positionsX <- as.vector(seq(0, 188, by=spacingX))
chrX <- as.vector(rep(1, length(positionsX)))
xX <- as.vector(0.01*distanceFromLast(chrX, positionsX))
LX <- length(xX) # number of X loci

## simulate allele freqs
alleleFreqs <- matrix(data=NA, nrow=S*(L+LX), ncol=K)
for(locus in 1:(L + LX)) {
  locus.rows <- seq(from=1+S*(locus-1), to=S*locus)
  for(pop in 1:K) {
    alleleFreqs[locus.rows, pop] <- afreqparams[, pop]/sum(afreqparams[, pop])
    ##  alleleFreqs[locus.rows, pop] <- rdirichlet(afreqparams[, pop])
  }
}                                             

genotypes <- character(L+LX)
avM <- numeric(N)
male <- rbinom(N, 1, 0.5)
popM <- popadmixparams[2] / sum(popadmixparams) # mean admixture proportions
for(individual in 1:N) {
  M1 <- rbeta(1, popadmixparams[1], popadmixparams[2])
  # M2 <- rbeta(1, popadmixparams[1], popadmixparams[2])# random mating
  M2 <- M1 #assortative mating
  avM[individual] <- 1 - 0.5*(M1 + M2)
  obs <- simulateAutosomalGenotypes(M1, M2, rho, x, L, alleleFreqs, S)
  ##make some genotypes missing
  for(locus in 1:L) if(runif(n=1) < 0.1)
    obs[locus] <- "0,0"
  obs <- c(obs, simulateXGenotypes(M1, M2, rhoX, xX, LX, alleleFreqs[-(1:(S*L)), ], S,
                                   as.logical(male[individual])))
  genotypes <- rbind(genotypes, obs)
}

## simulate outcome
alpha <- -beta*popM
if(model==1) {# linear regression
  outcome <- numeric(N)
  for(individual in 1:N) {
    outcome[individual] <- rnorm(1, mean=(alpha + beta*avM[individual]), sd=1) 
  }
  ofam <- gaussian
} else if(model==2) {# binary outcome
  outcome <- integer(N)
  for(individual in 1:N) {
    outcome[individual] <- rbinom(1, 1, 1 / (1+exp(-alpha - beta*avM[individual])))  
  }
  ofam <- binomial
} else { # Cox regression
  outcome <- simulateCox(N, N, beta, avM)
  ofam <- poisson
}

#reg.true <-summary.glm(glm(outcome ~ avM, family = ofam))

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
loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-(1:2)]),  rep(S,L+LX), c(x, xX), chrlabels, 
                   row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "Distance", "Chr")
write.table(loci, file="simdata/loci.txt", row.names=FALSE, sep="\t")

## write priorallelefreqs file
priorallelefreqs <- matrix(data=NA, nrow=S*(L+LX), ncol=K)
for( locus in 1:(L+LX) ) {
  locus.rows <- seq(from=1+S*(locus-1), to=S*locus)
  priorallelefreqs[locus.rows, ] <- afreqparams
}
locusnames <- rep(loci[, 1], each=S)
priorallelefreqs <- data.frame(locusnames, priorallelefreqs)
dimnames(priorallelefreqs)[[2]] <- c("Locus", "Pop1", "Pop2")
write.table(priorallelefreqs, file="simdata/priorallelefreqs.txt", sep="\t", row.names=FALSE)
