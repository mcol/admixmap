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
    if(is.na(x[locus])) { ## distance from last missing
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

simulateGenotypes <- function(sex, M1, M2, rho, x, L, Xchr.L, alleleFreqs) {
  # returns a vector of genotypes
  maternalGamete <- simulateHaploidAlleles(M2,rho,x,L+Xchr.L, alleleFreqs)  
  paternalGamete <- simulateHaploidAlleles(M1,rho, x, L+Xchr.L, alleleFreqs)
  diploidAlleles <- paste(paternalGamete, maternalGamete, sep=",") 
  if(sex==1) { ## haploid at X chr loci
    diploidAlleles[(L+1):(L+Xchr.L)] <- gsub(",[12]", "", diploidAlleles[(L+1):(L+Xchr.L)])
  }
  return(diploidAlleles)
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

PolyaLogL <- function(eta, mu, n) {
  ## K outcomes, N observations with same mu
  ## mu is vector of length K, n is matrix of counts with K rows and N cols
  LogL <- 0
  K <- length(mu)
  N <- dim(n)[1]
  for(i in 1:N) { # outer sum over N observations
    LogL <- LogL + lgamma(eta) - lgamma(sum(n[i, ]) + eta)
    for(k in 1:K){ # inner sum over outcomes 1:K
      LogL <- LogL + lgamma(n[i, k] + eta*mu[k]) - lgamma(eta*mu[k]) # outcome 1
    }
  }
  return(LogL)
}

PolyaLogLOverLoci <- function(eta, mu, counts) {
  ## K outcomes, M experiments with different mu each of which has N observations 
  ## mu is matrix with M rows, K cols
  ## counts is 3-way array with dim N, M, K
  M <- dim(mu)[1]
  LogL <- 0
  for(m in 1:M) { # M experiments with different mu
    LogL <- LogL + PolyaLogL(eta, mu[m, ], counts[m,,])
  }
  return(LogL)
}
    
  
##########################################################################
## Start of script
numChr <- 22
## chromosome lengths in cM
chr.L <- c(292,272,233,212,197,201,184,166,166,181,156,169,117,128,110,130,128,
           123,109,96,59,58,120)  # last value is length of X chr
N <- 200
sex <- 2 - rbinom(N, 1,  0.5)
NumSubPops <- 2 # num subpopulations
popadmixparams <- c(2, 6) # population admixture params for pop1, pop2
rho <- 6 # sum-of-intensities
spacing <- 40 # 40 cM spacing gives 99 autosomal loci
eta <- 5 # allele freq dispersion parameter #10 is upper limit with 200 obs and admixmparams Di(1,2)
beta <- 2 # regression slope for effect of admixture
gamma <- 0.4 # effect of allele 2 at candidate locus: standardized effect size if linear reg
                                        # log odds ratio if logistic reg
logistic <- TRUE # logistic or linear

## assign map distances
x <- numeric(0)
chr <- numeric(0)
length <- sum(chr.L)
chr.labels <- c(as.character(1:22), "X")
for(chromosome in 1:23) {
  positions <- seq(0, chr.L[chromosome], spacing)
  x <- c( x, positions) 
  chr <- c(chr, rep(chr.labels[chromosome], length(positions)))
}
print(data.frame(x, chr))

chrnum <- chr
chrnum[chr=="X"] <- "23"
dist <- distanceFromLast(as.integer(chrnum), x)
Xchr.L <- length(positions) # number of X chr loci
L <- length(x) - Xchr.L # number of autosomal loci

null.results <- data.frame(matrix(data=NA, nrow=0, ncol=4))
candidate.results <- data.frame(matrix(data=NA, nrow=0, ncol=4))
results.colnames <- c("f.signed", "crude.p", "gc.p", "adj.p")
dimnames(null.results)[[2]] <- results.colnames
dimnames(candidate.results)[[2]] <- results.colnames
admixmap <- TRUE
numsims <- 1

## simulate allele freqs
mu <- numeric(L+Xchr.L) # ancestral freqs allele 1
alleleFreqs <- matrix(data=NA, nrow=2*(L+Xchr.L), ncol=NumSubPops)
for(locus in 1:(L+Xchr.L)) {
  mu[locus] <- rbeta(1, 2, 2)
  alleleFreqs[2*locus - 1, ] <- rbeta(2, mu[locus]*eta, (1 - mu[locus])*eta)
                                        # freqs allele 1 in each of NumSubPops subpops
  alleleFreqs[2*locus, ] <- 1 - alleleFreqs[2*locus - 1, ] # freqs allele 2
}
alleleFreqs[,1] <- 0.8
alleleFreqs[,2] <- 0.2

p1 <- alleleFreqs[seq(2, 2*(L+Xchr.L), by=2), 1]
p2 <- alleleFreqs[seq(2, 2*(L+Xchr.L), by=2), 2]
## positive f-value if freq allele 2 higher in pop1 than pop2
f.signed <- sign(p1-p2)*(p1 - p2)^2 / ((p1+p2)*(2 - p1 - p2)) 
## if f.signed is negative, true effect is in opposite direction to confounding effect

## choose candidate at random
candidate <- 1 + floor(L*runif(1)) # returns number between 1 and L
## choose a candidate locus from centile of signed f-values
##candidate <- match(floor(0.05*L),   rank(f.signed)) # 5th centile of f-value: negative confounding
##candidate <- match(floor(0.95*L), rank(f.signed)) # 95th centile of f-value: positive confounding
# cat("Candidate locus", candidate, "with signed f-value", f.signed[candidate], "\n")

## simulate genotypes and outcome 
genotypes <- character(L+Xchr.L)
outcome <- numeric(N)
avM <- numeric(N)
g <- numeric(N)
popM <- popadmixparams[1] / sum(popadmixparams) # mean admixture proportions
popg <- 2*(popM*alleleFreqs[2*candidate, 1] + (1 - popM)*alleleFreqs[2*candidate, 2])
                                        # mean number of copies allele 2 at candidate locus
g.positive <- 2*popM > popg # positive confounding of g by admixture 

for(individual in 1:N) {
  M1 <- 1 - rbeta(1, popadmixparams[1], popadmixparams[2]) ## M1 is prob pop 1
  ## M2 <- rbeta(1, 4, 1)# random mating
  M2 <- M1 #assortative mating
  avM[individual] <- 1 - 0.5*(M1 + M2)
  obs <- simulateGenotypes(sex[individual], M1, M2, rho, dist, L, Xchr.L, alleleFreqs)
  ## recode genotype at candidate locus
  if(obs[candidate] == "1,1") {
    g[individual] <- 0
  } else if(obs[candidate] == "1,2" | obs[candidate] == "2,1") {
    g[individual] <- 1 
  } else if(obs[candidate] == "2,2") {
    g[individual] <- 2
  }
  ##make some genotypes missing
  ##for(locus in 1:L) if(runif(n=1) < 0.1)
  ##    obs[locus]<-"0,0"
  genotypes <- rbind(genotypes, obs)
  
  ## simulate outcome
  alpha <- -beta*popM - gamma*popg 
  if(logistic) { # logistic regression with approx equal numbers of cases and controls
    outcome[individual] <-
      rbinom(1, 1, 1 / (1+exp(-(alpha + beta*avM[individual] + gamma*g[individual]))))  
    ofam <- binomial
  } else { # linear regression
    outcome[individual] <-
      rnorm(1, mean=(alpha + beta*avM[individual]+ gamma*g[individual]), sd=1) 
    ofam <- gaussian
  }
}
genotypes <- genotypes[-1, ] # drop empty first row

## use mkdirs.default in R.utils package to create dir

## write outcome variable to file
outcome.table <- data.frame(outcome, row.names=NULL) 
write.table(outcome.table, file="data/outcome.txt", row.names=FALSE, col.names=TRUE)
## write true admixture proportions to file
Mvector.table <- data.frame(avM, row.names=NULL)
write.table(outcome, file="data/Mvalues.txt", row.names=FALSE,
            col.names=TRUE)
id = as.character(seq(1:N))

## write genotypes
row.names(genotypes) <- NULL
genotypes <- data.frame(id, genotypes) #, row.names=NULL)
write.table(genotypes, file="data/genotypes.txt", sep="\t", row.names=FALSE, quote=FALSE)

## write locus file
x[is.na(x)] <- NA
loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-1]),
                   rep(2, L + Xchr.L),  dist, chr, row.names=NULL)
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
locusnames <- paste("X", seq(1:(L+Xchr.L)), sep="")
priorallelefreqs <- data.frame(locusnames, priorallelefreqs)
write.table(priorallelefreqs, file="data/priorallelefreqs.txt", row.names=FALSE, sep="\t",
            quote=FALSE)


  
