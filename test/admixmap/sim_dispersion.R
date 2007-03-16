## script to generate data from a dispersion model


ro <- 5
x <- 0.03
L <- 100
I <- 400

simulateHaploidAlleles <- function(M,ro,x,L) {
  T <- t(matrix(data=c(M + (1-M)*exp(-ro*x),  1-M - (1-M)*exp(-ro*x),
                M - M*exp(-ro*x),      1-M + M*exp(-ro*x)    ),
              nrow=2, ncol=2))
  gameteAncestry <- numeric(L)
  randAnc <- runif(L)
  y <- numeric(L)
  alleleFreqs <- numeric(2)
  gameteAncestry[1] <- ifelse(randAnc[1] > M, 2, 1) # M is prob pop 1.
  alleleFreqs[1] <- rbeta( 1, 80, 20 )
  alleleFreqs[2] <- rbeta( 1, 200, 800 )
  y[1] <- ifelse(runif(1) > alleleFreqs[gameteAncestry[1]], 2, 1)
  for(locus in 2:L) {
# allele frequencies in admixed population
    alleleFreqs[1] <- rbeta( 1, 80, 20 )
    alleleFreqs[2] <- rbeta( 1, 200, 800 )
    gameteAncestry[locus] <- ifelse(randAnc[locus] > T[gameteAncestry[locus-1], 1],
                                    2, 1)
    y[locus] <- ifelse(runif(1) > alleleFreqs[gameteAncestry[locus]], 2, 1)
  }
  return(y)
}

simulateGenotypes <- function(Mfather,Mmother,ro,x,L) {
  paternalGamete <- simulateHaploidAlleles(Mfather,ro,x,L)
  maternalGamete <- simulateHaploidAlleles(Mmother,ro,x,L)
  simulateGenotypes <- paste(paternalGamete, ",", maternalGamete, sep="")
}

genotypes <- character(L)
for(individual in 1:I) {
  M <- rbeta(1, 4, 1)
  obs <- simulateGenotypes(M, M, ro, x, L)
  genotypes <- rbind(genotypes, obs)
}
genotypes <- genotypes[-1,]
id=seq(1:I)
genotypes <- data.frame(id, genotypes, row.names=NULL)
write.table(genotypes, file="c:/cvs2/genepi/sims/dispersion/data/genotypes.txt", row.names=FALSE)

loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-1]),
                   rep(2,L), c(100, rep(x, L-1)),
                   row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "Distance")
write.table(loci, file="c:/cvs2/genepi/sims/dispersion/data/loci.txt", row.names=FALSE)

ewalleleCounts <- matrix(nrow=2,ncol=2)
allelecounts <- numeric(2)
locusnames <- character(0)

# allele frequencies in parental population
newalleleCounts <- matrix(nrow=2,ncol=2)
  alleleFreqs.afr.par <- rbeta( L, 80, 20 )
  alleleFreqs.eur.par <- rbeta( L, 200, 800 )

for(locus in 1:L) {
# simulate allele counts in parental population
  newalleleCounts[1,1] <- rbinom( 1, 1000, alleleFreqs.afr.par[locus] )
  newalleleCounts[2,1] <- 1000 - newalleleCounts[1,1]
  newalleleCounts[1,2] <- rbinom( 1, 1000, alleleFreqs.eur.par[locus] )
  newalleleCounts[2,2] <- 1000 - newalleleCounts[1,2]

  allelecounts <- rbind(allelecounts, newalleleCounts)
  locusnames <- c(locusnames, rep(as.vector(loci[locus,1]), 2))
}
allelecounts <- data.frame(locusnames, allelecounts[-1,], row.names=NULL)
write.table(allelecounts, file="c:/cvs2/genepi/sims/dispersion/data/allelecounts.txt", row.names=FALSE)
