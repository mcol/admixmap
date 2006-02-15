## Load the R MPI package if it is not already loaded
if (!is.loaded("mpi_initialize")) { 
#  library("Rmpi") 
} 

## Spawn as many slaves as possible 
#mpi.spawn.Rslaves() 
## In case R exits unexpectedly, have it automatically clean up 
## resources taken up by Rmpi (slaves, memory, etc...) 
.Last <- function() { 
  if (is.loaded("mpi_initialize")) { 
    if (mpi.comm.size(1) > 0){ 
      print("Please use mpi.close.Rslaves() to close slaves.") 
      mpi.close.Rslaves() 
    } 
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize") 
  } 
} 

## Tell all slaves to return a message identifying themselves 
#mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))

## Tell all slaves to close down, and exit the program 
#mpi.close.Rslaves() 
#mpi.quit()

## script to simulate data from stratified population for admixmap
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

simulateSamples <- function(N, numsims, NumSubPops, popadmixparams, rho, eta, beta, gamma, pthreshold,
                            logistic, x, L) {
  null.results <- data.frame(matrix(data=NA, nrow=0, ncol=5))
  candidate.results <- data.frame(matrix(data=NA, nrow=0, ncol=5))
  results.colnames <- c("f.signed", "crude.p", "gc.p", "adj2.p", "adj.p")
  dimnames(null.results)[[2]] <- results.colnames
  dimnames(candidate.results)[[2]] <- results.colnames
  for(sims in 1:numsims) {
    ## simulate correlated allele freqs
    mu <- numeric(L) # ancestral freqs allele 1
    alleleFreqs <- matrix(data=NA, nrow=2*L, ncol=NumSubPops)
    for(locus in 1:L) {
      mu[locus] <- rbeta(1, 2, 2)
      alleleFreqs[2*locus - 1, ] <- rbeta(2, mu[locus]*eta, (1 - mu[locus])*eta)
                                        # freqs allele 1 in each of NumSubPops subpops
      alleleFreqs[2*locus, ] <- 1 - alleleFreqs[2*locus - 1, ] # freqs allele 2
    }
    p1 <- alleleFreqs[seq(2, 2*L, by=2), 1]
    p2 <- alleleFreqs[seq(2, 2*L, by=2), 2]
    ## p1 and p2 are freqs of allele 2 in pop1 and pop2
    f.signed <- sign(p2-p1)*(p1 - p2)^2 / ((p1+p2)*(2 - p1 - p2)) 
    ## if f.signed is negative, freq allele 2 (trait-raising allele) is higher in pop2
    ## trait value increases with proportionate admixture from pop1
    ## so true effect is in opposite direction to confounding effect
    
    ## choose candidate at random
    candidate <- 1 + floor(L*runif(1)) # returns number between 1 and L
    cat("Candidate locus", candidate, "with signed f-value", f.signed[candidate], "\n")
    
    ## simulate genotypes and outcome 
    genotypes <- character(L)
    outcome <- numeric(N)
    avM <- numeric(N)
    popM <- popadmixparams[1] / sum(popadmixparams) # mean admixture proportions
    popg <- 2*(popM*alleleFreqs[2*candidate, 1] + (1 - popM)*alleleFreqs[2*candidate, 2])
                                        # expected number of copies allele 2 at candidate locus
    g.positive <- 2*popM > popg # positive confounding of g by admixture
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
    }
    genotypes <- genotypes[-1, ]
    
    ## recode genotypes as 0, 1, 2
    genotypes.r <- matrix(data=NA, nrow=dim(genotypes)[1], ncol=dim(genotypes)[2])
    for(i in 1:dim(genotypes)[1]) {
      for(j in 1:dim(genotypes)[2]) {
        if(genotypes[i, j] == "1,1") {
          genotypes.r[i, j] <- 0
        } else if(genotypes[i, j] == "1,2" | genotypes[i, j] == "2,1") {
          genotypes.r[i, j] <- 1 
        } else if(genotypes[i, j] == "2,2") {
          genotypes.r[i, j] <- 2
        }
      }
    }
    g <- genotypes.r[, candidate]
    
    for(individual in 1:N) {
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
    
    ## write outcome variable to file
    outcome.table <- data.frame(outcome, row.names=NULL) 
    write.table(outcome.table, file="data/outcome.txt", row.names=FALSE, col.names=TRUE)
    ## write true admixture proportions to file
    Mvector.table <- data.frame(avM, row.names=NULL)
    write.table(outcome, file="data/Mvalues.txt", row.names=FALSE,
                col.names=TRUE)
    ## write genotypes to file
    genotypes.gc <- genotypes
    for(col in 1:dim(genotypes)[2]) {
      genotypes.gc[, col] <- gsub(",", "\ ", as.vector(genotypes.gc[, col]))
    }  
    id = as.character(seq(1:N))
    ## write for ADMIXMAP
    genotypes <- data.frame(id, genotypes, row.names=NULL)
    write.table(genotypes, file="data/genotypes.txt", sep="\t", row.names=FALSE)
    ## write for GC
    gc.genotypes <- data.frame(id, outcome, genotypes.gc)
    write.table(gc.genotypes, file="data/gcgenotypes.txt", row.names=FALSE, col.names=FALSE,
                quote=F, sep=" ")
    ## write locus file
    x[is.na(x)] <- 100
    loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-1]),  rep(2,L),  x, row.names=NULL)
    dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "Distance")
    write.table(loci, file="data/loci.txt", row.names=FALSE)
    ## write allelefreqs files
    trueallelefreqs <- data.frame(rep(loci[,1], each=2), alleleFreqs)
    dimnames(trueallelefreqs)[[2]] <- c("Locus", "Pop1", "Pop2")
    write.table(trueallelefreqs, file="data/trueallelefreqs.txt",
                row.names=FALSE)                               
    
    ## run admixmap with outcome var and single population
    system("../test/admixmap argsSinglePop.txt")
    Sys.putenv("RESULTSDIR" = "SinglePopResults")
    source("../test/AdmixmapOutput.R")
    ## run genomic control analysis - must set L and pthreshold in gcdriver.txt
    source("gcf.R")
    ## run admixmap with no outcomevar and two populations
    system("../test/admixmap argsNoOutcome.txt")
    Sys.putenv("RESULTSDIR" = "NoOutcomeResults")
    source("../test/AdmixmapOutput.R")
    ## run admixmap with outcome var and two populations
    system("../test/admixmap argsTwoPops.txt")
    Sys.putenv("RESULTSDIR" = "TwoPopsResults")
    source("../test/AdmixmapOutput.R")
    
    ## calculate score tests for two-step structured association approach
    avM.nOutcome <- read.table("NoOutcomeResults/IndAdmixPosteriorMeans.txt", header=T)[, 2]
    r.nOutcome <- glm(outcome ~ avM.nOutcome, family=ofam)
    resid.nOutcome <- resid(r.nOutcome)
    residvar <- (sum(resid.nOutcome^2) - sum(resid.nOutcome)^2)/(N-2)
    nOutcome.pvalues <- numeric(L)
    r.pvalues <- numeric(L)
    for(locus in 1:L) {
      a <- genotypes.r[, locus] - mean(genotypes.r[, locus])
      score <- sum(a * resid.nOutcome) / residvar
      info <- sum(a^2) / residvar
      nOutcome.pvalues[locus] = 2*pnorm(-abs(score / sqrt(info)))
    }
    
    ## read adjusted and unadjusted p-values from file
    crude.pvalues <- read.table(file="SinglePopResults/TestsAllelicAssociationFinal.txt", header=T)[, 7]
    gc.pvalues <- read.table(file="GCTests.txt", header=T)[, 2]
    adj.pvalues <- numeric(L)
    adj.pvalues <- read.table(file="TwoPopsResults/TestsAllelicAssociationFinal.txt", header=T)[, 7]
    ## print posterior quantiles for params 
    reg.quantiles <- read.table(file="TwoPopsResults/PosteriorQuantiles.txt", header=T)
    print(reg.quantiles)
    ## plot posterior means of individual admixture against true values 
    avM.estimates <- read.table("TwoPopsResults/IndAdmixPosteriorMeans.txt", header=T)[, 2]
    reg.estimates <- summary.glm(glm(outcome ~ avM.estimates, family=ofam))
    plot(avM, avM.estimates, xlim=c(0,1), ylim=c(0,1))
    lines(c(0,1),c(0,1), type="l")
    
    ## bind p-values into a table
    all.results <- data.frame(f.signed, crude.pvalues, gc.pvalues, nOutcome.pvalues, adj.pvalues)
    dimnames(all.results)[[2]] <- results.colnames
    ## append to tables null.results and candidate.results
    null.results <- rbind(null.results, all.results[-candidate, ])
    candidate.results <- rbind(candidate.results, all.results[candidate, ])
    
    ## plot adjusted against crude pvalues
    postscript("TwoPopsResults/PValues.ps")
    plotchars <- numeric(L)
    plotchars[1:L] <- 1
    plotchars[candidate] <- 2
    plotcols <- character(L)
    plotcols[1:L] <- "black"
    plotcols[candidate] <- "red"
    plot(-log10(crude.pvalues), -log10(gc.pvalues), xlim=c(0,7), ylim=c(0,7),
         xlab="-log10 genomic control p-values", ylab="-log10 crude p-values",
         pch=plotchars, col=plotcols)
    plot(-log10(crude.pvalues), -log10(adj.pvalues), xlim=c(0,7), ylim=c(0,7),
         xlab="-log10 one-step adjusted p-values", ylab="-log10 crude p-values",
         pch=plotchars, col=plotcols)
    plot(-log10(nOutcome.pvalues), -log10(adj.pvalues), xlim=c(0,7), ylim=c(0,7),
         xlab="-log10 one-step adjusted p-values", ylab="-log10 two-step adjusted p-values",
         pch=plotchars, col=plotcols)
    dev.off()
  } # end simulations loop
  return(list(null.results, candidate.results))
}

##########################################################################
## start of script

numChr <- 22
## chromosome lengths in cM
chr.L <- c(292,272,233,212,197,201,184,166,166,181,156,169,117,128,110,130,128,123,109,96,59,58)
N <- 500
numsims <- 1
NumSubPops <- 2 # num subpopulations
popadmixparams <- c(1, 2) # population admixture params for pop1, pop2
rho <- 6 # sum-of-intensities
spacing <- 25 # 40 cM spacing gives 99 loci, 30 cM spacing 128 loci, 25 cM 151 loci
eta <- 51 # allele freq dispersion parameter 
beta <- 2.5 # regression slope for effect of admixture
gamma <- 0.25 # effect of allele 2 at candidate locus: standardized effect size if linear reg
                                        # log odds ratio if logistic reg
pthreshold <- 0.01
logistic <- F # logistic or linear

## assign map distances
x <- numeric(0)
chr <- numeric(0)
length <- sum(chr.L)
for(chromosome in 1:22) {
  positions <- seq(0, chr.L[chromosome], spacing)
  x <- c( x, positions) 
  chr <- c(chr, rep(chromosome, length(positions)))
}
x <- 0.01*distanceFromLast(chr, x)
L <- length(x) # number of loci: 151 if spacing is 25 cM

results <- simulateSamples(N, numsims, NumSubPops, popadmixparams, rho, eta, beta, gamma, pthreshold,
                            logistic, x, L) 
null.results <- results[[1]]
candidate.results <- results[[2]]
write.table(null.results, file="TwoPopsResults/NullResults.txt", col.names=T, row.names=F, sep="\t")
write.table(candidate.results, file="TwoPopsResults/CandidateResults.txt", col.names=T, row.names=F,
            sep="\t")

results.colnames <- c("f.signed", "crude.p", "gc.p", "adj2.p", "adj.p")
## convert p-values to error rates
type1.error <- data.frame(null.results$f.signed,
                          null.results$crude.p < pthreshold,
                          null.results$gc.p < pthreshold,
                          null.results$adj2.p < pthreshold,
                          null.results$adj.p < pthreshold)
type2.error <- data.frame(candidate.results$f.signed,
                          candidate.results$crude.p > pthreshold,
                          candidate.results$gc.p > pthreshold,
                          candidate.results$adj2.p > pthreshold,
                          candidate.results$adj.p > pthreshold)
dimnames(type1.error)[[2]] <- results.colnames
dimnames(type2.error)[[2]] <- results.colnames

## plot type 1 and type 2 error rates by quintile of signed f-value
groups <- 5
f.gr <- 1 + floor(groups*rank(type1.error$f.signed)/(1+dim(type1.error)[1]))
t1.crude <- tapply(type1.error$crude.p, f.gr, mean, na.rm=T)
t1.gc <- tapply(type1.error$gc.p, f.gr, mean, na.rm=T)
t1.adj2 <-  tapply(type1.error$adj2.p, f.gr, mean, na.rm=T)
t1.adj <-  tapply(type1.error$adj.p, f.gr, mean, na.rm=T)

f2.gr <- 1 + floor(groups*rank(type2.error$f.signed)/(1+dim(type2.error)[1]))
t2.crude <- tapply(type2.error$crude.p, f2.gr, mean, na.rm=T)
t2.gc <- tapply(type2.error$gc.p, f2.gr, mean, na.rm=T)
t2.adj2 <-  tapply(type2.error$adj2.p, f2.gr, mean, na.rm=T)
t2.adj <-  tapply(type2.error$adj.p, f2.gr, mean, na.rm=T)
Fst <- 1/(eta-1)

par(bty="l")

postscript("Type1ErrorRates.ps")
plotchars <- c(1, 2, 15, 16)
plotcols <- c("black", "red", "green", "blue")
legend.x <- rep(2, 4)
legend.labels <- c("Crude", "Genomic control", "Two-step structured association",
                   "One-step structured association")
plot(dimnames(t1.crude)[[1]], t1.crude, ylim=c(0, 0.1), type="b", 
     xlab=expression(paste("Quintile of standardized allele frequency differential with ",
         F[ST], "= 0.02")), # should substitute value of Fst 
     ylab="Type 1 error rate", pch=plotchars[1], col=plotcols[1])
points(dimnames(t1.gc)[[1]], t1.gc,  pch=plotchars[2], col=plotcols[2])
points(dimnames(t1.adj2)[[1]], t1.adj2,  pch=plotchars[3], col=plotcols[3])
points(dimnames(t1.adj)[[1]], t1.adj,  pch=plotchars[4], col=plotcols[4])
legend.y <- seq(0.1, 0.07, by=-0.01)
text(legend.x, legend.y, labels=legend.labels, adj=c(0, 0.5))   
points(legend.x - 0.2, legend.y, pch=plotchars, col=plotcols)
dev.off()

postscript("Type2ErrorRates.ps")
plot(dimnames(t2.crude)[[1]], t2.crude, ylim=c(0, 1), type="b",
     xlab=expression(paste("Quintile of standardized allele frequency differential with ",
         F[ST], "= 0.02")), # should substitute value of Fst 
     ylab="Type 2 error rate", pch=plotchars[1], col=plotcols[1])
points(dimnames(t2.gc)[[1]], t2.gc,  pch=plotchars[2], col=plotcols[2])
points(dimnames(t2.adj2)[[1]], t2.adj2,  pch=plotchars[3], col=plotcols[3])
points(dimnames(t2.adj)[[1]], t2.adj,  pch=plotchars[4], col=plotcols[4])
legend.y <- seq(0.9, 0.75, by=-0.05)
text(legend.x, legend.y, labels=legend.labels, adj=c(0, 0.5))   
points(legend.x - 0.2, legend.y, pch=plotchars, col=plotcols)
dev.off()

print(apply(type1.error[, -1], 2, mean, na.rm=T))
print(apply(type2.error[, -1], 2, mean, na.rm=T))
