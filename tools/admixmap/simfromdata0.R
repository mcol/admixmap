#script to simulate test data for admixmap from data
## currently only for diallelic, simple loci (no composite loci) and 2 populations
# modify these to match the data
#NumLoci <- 40
NumInd <- 100
NPops <- 2
inputdir <- "c:/cvs/genepi/test/data"
locusfile <- "loci.txt"
ingenotypesfile <- "genotypes.txt"
inallelefreqfile <- "priorallelefreqs.txt"
outputdir <- "c:/cvs/genepi/sims/testsimdata"
outgenotypesfile <- "simgenotypes.txt"
historicallelefreqfile <- "histAllelefreqs.txt"
priorallelefreqfile <- "simpriorallelefreqs.txt"
randommatingmodel <- FALSE
fixedallelefreqs <- TRUE
missinggenotypes <- FALSE
globalrho <- FALSE
simulateHaploidAlleles <- function(M,ro,x,Loci, AlleleFreqs) {
                                        #x = locus spacings
                                        #M = prob of pop 1
                                        #rho is sumintensities
  gameteAncestry <- numeric(Loci)
  randAnc <- runif(Loci)
  gameteAncestry[1] <- ifelse(randAnc[1] > M, 2, 1)
                                        # M is prob pop 1.gameteAncestry[1]
  HaploidAlleles <- ifelse(runif(1) > AlleleFreqs[1,1,gameteAncestry[1]], 2, 1)
  for(locus in 2:Loci) {
    T <- t(matrix(data=c(M + (1-M)*exp(-ro*x[locus]),  1-M - (1-M)*exp(-ro*x[locus]),
                    M - M*exp(-ro*x[locus]), 1-M + M*exp(-ro*x[locus]) ), nrow=2, ncol=2))
    gameteAncestry[locus] <- ifelse(randAnc[locus] > T[gameteAncestry[locus-1], 1], 2, 1)
    HaploidAlleles <- c(HaploidAlleles,ifelse(runif(1) > AlleleFreqs[locus,1,gameteAncestry[locus]], 2, 1))  }
  return(HaploidAlleles)
}
simulateGenotypes <- function(M1,M2, rho1, rho2,xx,L, AlleleFreqs) {
  simulateGenotypes <- NULL
  for(chr in 1:length(L))#loop over chromosomes
    {    if(chr==1)seq <- 1:L[1] else seq <- (sum(L[1:(chr-1)])+1):sum(L[1:chr])
         x <- xx[seq]
         paternalGamete <- simulateHaploidAlleles(M1,rho1,x,L[chr], AlleleFreqs[seq,,])
         maternalGamete <- simulateHaploidAlleles(M2,rho2,x,L[chr], AlleleFreqs[seq,,])
         simulateGenotypes <- c(simulateGenotypes,paste(paternalGamete, ",", maternalGamete, sep=""))
       }
  return(simulateGenotypes)
}

## function to tabulate allele counts
countAlleles <- function(genotypes) {
  obs <- table(genotypes, exclude="")
  genotype.labels <- dimnames(obs)[[1]]
  obs11 <- ifelse("1,1" %in% genotype.labels, obs[match("1,1", genotype.labels)], 0)
  obs12 <- ifelse("1,2" %in% genotype.labels, obs[match("1,2", genotype.labels)], 0)
  obs22 <- ifelse("2,2" %in% genotype.labels, obs[match("2,2", genotype.labels)], 0)
  count1 <- 2*obs11 + obs12
  count2 <- 2*obs22 + obs12
  return(c(count1, count2))
}
## read genotypes
tablecheckGenotypes <- function(genotypes, freqs){
  counts <- data.frame(0, 0)
  for(locus in 1:dim(genotypes)[2]) {
    counts <- rbind(counts, countAlleles(genotypes[, locus]))  }
  counts <- counts[-1, ]
  freqs1.AfrAm <- counts[, 1] / apply(counts, 1, sum)
  ## read allelefreqs table
  LocusNames <- freqs[seq(from=1, to=(dim(freqs)[1] - 1), by=2), 1]
  freqs.1 <- freqs[seq(from=1, to=(dim(freqs)[1] - 1), by=2), -1]
  freqs.2 <- freqs[seq(from=2, to=dim(freqs)[1], by=2), -1]
  freqs1 <- freqs.1 / (freqs.1 + freqs.2)
  ## plot observed against expected allele freqs
  freqs1.Expected <- 0.8*freqs1[, 1] + 0.2*freqs1[, 2]
  plot(freqs1.Expected, freqs1.AfrAm)
}

#################### Start of script
locus.table <- read.table(paste(inputdir,locusfile,sep="/"),header=TRUE, row.names=NULL, colClasses=c("character", "numeric", "numeric"))
NumLoci <- nrow(locus.table)
distances <- locus.table[,3]#markerspacings
#L <- NumLoci#number of loci
#determine number of chromosomes and loci on each
NumChromosomes <- sum(distances==100)
L <- numeric(NumChromosomes)
c <- 0
for(i in 1:length(distances)){
  if(distances[i]==100)c <- c+1
  L[c] <- L[c]+1
}
allelefreq.table <- read.table(paste(inputdir, inallelefreqfile, sep="/"), header=TRUE, row.names=NULL, colClasses=c("character", rep("numeric",NPops)))
#TODO: account for input historicallelefreqs by adding 0.5 to allelecounts
alleleCounts <- array(data=rep(0,NumLoci*2*NPops),dim=c(NumLoci, 2, NPops))
alleleFreqs <- array(data=rep(0,NumLoci*2*NPops),dim=c(NumLoci, 2, NPops))
#create allelefreqs
for(locus in 1:NumLoci){
  alleleCounts[locus,,] <- as.matrix(allelefreq.table[(locus*2-1):(locus*2),2:(NPops+1)]-0.5)
#alleleFreqs <- matrix(c(1, 0, 0, 1), nrow=2, ncol=2)
#matrix(c(0.8, 0.2, 0.2, 0.8), nrow=2, ncol=2)
#alleleFreqs[locus,,] <- matrix(c(0.65, 0.0, 1.0, 0.35), nrow=2, ncol=2)
                          
  for(k in 1:NPops){
    if(fixedallelefreqs){
      alleleFreqs[locus,,k] <- alleleCounts[locus,,k] /sum(alleleCounts[locus,,k])    }else{# generate allelefreqs from prior
        freq <- rbeta(1,allelefreq.table[locus*2-1,1+k],allelefreq.table[locus*2,1+k] )
        alleleFreqs[locus,,k] <- c(freq, 1.0 - freq)    }  }}
#generate probs of Pop1 (African)
#use parameters 8 and 2 here for African/European
#TODO: allow user-specified proportions
#generate genotypes
genotypes <- NULL
for(individual in 1:NumInd) {
  M1 <- rbeta(1, 8, 2)#paternal
  if(randommatingmodel) M2 <- rbeta(1, 8, 2) else M2 <- M1#maternal
  rhopat <- rgamma(1,6,1)#6
  if(globalrho)rhomat <- rhopat else rhomat <- rgamma(1,6,1)
  obs <- simulateGenotypes(M1, M2, rhopat, rhomat, distances, L, alleleFreqs)
#make some (~20%) genotypes missing
#TODO: reproduce pattern of missing data in input genotypes file
#if(missinggenotypes)for(locus in 1:L)if(runif(n=1) < 0.2)obs[locus]<-""
  genotypes <- rbind(genotypes, obs)
}

if(missinggenotypes){
  ingenotypes  <- read.table(paste(inputdir,ingenotypesfile,sep="/"), header=TRUE, nrow=NumInd)[, -1]
  genotypes[ingenotypes==""] <- ""
  rm(ingenotypes)
}
# genotypes file
#genotypes <- genotypes[-1,]
id <- seq(1:NumInd)#indiv ids
genotypes <- data.frame(id, genotypes, row.names=NULL)
dimnames(genotypes)[[2]][-1] <- locus.table[,1]#locus names
write.table(genotypes, file=paste(outputdir,outgenotypesfile, sep="/"), row.names=FALSE)

#write locus file
write.table(locus.table, file=paste(outputdir, "loci.txt", sep="/"), row.names=F)

#create allelefreqs files
#old format allelefreqfile
#allelefreqs <- data.frame(as.vector(locus.table[,1]), alleleFreqs[,1,1], alleleFreqs[,1,2])
#dimnames(allelefreqs)[[2]] <- c("Locus", "Pop1", "Pop2")
#write.table(allelefreqs, file=paste(outputdir, outallelefreqfile, sep="/"), row.names=FALSE)

allelefreqs <- numeric(2)
locusnames <- character(0)

for(locus in 1:NumLoci) {
  allelefreqs <- rbind(allelefreqs, alleleCounts[locus,,])
  locusnames <- c(locusnames, rep(as.vector(locus.table[locus,1]), 2))
}
historicallelefreqs <- data.frame(locusnames, allelefreqs[-1,], row.names=NULL)
priorallelefreqs <- data.frame(locusnames, allelefreqs[-1,]+0.5, row.names=NULL)

#priorallelefreqfile
write.table(priorallelefreqs, file=paste(outputdir, priorallelefreqfile, sep="/"), row.names=FALSE)
#historicallelefreqfile
write.table(historicallelefreqs, file=paste(outputdir, historicallelefreqfile, sep="/"), row.names=FALSE)

checkGenotypes(genotypes, historicallelefreqs)#TODO: generate outcomevar and covariates files
