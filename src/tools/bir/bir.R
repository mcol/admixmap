## script to calculate bayesian info reward
## given allele counts in individuals not masked, true values of masked genotypes,
## and predictive probs

###################################################################################

priorfromcounts <- function(counts) {
  ## argument is integer vector of length 2 containing counts of each allele
  ## returns probability vector of length 3 
  ## draw 10000 samples from posterior distribution of population allele freqs
  ## given reference prior + counts
  p <- rbeta(10000, counts[1] + 0.5, counts[2] + 0.5)
  ## average genotype freqs assuming HWE
  gprobs11 <- mean(p^2)
  gprobs22 <- mean((1 - p)^2)
  return(c(gprobs11, 1 - gprobs11 - gprobs22, gprobs22))
}

info.reward2 <- function(prior, predictive, t) {
   ## prior is a probability vector of length K
   ## predictive is a probability vector of length K
   ## t is the correct class: integer between 1 and K
   K <- length(prior)

   ## fix this later to normalize by dividing by max value given the prior
   
   ## normalize prior probability vector
   prior <- prior / sum(prior)
   if (sum(predictive) == 0) { ## no prediction, return zero
     i.reward <- 0
   } else {
     predictive <- predictive / sum(predictive)
     ## i.plus evaluates to -Inf if predictive[t] is zero
     i.plus <- log(predictive[t] / prior[t])
     ## i.minus evaluates to -Inf if any element of predictive[-t] is 1
     i.minus <- sum( log( (1 - predictive[-t]) / (1 - prior[-t]) ) )
     i.reward <- (i.plus + i.minus) / K
   }
   ## max value of info.reward2 given prior is for a confident prediction of the
   ## true value, such that predictive[t]= 1 - predictive[-t] = 1
   i.reward.max <- - log(prior[t]) - sum(log(1 - prior[-t]))
   return(i.reward / i.reward.max)
}

bir.locus <- function(prior, predictiveprobs2d, truevalues) {
  ## this function loops over individuals and returns the mean
  ## bayesian information reward at the locus
  N <- dim(predictiveprobs2d)[2]
  bir <- numeric(N)
  for(i in 1:N) {
    bir[i] <- info.reward2(prior, predictiveprobs2d[, i], truevalues[i])
  }
  return(mean(bir))  
}
  
bir.loci <- function(counts2d, predictiveprobs3d, truevalues2d) {
  ## this function loops over loci and returns the mean bayesian info reward
  ## counts - dim 2 x numloci
  ## predictiveprobs3d - dim 3 x numloci x N
  ## truevalues2d - dim N x numloci
  ## check array dimensions for consistency
  if(dim(predictiveprobs3d)[1] != 3 ||
     dim(predictiveprobs3d)[3] != dim(truevalues2d)[1] ||
     dim(predictiveprobs3d)[2] != dim(truevalues2d)[2] ||
     dim(counts2d)[1] !=2 ||
     dim(counts2d)[2] != dim(truevalues2d)[2]) {
    cat("array dimensions inconsistent \n")
    cat(" truevalues2d", dim(truevalues2d))
    cat(" counts2d", dim(counts2d))
    cat(" predictiveprobs3d", dim(predictiveprobs3d))
    return(1)
  } 

  masked.loci <- dimnames(predictiveprobs3d)[[2]]
  bir <- NULL
  for(locus in masked.loci) {
    prior <- priorfromcounts(counts2d[, locus])
    bir <- c(bir, bir.locus(prior, predictiveprobs3d[,locus,], truevalues2d[, locus]))
  }
  return(mean(bir))
  ##return bir
}

read.fastphase <- function(fpfilename, tested.gametes, tested.loci) {
  ## tested.gametes is a logical vector specifying which gametes are from diploid individuals
  ## whose missing genotypes are to be imputed
  ## returns a 3-way array of genotype probs
  fplines <- scan(file=fpfilename, comment.char=">", sep="", what="character", quiet=T)
  numloci <- nchar(fplines[1])

  numrows <- length(fplines); 
  numgametes <- length(tested.gametes) 
  numsamples <- numrows/numgametes

  testrows <- rep(tested.gametes, numsamples)
  fplines <- fplines[testrows] # should be of length 2*numtested.indivs*numsamples
  
  numtested.indivs <- length(tested.gametes[tested.gametes]) / 2
  g <- matrix(data=NA, nrow=numtested.indivs, ncol=numloci)
  for(i in 1:numtested.indivs) {
    g[i, ] <- as.integer(strsplit(fplines[2*i - 1], split="")[[1]]) +
      as.integer(strsplit(fplines[2*i],     split="")[[1]]) + 1
  }
  g <- g[, tested.loci]
  numtested.loci <- length(tested.loci[tested.loci])
  g <- array(g, dim=c(numsamples, numtested.indivs, numtested.loci))
  
  ## check that this fills correctly
  inv.numsamples <- 1 / numsamples
  gprobs <- inv.numsamples * apply(g, 2:3, table)
  
  ## gprobs should have dims 3, numtested.loci, numtested.indivs  
  return(aperm(gprobs, c(1, 3, 2)))
}

#######################################################################################

which.chr <- "chr22_5kloci"
#which.chr <- "chr22"

dataprefix <- paste("data", which.chr, sep="/")

## remove any old output file
system("[ -e birResults.txt ] && rm birResults.txt")

popnames <- c("YRI", "CEU", "JPTCHB")

for( pop in popnames){
## 1: hapmixmap results
  
  ## write pathnames to file
  system(paste("find results/", which.chr, "hapmixmap/", pop, "/ -name PPGenotypeProbs.txt >filenames.txt"))
  filenames <- scan(file="filenames.txt", what="character", quiet=T)

  ## set data directory for this population
  datadir <- paste(dataprefix, "hapmixmap", pop, sep="/")

  ## get true values from data directory
  truevalues2d <- dget(paste(datadir, "obs_masked_genotypes.txt",sep="/"))
  locusnames.truevalues <- dimnames(truevalues2d)[[2]]
  ## get observed counts from data directory
  counts2d <- dget(paste(datadir, "obs_allele_counts.txt", sep = "/"))
  
  ##
  ## loop over directories containing hapmixmap results
  ## 
  for(run in 1:length(filenames)) {
    run.name <- strsplit(strsplit(filenames[run], "/")[[1]][5], "_")[[1]]
    ##states <- run.name[1]
    ##arrival.priors <- run.name[2]
    ##dispersion.priors <- run.name[3]
    ##seed <- run.name[4]
    
    predictiveprobs3d <- dget(as.character(filenames[run]))
    ## check that locus names match
    if(length(locusnames.truevalues[locusnames.truevalues!=dimnames(predictiveprobs3d)[[3]]]) > 0) {
      print(paste(filenames[run],
                  "locus names mismatch between true values and hapmixmap predictive probs\n"))
    }
    bir.result <-  bir.loci(counts2d, predictiveprobs3d, truevalues2d)

    cat(run.name, pop, dim(predictiveprobs3d), bir.result, "\n")
    ## append results to file
    cat(run.name, pop, dim(predictiveprobs3d), bir.result, "\n", file="birResults.txt", append=T)
  }
  
  ##
  ## impute results
  ##
  impute.outfile <- paste("results", which.chr, "impute", pop, "out.txt", sep="/")

  ##read in IMPUTE output, using SNP names in col2 as row names and drop first 5 cols
  gprobs <- read.table(file=impute.outfile, row.names=2)[, -c(1:4)]
  gprobs <- gprobs[row.names(gprobs) %in% locusnames.truevalues, ]
  ## check that locus names match
  if(length(locusnames.truevalues[locusnames.truevalues!=row.names(gprobs)]) > 0) {
    print(paste(pop,
                "locus names mismatch between true values and impute predictive probs\n"))
  }
  numloci <- dim(gprobs)[1]
  N <- dim(gprobs)[2] / 3 # num masked individuals 
  predictiveprobs3d <- aperm(array(as.matrix(gprobs), dim=c(numloci, N, 3)), c(3, 1, 2))
  bir.result <-  bir.loci(counts2d, predictiveprobs3d, truevalues2d)
  cat("impute", pop, dim(predictiveprobs3d), bir.result, "\n")

  ##
  ## fastphase results
  ##  
  fpfilename <- paste("results", which.chr, "fastphase", pop, "fastphase_sampledHgivG.txt", sep="/")

  masked.loci <- scan(paste(dataprefix, "hapmixmap", pop, "masked_loci.txt", sep="/"), quiet=T)
  if(pop == "JPTCHB"){
    masked.gametes <- c(rep(F, 150), rep(T, 30))
  }else{
    masked.gametes <- c(rep(F, 100), rep(T, 20))
  }
  ##predictiveprobs3d <- read.fastphase(fpfilename, masked.gametes, masked.loci)
  ##bir.result <-  bir.loci(counts2d, predictiveprobs3d, truevalues2d)
  ##cat("fastphase", pop, dim(predictiveprobs3d), bir.result, "\n")
  

  
}
