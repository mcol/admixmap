## script to calculate bayesian info reward
## given allele counts in individuals not masked, true values of masked genotypes,
## and predictive probs

###################################################################################

entropy <- function(q) {
  invsum <- 1 / sum(q)
  q <- q * invsum
  u <- q  * log(q)
  u[q == 0] <- 0
  return(-sum(u))
}

entropy.array <- function(probs) {
  ## probs is an array of dim K x N
  ## returns array of length N
  N <- dim(probs)[2]
  h <- numeric(N)
  for(i in 1:N) {
    h[i] <- entropy(probs[, i])
  }
  return(h) 
}

entropy.array3d <- function(probs) {
  ## probs is an array of dim K x N x M
  ## returns array of length N, each element is mean of M entropies 
  N <- dim(probs)[2]
  h <- numeric(N)
  for(i in 1:N) {
    h[i] <- mean(entropy.array(probs[, i, ]), na.rm=T)
  }
  return(h)
}

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

   ## normalize prior probability vector
   prior <- prior / sum(prior)
   if (is.na(sum(predictive)) | sum(predictive) == 0) { ## no prediction, return zero
     i.reward <- 0
   } else { # normalize predictive prob vector
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
  ## loops over loci and returns the mean bayesian info reward
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
    cat(" predictiveprobs3d", dim(predictiveprobs3d), "\n")
    return(1)
  } 

  bir <- NULL
  for(locus in 1:dim(counts2d)[2]) { 
    prior <- priorfromcounts(counts2d[, locus])
    bir <- c(bir, bir.locus(prior, predictiveprobs3d[,locus,], truevalues2d[, locus]))
  }
  return(mean(bir))
}

read.fastphase <- function(fpfilename, masked.gametes, masked.loci) {
  ## masked.gametes is a logical vector specifying which gametes are from diploid individuals
  ## whose missing genotypes are to be imputed
  ## returns a 3-way array of genotype probs
  masked.gametes <- rev(masked.gametes) ## tested individuals are written first by fp
  numgametes <- length(masked.gametes) 
  fplines <- scan(file=fpfilename, comment.char=">", sep="", what="character", quiet=T)
  ## fplines is vector: each line specifies a gamete
  numrows <- length(fplines); 
  numsamples <- numrows/numgametes
  test.lines <- rep(masked.gametes, numsamples)
  fplines <- fplines[test.lines] ## subset rows for tested gametes
  ## vector now has length 2 * numtested.indivs * numsamples

  ## split haplotype strings into elements
  numloci <- length(strsplit(fplines[1], split="")[[1]])
  haps <- matrix(data=NA, nrow=length(fplines), ncol=numloci)
  for(line in 1:length(fplines)) {
    haps[line, ] <- as.integer(strsplit(fplines[line], split="")[[1]])
  }
  haps <- haps[, masked.loci] ## subset tested loci

  numtested.indivs <- length(masked.gametes[masked.gametes]) / 2
  nummasked.loci <- length(masked.loci[masked.loci])
  numloci <- nchar(fplines[1])
  ## read haps into a 3-way array of diploid genotypes
  g <- array(data=NA, dim=c(numsamples, numtested.indivs, nummasked.loci))
  for(j in 1:numsamples) {
    for(i in 1:numtested.indivs) { # code genotypes as 1, 2, 3
      g[j, i, ] <- haps[(j - 1)*2*numtested.indivs + 2*i - 1, ] +
                   haps[(j - 1)*2*numtested.indivs + 2*i, ]     
    }
  }
  g <- g + 1
  
  inv.numsamples <- 1 / numsamples
  ## calculate genotype probs
  gprobs <- array(data=NA, dim=c(3, nummasked.loci, numtested.indivs))
  for(i in 1: numtested.indivs) {
    for(j in 1:nummasked.loci) {
      gprobs[, j, i] <- table(factor(g[, i, j], levels=1:3))
      
    }
  }
  gprobs <- inv.numsamples * gprobs
  return(gprobs)
}

#######################################################################################

which.chr <- "chr22_5kloci"
#which.chr <- "chr22"

dataprefix <- paste("data", which.chr, sep="/")
resultsprefix <- "/exports/work/scratch/pmckeigu/"

## remove any old output file
system("[ -e birResults.txt ] && rm birResults.txt")

popnames <- c("YRI", "CEU", "JPTCHB")

for( pop in popnames){
  ## 1: hapmixmap results
  
  ## write pathname of PPGenotypeProbs.txt to filenames.txt
  system(paste("find ", resultsprefix, which.chr, "/hapmixmap/", pop,
               " -name PPGenotypeProbs.txt > filenames.txt", sep=""))
  ## read the filenames 
  filenames <- scan(file="filenames.txt", what="character", quiet=T)
  
  ## set data directory for this population
  datadir <- paste(dataprefix, "hapmixmap", pop, sep="/")
  
  ## get true values from data directory
  truevalues2d <- dget(paste(datadir, "obs_masked_genotypes.txt",sep="/"))
  locusnames.truevalues <- dimnames(truevalues2d)[[2]]
  ## get observed counts from data directory
  counts2d <- dget(paste(datadir, "obs_allele_counts.txt", sep = "/"))
  
  priors <- matrix(data=NA, nrow=3, ncol=dim(counts2d)[2])
  for(locus in 1:dim(priors)[2]) {
    priors[, locus] <- priorfromcounts(counts2d[, locus])
  }
  priors.entropy <- entropy.array(priors)
  avpriorh <- mean(priors.entropy)
  cat(pop, "prior", "-", "-", "-", "-", "-", "-", "-", avpriorh, "\n")
  cat(pop, "prior", "-", "-", "-", "-", "-", "-", "-", avpriorh, "\n", file="birResults.txt", append=T)
  
  if(length(filenames) > 0) { # loop over directories containing hapmixmap results
    for(run in 1:length(filenames)) {
      run.name <- strsplit(strsplit(filenames[run], "/")[[1]][9], "_")[[1]]
      predictiveprobs3d <- dget(as.character(filenames[run]))
      ## check that locus names match
      if(length(locusnames.truevalues[locusnames.truevalues!=dimnames(predictiveprobs3d)[[3]]]) > 0) {
        print(paste(filenames[run],
                    "locus names mismatch between true values and hapmixmap predictive probs"))
      }
      bir.result <-  bir.loci(counts2d, predictiveprobs3d, truevalues2d)
      predictive.entropy <- entropy.array3d(predictiveprobs3d)
                                        #plot(priors.entropy, predictive.entropy)
      avh <- mean(predictive.entropy)
      cat(pop, run.name, dim(predictiveprobs3d), bir.result, avh/avpriorh, "\n")
      ## append results to file
      cat(pop, run.name, dim(predictiveprobs3d), bir.result, avh/avpriorh, "\n", file="birResults.txt", append=T)
    }
    ## loop over hapmixmap runs merging runs that differ only in the seed
    ## first create a vector of filenames with seed stripped off
    filepath <- gsub("_seed-[[:digit:]]", "", filenames)     
    filepath <- unique(filepath)
    filepath <- gsub("/PPGenotypeProbs.txt", "", filepath)     
    
    for(run in 1:length(filepath)) {
      run.name <- strsplit(strsplit(filepath[run], "/")[[1]][9], "_")[[1]]
      
      ## get vector of filenames for which options match current run.name
      findcmd <- paste("find ", filepath[run], "* -name PPGenotypeProbs.txt > seeds.txt", sep="")
      system(findcmd)
      seeds <- as.character(read.table(file="seeds.txt", header=F)[, 1])
      ## loop over these runs to average predictiveprobs3d
      predictiveprobs3d <- dget(seeds[1])
      for(s in 2:length(seeds)) {
        predictiveprobs3d <- predictiveprobs3d + dget(as.character(seeds[s]))
      }
      predictiveprobs3d <- 1/length(seeds) * predictiveprobs3d
      
      ## check that locus names match
      if(length(locusnames.truevalues[locusnames.truevalues!=dimnames(predictiveprobs3d)[[3]]]) > 0) {
        print(paste(filenames[run],
                    "locus names mismatch between true values and hapmixmap predictive probs\n"))
      }
      bir.result <-  bir.loci(counts2d, predictiveprobs3d, truevalues2d)
      predictive.entropy <- entropy.array3d(predictiveprobs3d)
      avh <- mean(predictive.entropy)
      cat(pop, run.name, dim(predictiveprobs3d), bir.result, avh/avpriorh, "\n")
      ## append results to file
      cat(pop, run.name, "-", dim(predictiveprobs3d), bir.result, avh/avpriorh, "\n", file="birResults.txt", append=T)
    }
  }
  
  ## impute results
  ##
  impute.outfile <- paste("results", which.chr, "impute", pop, "out.txt", sep="/")
  ##read in IMPUTE output, using SNP names in col2 as row names and drop first 5 cols
  gprobs <- read.table(file=impute.outfile, as.is=T, row.names=2)[, -c(1:4)]
  gprobs <- gprobs[row.names(gprobs) %in% locusnames.truevalues, ]
  ## check that locus names match
  if(length(locusnames.truevalues[locusnames.truevalues!=row.names(gprobs)]) > 0) {
    print(paste(pop,
                "locus names mismatch between true values and impute predictive probs"))
  }
  numloci <- dim(gprobs)[1]
  N <- dim(gprobs)[2] / 3 # num masked individuals 
  predictiveprobs3d <- array(as.numeric(as.matrix(gprobs)), dim=c(numloci, 3, N))
  predictiveprobs3d <- aperm(predictiveprobs3d, c(2, 1, 3))
  
  ## flatten prob array to eliminate zeros
  predictiveprobs3d <- 0.0001 + predictiveprobs3d
  
  bir.result <-  bir.loci(counts2d, predictiveprobs3d, truevalues2d)
  h <- entropy.array3d(predictiveprobs3d) # returns NaN where all elements are 0
  h[is.nan(h)] <- priors.entropy[is.nan(h)] # replace these elements with prior entropy
  avh <- mean(h)
  cat(pop, "impute", "-", "-", "-", dim(predictiveprobs3d), bir.result, avh/avpriorh, "\n")
  cat(pop, "impute", "-", "-", "-", dim(predictiveprobs3d), bir.result, avh/avpriorh, "\n",
      file="birResults.txt", append=T)

  ## fastphase results
  ##
  tempdir <- "/exports/work/scratch/pmckeigu"
  fpfilename <- paste(tempdir, which.chr, "fastphase", pop, "_sampledHgivG.txt", sep="/")

  all.locusnames <- read.table(file=paste(dataprefix, "hapmixmap", pop, "train_loci.txt", sep="/"),
                               na.strings="#", header=T, as.is=T, fill=T)[, 1]
  masked.loci <- match(locusnames.truevalues, all.locusnames)
  if(pop == "JPTCHB"){
    masked.gametes <- c(rep(F, 150), rep(T, 30))
  }else{
    masked.gametes <- c(rep(F, 100), rep(T, 20))
  }
  predictiveprobs3d <- read.fastphase(fpfilename, masked.gametes, masked.loci)
  predictiveprobs3d <- 0.0001 + predictiveprobs3d 
  bir.result <-  bir.loci(counts2d, predictiveprobs3d, truevalues2d)
  avh <- mean(entropy.array3d(predictiveprobs3d))
  cat(pop, "fastphase", "-", "-", "-", dim(predictiveprobs3d), bir.result, avh/avpriorh, "\n")
  cat(pop, "fastphase", "-", "-", "-", dim(predictiveprobs3d), bir.result, avh/avpriorh, "\n",
      file="birResults.txt", append=T)
}
