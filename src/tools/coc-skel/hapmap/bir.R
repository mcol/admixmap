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
  ## predictiveprobs3d - dim 3 x N x numloci
  ## truevalues2d - dim N x numloci
  ## check array dimensions for consistency
  if(dim(predictiveprobs3d)[1] != 3 ||
     dim(predictiveprobs3d)[2] != dim(truevalues2d)[1] ||
     dim(predictiveprobs3d)[3] != dim(truevalues2d)[2] ||
     dim(counts2d)[1] !=2 ||
     dim(counts2d)[2] != dim(truevalues2d)[2]) {
    cat("inconsistent array dimensions\n")
    cat("truevalues2d", dim(truevalues2d))
    cat("counts2d", dim(counts2d))
    cat("predictiveprobs3d", dim(predictiveprobs3d))
    return(1)
  } 
  bir <- numeric(numloci)
  for(locus in 1:numloci) {
    prior <- priorfromcounts(counts2d[, locus])
    bir[locus] <- bir.locus(prior, predictiveprobs3d[, , locus], truevalues2d[, locus])
  }
  return(mean(bir))
}


#######################################################################################

## from console on walton, create a tar file containing genotype probs
## find -name "PPGenotypeProbs.txt" | xargs tar rvf PPGenotypeProbs.tar
## from local console, download the tar file
## scp -r -P 8022 pmckeigue@walton-local:/ichec/work/ndlif006b/maciej/PPGenotypeProbs.tar .

popnames <- c("Afr", "Eur", "Asian")
# loop over 3 populations
for(pop in 1:3) {
  ## write pathnames to file
  ## irw-2 runs appear to be comparisons of different priors
  system(paste("find * -name PPGenotypeProbs.txt | grep irw-2 | grep",
               popnames[pop], "> filenames.txt"))
  filenames <- read.table(file="filenames.txt", header=F, as.is=T)[, 1]
  runs.table <- read.table(file="filenames.txt", header=F, sep="/", as.is=T)[, -c(1, 6)]
  priors2 <- t(matrix(unlist(strsplit(runs.table[, 1], "afpp-")), nrow=2))
  priors2 <- t(matrix(unlist(strsplit(runs.table[, 1], "afpp-")), nrow=2))
  arrival.priors <- t(matrix(unlist(strsplit(priors2[, 1],    "arp-" )), nrow=2))[, 2]
  dispersion.priors <- t(matrix(unlist(strsplit(priors2[, 2],    "-s" )), nrow=2))[, 1]
  param.priors <- data.frame(arrival.priors, dispersion.priors, stringsAsFactors=F)


  ## set data directory for this population
  datadir <- paste("chr22-5k-data", popnames[pop], sep="/")

  ## get true values and priors
  truevalues <- as.matrix(dget(paste(datadir, "mi_cc_observed_dput.txt", sep="/")))
  N <- dim(truevalues)[1]
  numloci <- dim(truevalues)[2]
  truevalues <- as.vector(truevalues)
  truevalues2d <- integer(N * numloci)
  truevalues2d[truevalues=="1,1"] <- 1
  truevalues2d[truevalues=="1,2" | truevalues=="2,1"] <- 2
  truevalues2d[truevalues=="2,2"] <- 3
  truevalues2d <- matrix(data=truevalues2d, nrow=N)
  
  ## these freqs appear to be calculated from 60 diploid individuals: should have
  ## excluded the 10 indivs with masked genotypes
  gfreqs <- dget(paste(datadir,"genotype-freqs.R", sep="/"))
  counts1 <-  2*gfreqs[, 1] + gfreqs[, 2]
  counts2 <-  gfreqs[, 2] + 2*gfreqs[, 3]
  counts2d <- 60 * rbind(counts1, counts2)

  for(run in 1:length(filenames)) {
    predictiveprobs3d <- aperm(dget(filenames[run]),c(1, 3, 2))
    bir.result <-  bir.loci(counts2d, predictiveprobs3d, truevalues2d)
    priors <- c(param.priors[run, 1], param.priors[run, 2])
    cat(priors, dim(predictiveprobs3d), bir.result, "\n")
    ## append results to file
    #cat(priors, dim(predictiveprobs3d), bir.result, "\n", file="birResults.txt", append=T)
  }
}
