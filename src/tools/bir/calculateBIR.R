# vim:set ft=r:
#
# authors: Maciej Blizinski, David O'Donnell, Paul McKeigue
# 
# Reads PPGenotypeProbs.txt file and calculates bayesian information
# reward.
#
# Should be called like this:
#
# R --vanilla --args data_dir output_dir \
#   <bayesian-information-reward.R [>logfile]


####################################
## retrieve command-line arguments
#####################################
## args is the full command used to invoke R
## typically /path/to/Rterm <R args like --vanilla> [--args [args to script]] 
args <- commandArgs();
##check if --args used. This avoids a problem with earlier versions of R
##or when script run without '--args'. args.pos is the index of '--args' in the args array.
args.pos <- match("--args", args)

if(is.null(args) || is.na(args.pos) || length(args) != args.pos +2){
  cat ("Incorrect syntax: should be\n",
##       "R CMD BATCH --vanilla --args datadir=... resultsdir=... bayesian-information-reward.R\n",
##  "or\n",
  "R --vanilla --args {datadir} {resultsdir} <bayesian-information-reward.R\n")
  q("no")
}

data.dir <- args[args.pos + 1]
results.dir <- args[args.pos + 2]
##data.dir <- get.option("datadir", args[args.pos+1])
##results.dir <- get.option("resultsdir", args[args.pos+2])

#########################
## function definitions
#########################
priorfromcounts <- function(counts) {
  ## argument is integer vector of length 2 containing counts of each allele
  ## returns probability vector of length 3 
  ## draw 10000 samples from posterior distribution of population allele freqs
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
   ## true value.  predictive[t]= 1-predictive[-t] = 1
   ## evaluates to  - log(prior[t]) + log(1 - prior[-t]) 
   return(i.reward)
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
  ## this function loops over loci and returns the bayesian info reward as a vector
  ## counts - dim 2 x numloci
  ## predictiveprobs3d - dim 3 x numloci x N
  ## truevalues2d - dim N x numloci
  
  ## check array dimensions for consistency
  if(dim(predictiveprobs3d)[1] != 3 ||
     dim(predictiveprobs3d)[3] != dim(truevalues2d)[1] ||
     dim(predictiveprobs3d)[2] != dim(truevalues2d)[2] ||
     dim(counts2d)[1] !=2 ||
     dim(counts2d)[2] != dim(truevalues2d)[2]) {
    print("inconsistent array dimensions")
    return(1)
  }

  masked.loci <- dimnames(predictiveprobs3d)[[2]]
  bir <- NULL
  for (locus in masked.loci) {
    if(counts2d[1,locus]>0 && counts2d[2,locus]>0){##skip monomorphic loci
      prior <- priorfromcounts(counts2d[, locus])
      bir <- c(bir, bir.locus(prior, predictiveprobs3d[,locus ,], truevalues2d[, locus]))
    }else{
      bir <- c(bir, 0)
    }
  }
  return(bir)
}

# Debug, to look at the data by hand
gp.dbg <- function(indiv, loc) {
	print(c("original gentype:", obs[indiv,loc]))
	print(GP[,loc,indiv])
	print(c("sum:", sum(GP[,indiv,loc])))
}

#######################################
## script starts here
######################################

# Dimensions of Genotype Probs are:
#
# 1. Genotype
# 2. Locus
# 3. Individual
GP <- dget(paste(results.dir, "PPGenotypeProbs.txt", sep = "/"))
obs.genotypes <- dget(paste(data.dir, "obs_masked_genotypes.txt", sep = "/"))
#prior <- dget(paste(data.dir, "genotype_freqs.txt", sep = "/"))
allele.counts <- dget(paste(data.dir, "obs_allele_counts.txt", sep = "/"))

# Array for Bayesian information reward
BIR <- bir.loci(allele.counts, GP, obs.genotypes)
  
write(mean(BIR, na.rm = TRUE),
      file = paste(results.dir, "mean-bayesian-information-reward.txt", sep = "/"),
      )

                      
write.table(data.frame(locus=(dimnames(GP)[[2]]), BIR),
      file = paste(results.dir, "bayesian-information-reward-by-locus.txt",
        sep = "/"),
      row.names=F, col.names=T)

#write.table(mean(BIR[,2], na.rm = TRUE),
#            file = paste(results.dir, "mean-bayesian-information-reward-no-uncert.txt",
#              sep = "/"),
#            col.names = FALSE, row.names = FALSE)

