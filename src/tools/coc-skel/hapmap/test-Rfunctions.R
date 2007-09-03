

priorfromcounts <- function(counts) {
  ## argument is integer vector of length 2 containing counts of each allele
  ## returns probability vector of length 3 
  ## fix this later to marginalize over population freqs
  p <- (counts[1] + 0.5) / (sum(counts) + 1)
  q <- 1 - p
  return(c(p^2, 2*p*q, q^2))
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
    print("inconsistent array dimensions")
    return(1)
  } 
  bir <- numeric(numloci)
  for(locus in 1:numloci) {
    prior <- priorfromcounts(counts2d[, locus])
    bir[locus] <- bir.locus(prior, predictiveprobs3d[, , locus], truevalues2d[, locus])
  }
  return(mean(bir))
}


## call bir.loci with following arguments:
## counts2d : 2 x numloci matrix of allele counts
## predictiveprobs3d: 3 x N x numloci array of predictive probs
## truevalues2d: N x numloci array of true values

N <- 2
numloci <- 4

counts2d <- matrix(data = c(10, 12, 10, 30,
                     10, 20, 10, 20), nrow=2)
predictiveprobs3d <- array(data = c(0.5, 0.25, 0.25, 0.5, 0.25, 0.25,
                             0.5, 0.25, 0.25, 0.5, 0.25, 0.25,
                             0.5, 0.25, 0.25, 0.5, 0.25, 0.25,
                             0.5, 0.25, 0.25, 0.5, 0.25, 0.25), dim=c(3, 2, 4))
truevalues2d <- matrix(data=c(1, 2, 3, 1, 2, 3, 2, 3), nrow=2)

cat("BIR is ", bir.loci(counts2d, predictiveprobs3d, truevalues2d), "\n")
