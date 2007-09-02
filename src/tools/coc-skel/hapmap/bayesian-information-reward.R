# vim:set ft=r:
#
# author: Maciej Blizinski
# 
# Reads PPGenotypeProbs.txt file and calculates bayesian information
# reward.
#
# Should be called like this:
#
# R --vanilla --args datadir=data_dir resultsdir=output_dir \
#   <bayesian-information-reward.R [>logfile]

# Read some functions from external script
source("maskGenotypesFunctions.R")

##retrieve command-line arguments
##args is the full command used to invoke R
##typically /path/to/Rterm <R args like --vanilla> [--args [args to script]] 
args <- commandArgs();
##check if --args used. This avoids a problem with earlier versions of R
##or when script run without '--args'. args.pos is the index of '--args' in the args array.
args.pos <- match("--args", args)

print(args)
print(args.pos)

if(is.null(args) || is.na(args.pos) || length(args) != args.pos +2){
  cat ("Incorrect syntax: should be\n",
##       "R CMD BATCH --vanilla --args datadir=... resultsdir=... bayesian-information-reward.R\n",
##  "or\n",
  "R --vanilla --args datadir=... resultsdir=... <bayesian-information-reward.R\n")
  q("no")
}

##data.dir <- args[args.pos + 1]
##results.dir <- args[args.pos + 2]
data.dir <- get.option("datadir", args[args.pos+1])
results.dir <- get.option("resultsdir", args[args.pos+2])

# Dimensions of Genotype Probs are:
#
# 1. Genotype
# 2. Locus
# 3. Individual
GP <- dget(paste(results.dir, "PPGenotypeProbs.txt", sep = "/"))
obs.genotypes <- dget(paste(data.dir, "obs_masked_genotypes.txt", sep = "/"))
prior <- dget(paste(data.dir, "genotype_freqs.txt", sep = "/"))

# Debug, to look at the data by hand
gp.dbg <- function(indiv, loc) {
	print(c("original gentype:", orig[indiv,loc]))
	print(GP[,indiv,loc])
	print(c("sum:", sum(GP[,indiv,loc])))
}

# This, when called as a function, works horribly slowly.
# mi <- get.coc.table(GP, orig, genotypes)

# Array for Bayesian information reward

BIR <- NULL
masked.loci <- dimnames(GP)[[2]]
for (locus.name in masked.loci) {
  if (is.na(GP[1,locus.name,1])) {
    next
  }
      
  ## print(c("Number of unique genotypes: ", length(unique(na.omit(orig[, locus.no])))))
  BIR <- rbind(BIR, bir.locus2(locus.name, GP, obs.genotypes, prior))      
  ##}
}
dimnames(BIR) <- list(masked.loci, c("BIR","BIR, no uncertainity","Unique genotypes"))


write.table(BIR,
            file = paste(results.dir, "bayesian-information-reward-by-locus.txt",
              sep = "/"),
            row.names = TRUE, col.names = NA)
dput(BIR, file = paste(results.dir, "bayesian-information-reward-by-locus-dput.txt",
            sep = "/"))

write.table(mean(BIR[,1], na.rm = TRUE),
            file = paste(results.dir, "mean-bayesian-information-reward.txt", sep = "/"),
            col.names = FALSE, row.names = FALSE)

write.table(mean(BIR[,2], na.rm = TRUE),
            file = paste(results.dir, "mean-bayesian-information-reward-no-uncert.txt",
              sep = "/"),
            col.names = FALSE, row.names = FALSE)

