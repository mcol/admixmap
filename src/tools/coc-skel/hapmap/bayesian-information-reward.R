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

print (data.dir)
print (results.dir)
q("no")

genotypes <- c("1,1", "1,2", "2,2")

# Dimensions of Genotype Probs are:
#
# 1. Genotype
# 2. Locus
# 3. Individual
GP <- dget(paste(results.dir, "PPGenotypeProbs.txt", sep = "/"))
orig <- dget(paste(data.dir, "mi_cc_observed_dput.txt", sep = "/"))
v <- as.vector(as.matrix(orig))
v[which(v == "2,1")] <- "1,2"
orig <- matrix(v, nrow = dim(orig)[1], ncol = dim(orig)[2])
prior <- dget(paste(data.dir, "genotype-freqs.R", sep = "/"))

# Debug, to look at the data by hand
gp.dbg <- function(indiv, loc) {
	print(c("original gentype:", orig[indiv,loc]))
	print(GP[,indiv,loc])
	print(c("sum:", sum(GP[,indiv,loc])))
}

# This, when called as a function, works horribly slow.
# mi <- get.coc.table(GP, orig, genotypes)
nan_idx <- as.vector(which(is.na(GP[1,,1])))
# Array for mutual information.
BIR <- matrix(nrow = length(GP[1,,1]), ncol = 3)
colnames(BIR) <- c(
	"BIR",
	"BIR, no uncertainity",
	"Unique genotypes")

for (locus.no in 1:length(GP[1,,1])) {
	if (length(nan_idx) > 0) {
		# FIXME: This is an ad-hoc fix that only works for 1 NA
		# locus. If there are more NA loci, it will fail.
		if (locus.no == nan_idx[1]) {
			next
		}
	}
	# print(c("Number of unique genotypes: ", length(unique(na.omit(orig[, locus.no])))))
	joint.pd <- joint.prob(GP, locus.no, genotypes)
	joint.pd.no.uncert <- joint.prob(GP, locus.no, genotypes, throw_away_uncertainity = TRUE)
	# BIR[locus.no, ] <- coc.locus(locus.no, GP, orig, genotypes)
	BIR.tmp <- bir.locus(locus.no, GP, orig, genotypes, prior)
	BIR[locus.no, 1] <- BIR.tmp[1]
	BIR[locus.no, 2] <- BIR.tmp[2]
	BIR[locus.no, 3] <- BIR.tmp[3]
}

# BIR <- get.coc.table (GP, orig, genotypes)

# Unique loci in original data. Numer 1 indicates that all individuals
# had the same genotype in specific locus.

# for (locus.id in 1:length(dimnames(GP)[[2]]) {

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

