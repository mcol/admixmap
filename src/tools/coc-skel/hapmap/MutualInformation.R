# vim:set ft=r:
#
# author: Maciej Blizinski
# 
# Reads PPGenotypeProbs.txt file and calculates the mutual information.
#
# Should be called in batch mode, exactly like this:
#
# R CMD BATCH --no-save --no-restore \
# --chromosome=Chr22 --population=Eur --states=4 MutualInformation.R

# Reading some functions.
source("maskGenotypesFunctions.R")

args <- commandArgs()
chromosome <- get.option("chromosome", args[7])
# chromosome <- "Chr22"
population <- get.option("population", args[8])
# population <- "Eur"
states <- get.option("states", args[9], int = TRUE)
# states <- "8"

genotypes <- c("1,1", "1,2", "2,2")

results.dir <- paste(population, "/", chromosome, "Results", states, "States2", sep = "")

print(results.dir)

# Dimensions of Genotype Probs are:
#
# 1. Genotype
# 2. Locus
# 3. Individual
GP <- dget(paste(results.dir, "PPGenotypeProbs.txt", sep = "/"))
orig <- dget(paste(population, "chr22data", "mi_cc_observed_dput.txt", sep = "/"))
v <- as.vector(as.matrix(orig))
v[which(v == "2,1")] <- "1,2"
orig <- matrix(v, nrow = dim(orig)[1], ncol = dim(orig)[2])


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
mi <- matrix(nrow = length(GP[1,,1]), ncol = 3)
colnames(mi) <- c(
	"Mutual information",
	"Mutual information, no uncertainity",
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
	# mi[locus.no, ] <- coc.locus(locus.no, GP, orig, genotypes)
	mi.tmp <- coc.locus(locus.no, GP, orig, genotypes)
	mi[locus.no, 1] <- mi.tmp[1]
	mi[locus.no, 2] <- mi.tmp[2]
	mi[locus.no, 3] <- mi.tmp[3]
	# print(as.matrix(cbind(c("Locus no.", locus.no), as.vector(mi[locus.no, ]))))
	# print((cbind(t(c("Locus no.", locus.no)), t(as.vector(mi[locus.no, ])))))
}

# mi <- get.coc.table (GP, orig, genotypes)

# Unique loci in original data. Numer 1 indicates that all individuals
# had the same genotype in specific locus.

# for (locus.id in 1:length(dimnames(GP)[[2]]) {

write.table(mi,
	file = paste(results.dir, "coefficient-of-constraint-by-locus.txt", sep = "/"),
	row.names = TRUE, col.names = NA)
dput(mi, file = paste(results.dir, "coefficient-of-constraint-by-locus-dput.txt", sep = "/"))
write.table(mean(mi[,1], na.rm = TRUE),
	file = paste(results.dir, "mean-coefficient-of-constraint.txt", sep = "/"),
	col.names = FALSE, row.names = FALSE)
write.table(mean(mi[,2], na.rm = TRUE),
	file = paste(results.dir, "mean-coefficient-of-constraint-no-uncert.txt", sep = "/"),
	col.names = FALSE, row.names = FALSE)

