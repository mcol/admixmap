# vim:set ft=r:
#
# author: Maciej Blizinski
# 
# Reads PPGenotypeProbs.txt file and calculates the mutual information.
#
# Should be called in batch mode, exactly like this:
#
# R CMD BATCH --no-save --no-restore --population=Eur --states=4 MutualInformation.R

# Reading some functions.
source("maskGenotypesFunctions.R")

args <- commandArgs()
chromosome <- get_option("chromosome", args[7])
# chromosome <- "Chr22"
population <- get_option("population", args[8])
# population <- "Asian"
states <- get_option("states", args[9], int = TRUE)
# states <- "4"

genotypes = c("1,1", "1,2", "2,2")

results.dir = paste(population, "/", chromosome, "Results", states, "States2", sep = "")

GP = dget(paste(results.dir, "PPGenotypeProbs.txt", sep = "/"))
orig = dget(paste(population, "chr22data", "mi_cc_observed_dput.txt", sep = "/"))
v = as.vector(as.matrix(orig))
v[which(v == "2,1")] <- "1,2"
orig = matrix(v, nrow = dim(orig)[1], ncol = dim(orig)[2])


# Debug, to look at the data by hand
gp.dbg <- function(indiv, loc) {
	print(c("original gentype:", orig[indiv,loc]))
	print(GP[,indiv,loc])
	print(c("sum:", sum(GP[,indiv,loc])))
}

# Array for mutual information.
mi <- matrix(nrow = length(GP[1,1,]), ncol = 3)
colnames(mi) <- c("Mutual information", "Mutual information, no uncertainity", "Unique genotypes")

# Unique loci in original data. Numer 1 indicates that all individuals
# had the same genotype in specific locus.

for (locus.id in dimnames(GP)[[3]]) {
	locus.no = as.numeric(locus.id)
	# print(c("Number of unique genotypes: ", length(unique(na.omit(orig[, locus.no])))))
	joint.pd <- joint.prob(GP, locus.no, genotypes)
	joint.pd.no.uncert <- joint.prob(GP, locus.no, genotypes, throw_away_uncertainity = TRUE)
	mi[locus.no, 1] <- mutual.information(joint.pd)
	mi[locus.no, 2] <- mutual.information(joint.pd.no.uncert)
	mi[locus.no, 3] <- length(unique(na.omit(orig[, locus.no])))
}

write.table(mi,
	file = paste(results.dir, "/mutual-information-by-locus.txt", sep = ""),
	row.names = TRUE, col.names = NA)
write.table(mean(mi[1]),
	file = paste(results.dir, "/mean-mutual-information.txt", sep = ""),
	col.names = FALSE, row.names = FALSE)
write.table(mean(mi[2]),
	file = paste(results.dir, "/mean-mutual-information-no-uncert.txt", sep = ""),
	col.names = FALSE, row.names = FALSE)

