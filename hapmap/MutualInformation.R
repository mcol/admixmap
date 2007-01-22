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
population <- get_option("population", args[7])
states <- get_option("states", args[8], int = TRUE)

genotypes = c("1,1", "1,2", "2,2")

results.dir = paste(population, "/", "Results", states, "States", sep = "")

# Evauluate entropy of a vector
entropy <- function(vec, na.rm = FALSE) {
	if (!na.rm) {
		na.fail(vec)
	}
	if (max(vec, na.rm = TRUE) > 1)
		stop("argument elements must be at most 1.")
	if (min(vec, na.rm = TRUE) < 0)
		stop("argument elements must be greater than 0.")
	vecsum <- sum(vec, na.rm = na.rm)
	# Since it will never sum to exactly 1, there is certain precision
	if (abs(vecsum - 1.0) > 1e-5) {
		stop("argument needs to sum to 1. Current sum is ", vecsum)
	}
	# Zeros do not contribute to entropy
	# H(a_1, a_2, ..., 0) = H(a_1, a_2, ...)
	vec <- vec[which(vec != 0)]
	return(sum(-vec * log(vec), na.rm = na.rm))
}

GP = dget(paste(results.dir, "PPGenotypeProbs.txt", sep = "/"))
orig = dget(paste(population, "chr22data", "genotypes5000_observed_dput.txt", sep = "/"))

# Debug, to look at the data by hand
gp.dbg <- function(indiv, loc) {
	print(c("original gentype:", orig[indiv,loc]))
	print(GP[,indiv,loc])
	print(c("sum:", sum(GP[,indiv,loc])))
}

# Returns a function that adds a prefix to a string
# NB: It doesn't return a string, but a function, which can be later
# called with an argument.
#
# Usage:
# > add.something <- add.prefix("something ")
# > add.something("strange")
# [1] "something strange"
# > add.something("nice")
# [1] "something nice"
#
add.prefix <- function(s) {
	return(function(x) (paste(s, x, sep = "")))
}

get.clean.cm <- function() {
	m = matrix(numeric(9), nrow = 3, ncol = 3)
	# t: true
	# c("t1,1", "t1,2", "t2,2")
	colnames(m) <- lapply(genotypes, add.prefix("t"))
	# p: predicted
	# c("p1,1", "p1,2", "p2,2")
	rownames(m) <- lapply(genotypes, add.prefix("p"))
	# It should look like this:
	# > m
	#           t1,1 t1,2 t2,2
	#      p1,1    0    0    0
	#      p1,2    0    0    0
	#      p2,2    0    0    0
	return(m)
}

# Calculates a product of marginal distributions
marginal.dist.prod <- function(my.matrix) {
	if (!is.matrix(my.matrix)) {
		stop("Argument should be a matrix.")
	}
	mcols = colSums(my.matrix, na.rm = TRUE)
	# mcols = mcols / sum(mcols)
	mrows = rowSums(my.matrix, na.rm = TRUE)
	# mrows = mrows / sum(mrows)
	# marginal distributions product
	md.prod = matrix(mrows, ncol = 1) %*% matrix(mcols, nrow = 1)
	# md.prod[which(md.prod == 0)] <- NA
	dimnames(md.prod) <- dimnames(my.matrix)
	# print(md.prod)
	return(md.prod)
}

# Calculate the mutual information as given in the paper, page 18.
# \sum_{x,y} (p_{x,y} log p_{x,y} - p_x p_y log(p_{x.} p_{.y}))
# Uses the entropy() function.

# Finds the most probable genotypes and assigns p = 1 to them, setting
# all the other genotypes to zero. Needed for fair comparison with
# fastPHASE which doesn't output any uncertainity.
sharpen.probs<- function(a) {
	for (indiv in dimnames(a)[[2]]) {
		vec = a[,indiv]
		m = max(vec)
		maxs_idx = which(vec == m)
		if (length(maxs_idx) != 1) {
			# There is more than one genotype with
			# probability equal to the max probability.
			warning("Can't find max element.", ", ", indiv, ", ", paste(vec, collapse = " "))
			# In such a case, we will throw away all but the
			# first probability.
			maxs_idx <- maxs_idx[1]
		}
		a[, indiv] <- c(0, 0, 0)
		a[maxs_idx, indiv] <- 1
	}       
	return(a)
}

joint.prob <- function(GP, locus.no, genotypes, throw_away_uncertainity = FALSE) {
	joint.pd <- get.clean.cm()
	indiv.total <- 0
	for (genotype in genotypes) {
		# message(c("Genotype start:", genotype))
		indiv.idx <- which(orig[, locus.no] == genotype)
		indiv.no <- length(indiv.idx)
		# Adding to total only individuals that have data
		indiv.total <- indiv.total + indiv.no
		# print(rowMeans(matrix(GP[, indiv.idx, locus.no], nrow = 3)))
		colname <- add.prefix("t")(genotype)
		# message(c("colname:", colname))
		indiv.probs <- matrix(GP[, indiv.idx, locus.no], nrow = 3)
		#     in1 in2 in3 in4 ...
		# 1,1 0.1 0.2
		# 1,2 0.2 0.3
		# 2,2 0.7 0.5
		# print(indiv.probs)
		if (throw_away_uncertainity) {
			indiv.probs = sharpen.probs(indiv.probs)
		}
		if (length(indiv.probs) != 0) {
			# Sum over the individuals
			joint.pd[, colname] <- (rowSums(indiv.probs))
		}
		# print(joint.pd)
		# message(c("Genotype end:", genotype))
	}
	joint.pd <- (joint.pd / indiv.total)
	return(joint.pd)
}

# Normalize columns, so the frequencies in the m  matrix are the same as
# in the s matrix

matrix.debug <- function(m, msg = "Matrix debug") {
	print(msg)
	print(dim(m))
	print(m)
	print(colSums(m))
	print(sum(m))
	print(entropy(m))
}

cols.normalize <- function(s, m) {
	# stop("halt here.")
	return(m)
}

mutual.information <- function(joint.pd) {
	# marginal distributions
	if (!is.matrix(joint.pd)) {
		stop("Argument should be a matrix.")
	}
	# TODO: Check if the matrix is square
	if (dim(joint.pd)[1] != dim(joint.pd)[2]) {
		stop("Matrix should be square.")
	}
	md.prod <- marginal.dist.prod(joint.pd)
	# mutual information
	# return(entropy(joint.pd) - entropy(md.prod))
	# return(entropy(md.prod) - entropy(joint.pd))
	# TODO: decide which way to do it: substract or divide
	return(entropy(joint.pd) / entropy(md.prod))
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

