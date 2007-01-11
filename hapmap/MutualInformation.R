# vim:set ft=r:
#
# FIXME: This file needs to be checked. Calculated mutual information is
# very close to zero in the results. For instance:
#
#  [1]  9.656637e-09  1.982607e-07 -1.319238e-06 -6.433266e-07 -6.290709e-09
#  [6] -4.539540e-07  0.000000e+00 -2.605956e-07  0.000000e+00  0.000000e+00
# [11]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -1.021872e-07
#
# I noticed that for many loci all the individuals have the same
# genotype. In my sample data:
#
# 1 genotype:  266 loci
# 2 genotypes: 174 loci
# 3 genotypes:  60 loci

region = "Eur"
states = 4

genotypes = c("1,1", "1,2", "2,2")

results.dir = paste(region, "/", "Results", states, "States", sep = "")

# Evauluate entropy of a vector
entropy <- function(vec) {
	if (max(vec, na.rm = TRUE) > 1)
		stop("vector elements must be at most 1.")
	if (min(vec, na.rm = TRUE) <= 0)
		stop("vector elements must be greater than 0.")
	return(sum(-vec * log(vec), na.rm = TRUE))
}

GP = dget(paste(results.dir, "PPGenotypeProbs.txt", sep = "/"))
orig = dget(paste(region, "chr22data", "genotypes5000_observed_dput.txt", sep = "/"))

# Debug, to look at the data by hand
gp.dbg <- function(indiv, loc) {
	print(c("original gentype:", orig[indiv,loc]))
	print(GP[,indiv,loc])
	print(c("sum:", sum(GP[,indiv,loc])))
}

# Returns a function that prepends a string
add.prefix <- function(s) {
	return(function(x) (paste(s, x, sep = "")))
}

get.clean.cm <- function() {
	m = matrix(nrow = 3, ncol = 3)
	# t: true
	# c("t1,1", "t1,2", "t2,2")
	colnames(m) <- lapply(genotypes, add.prefix("t"))
	# p: predicted
	# c("p1,1", "p1,2", "p2,2")
	rownames(m) <- lapply(genotypes, add.prefix("p"))
	# It should look like this:
	# > m
	#           t1,1 t1,2 t2,2
	#      p1,1   NA   NA   NA
	#      p1,2   NA   NA   NA
	#      p2,2   NA   NA   NA
	return(m)
}

# Calculate the mutual information as given in the paper, page 18.
# \sum_{x,y} (p_{x,y} log p_{x,y} - p_x p_y log(p_{x.} p_{.y}))
# Uses the entropy() function.
mutual.information <- function(locus.no, genotypes) {
	for (genotype in genotypes) {
		message(c("Genotype start:", genotype))
		joint.pd <- get.clean.cm()
		indiv.idx <- which(orig[, locus.no] == genotype)
		# print(rowMeans(matrix(GP[, indiv.idx, locus.no], nrow = 3)))
		colname = add.prefix("t")(genotype)
		message(c("colname:", colname))
		joint.pd[, colname] <- rowMeans(matrix(GP[, indiv.idx, locus.no], nrow = 3))
		# print(joint.pd)
		message(c("Genotype end:", genotype))
	}
	# print(joint.pd)
	# marginal distributions
	mcols = colSums(joint.pd, na.rm = TRUE)
	mcols = mcols / sum(mcols)
	mrows = rowSums(joint.pd, na.rm = TRUE)
	mrows = mrows / sum(mrows)
	# marginal distributions product
	md.prod = matrix(mrows, nrow = 3) %*% matrix(mcols, ncol = 3)
	dimnames(md.prod) <- dimnames(joint.pd)
	md.prod[which(md.prod == 0)] <- NA
	# mutual information
	# print(md.prod)
	return(sum(joint.pd * log(joint.pd) - md.prod * log(md.prod), na.rm = TRUE))
}

# Array for mutual information.
mi <- array()

# Unique loci in original data. Numer 1 indicates that all individuals
# had the same genotype in specific locus.
uniloc <- array()

for (locus.id in dimnames(GP)[[3]]) {
	locus.no = as.numeric(locus.id)
	# print(c("Number of unique genotypes: ", length(unique(na.omit(orig[, locus.no])))))
	uniloc[locus.no] <- length(unique(na.omit(orig[, locus.no])))
	mi[locus.no] <- mutual.information(locus.no, genotypes)
}

