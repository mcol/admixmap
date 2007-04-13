# author: Maciej Blizinski

# Converts data in hapmapmix format into fastPHASE format.
#
# Data shoud be given as two objects:
#
# 1. genotypes, a matrix where rows are individuals and columns are
# diploid genotypes defined by values "1,1", "1,2", "2,2"
# 
# 2. loci, a matrix where rows are loci and the second column stores
# distance between current and previous loci.
#
# Typical call: save.fastphase(genotypes, loci, "fastphase.inp")

# Decode diploid data, two functions
# Return first genotype
#
# It could be also done with substr(x, 1, 1), but this seemingly stupid
# way is actually good for returning NA values. It also assures that
# only "1,1", "1,2", "2,1" and "2,2" values are recognized and warnings
# are provided if a different value is spotted.

get.gt.genotype <- function(x, n) {
	if (!(n == 1 || n == 2)) {
		stop("Second argument should be either 1 or 2.")
	}
	if (is.na(x)) { return(NA) }
	else if (x == "1,1") { return(c(1,1)[n]) }
	else if (x == "1,2") { return(c(1,2)[n]) }
	else if (x == "2,1") { return(c(2,1)[n]) }
	else if (x == "2,2") { return(c(2,2)[n]) }
	else if (x == "0,0") { return(NA) }
	else {
		stop("Unrecognized value: ", x)
		return(NA)
	}
}

# An array of 3 lines per individual
indivs.fastphase <- function(fp) {
	il <- array() # Individual lines
	i = 1
	for (indiv in 1:length(fp[,1,1])) {
		# message("Individual: ", indiv)
		il[i] <- paste("# id ", indiv)
		il[i+1] <- paste(fp[indiv,,1], collapse = "")
		il[i+2] <- paste(fp[indiv,,2], collapse = "")
		i = i + 3
	}
	return(il)
}

save.fastphase <- function(genotype, loci, out.fastphase.file) {
	number.indivs = length(genotype[,1])
	number.loci = length(genotype[1,])
	stopifnot(number.loci == length(loci[,2]))

	# Three-dimensional array, dimensions:
	# 1. Individual
	# 2. Locus
	# 3. Genotype, for two fastPHASE data file lines
	fp <- array(dim = c(number.indivs, number.loci, 2))

	# Map the diploid data with get.first() and get.second() functions.
	# message("Finding first genotype line.")
	fp[, , 1] <- matrix(sapply(as.matrix(genotypes),  get.gt.genotype, 1),
		ncol = number.loci, nrow = number.indivs)
	# message("Finding second genotype line.")
	fp[, , 2] <- matrix(sapply(as.matrix(genotypes), get.gt.genotype, 2),
		ncol = number.loci, nrow = number.indivs)

	# Replace NA's with "?" character
	fp[which(is.na(fp))] <- "?"

	loci[1, 2] <- 0 # fix problem with the initial NA value, it doesn't hurt.

	cat(
		number.indivs,
		number.loci,
		# diffinv does the integration, so distances between loci are
		# converted to offset positions
		# UPDATE: Locus position line removed becuase it's not
		# used anyway and fastPHASE won't read the data
		# properly.
		#
		# UPDATE: This line is necessary when using the -M2
		# option as stated in fastPHASE manual, section 5.10,
		# page 17.
		paste("P", paste(diffinv(loci[,2]), collapse = " ")),
		indivs.fastphase(fp),
		sep = "\n",
		file = out.fastphase.file)
}
