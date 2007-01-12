##script to mask genotypes (set them to missing) in order to assess prediction of missing genotypes in HAPMIXMAP

source("maskGenotypesFunctions.R")

########################################################################
## Configuration
########################################################################

# Base name for all the input and output files.

io.file.path <- "Eur/chr22data"
io.file.basename <- "genotypes5000"

##note:assuming genotypes to be masked are all diploid 
missing.genotype <- "\"0,0\""

in.genotypes.file <- get.io.filename(".txt")
in.loci.file <- paste(io.file.path, "loci5000.txt", sep = "/")
out.genotypes.file <- get.io.filename("_masked.txt")
out.index.file <- get.io.filename("_index.txt")
out.observed.genotypes.file <- get.io.filename("_observed.txt")
out.observed.genotypes.dput.file <- get.io.filename("_observed_dput.txt")
out.fastphase.file <- get.io.filename("_fastphase.inp")

percent.missing.indivs <- 10
percent.missing.loci <- 10

#################################################################

genotypes <- read.table( # returns a data.frame
	in.genotypes.file,
	header = TRUE,
	na.strings = c("\"0,0\"", "0,0", "\"0\"", "0"),
	colClasses = "character",
	row.names = "Individ")

number.loci <- ncol(genotypes)
number.indivs <- nrow(genotypes)

missing.loci <- sort(
                   sample(
                   c(1:number.loci),
	           size = round((number.loci * percent.missing.loci / 100.0)),
	           replace = FALSE))

missing.indivs <- sort(
                       sample(
	               c(1:number.indivs),
	               size = round((number.indivs * percent.missing.indivs / 100.0)),
	               replace = FALSE))

# Save original genotypes
# Useful for debug
genepi.write.table(genotypes, get.io.filename("_original_full.txt"))

# Save genotypes to be masked as a table
# This file should be small.
genepi.write.table(genotypes[missing.indivs, missing.loci], out.observed.genotypes.file)
# Save genotypes to be masked as R object
dput(genotypes[missing.indivs, missing.loci], out.observed.genotypes.dput.file)

# Erase the appropriate data by inserting missing values
genotypes[missing.indivs, missing.loci] <- NA

# Save the (big) changed file.
genepi.write.table(genotypes, out.genotypes.file)

cat(
	paste("maskedindivs = ", missing.indivs),
	paste("maskedloci = ", missing.loci),
	sep = "\n",
	file = out.index.file)

########################################################################
#  Write the data in fastPHASE format
########################################################################

loci = read.table(
	in.loci.file,
	header = TRUE,
	sep = "\t",
	row.names = "SNPid")

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

# Three dimensional array, dimensions:
# 1. Individual
# 2. Locus
# 3. Genotype, for two fastPHASE data file lines
fp <- array(dim = c(number.indivs, number.loci, 2))

# Tiny data
# genotypes = genotypes[1:4, 1:5]

# Map the diploid data with get.first() and get.second() functions.
fp[, , 1] <- matrix(sapply(as.matrix(genotypes),  get.gt.genotype, 1),
	ncol = number.loci, nrow = number.indivs)
fp[, , 2] <- matrix(sapply(as.matrix(genotypes), get.gt.genotype, 2),
	ncol = number.loci, nrow = number.indivs)

# Replace NA's with "?" character
fp[which(is.na(fp))] <- "?"

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

loci[1, 2] <- 0 # fix problem with the NA value, it doesn't hurt.

cat(
	number.indivs,
	number.loci,
	# diffinv does the integration, so distances between loci are
	# converted to offset positions
	paste("P", paste(diffinv(loci[,2]), collapse = " ")),
	indivs.fastphase(fp),
	sep = "\n",
	file = out.fastphase.file)

