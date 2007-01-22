# script to mask genotypes (set them to missing) in order to assess
# prediction of missing genotypes in HAPMIXMAP
#
# This file should be called in batch mode, exactly like this:
#
# R CMD BATCH --no-save --no-restore \
# --population=Eur \
# --basename=genotypes5000 \
# --loci-file=Eur/chr22data/loci5000.txt \
# --percent-indivs=10 \
# --percent-loci=10 \
# --indiv-offset=100 \
# maskGenotypes.R
#
# All the options must be specified and they must be in the same order
# as above.

source("maskGenotypesFunctions.R")

########################################################################
## Configuration
########################################################################

##note:assuming genotypes to be masked are all diploid 
missing.genotype <- "\"0,0\""

# Read command-line arguments. The arguments are recognized by their
# location, but nevertheless their names are being checked.
args <- commandArgs()
population <- get_option("population", args[7])
# Base name for all the input and output files.
io.file.basename <- get_option("basename", args[8])
in.loci.file <- get_option("loci-file", args[9])
percent.missing.indivs <- get_option("percent-indivs", args[10], int = TRUE)
percent.missing.loci <- get_option("percent-loci", args[11], int = TRUE)
indiv.offset <- get_option("indiv-offset", args[12], int = TRUE)


io.file.path <- paste(population, "chr22data", sep = "/")

in.genotypes.file <- get.io.filename(".txt")
# in.loci.file <- paste(io.file.path, "loci5000.txt", sep = "/")
out.genotypes.file <- get.io.filename("_masked.txt")
out.index.file <- get.io.filename("_index.txt")
out.observed.genotypes.file <- get.io.filename("_observed.txt")
out.observed.genotypes.dput.file <- get.io.filename("_observed_dput.txt")
out.fastphase.file <- get.io.filename("_fastphase.inp")

#################################################################

message("Reading data with read.table(). This might take a while.")
message("The file being read is: ", in.genotypes.file)
genotypes <- read.table(
	in.genotypes.file,
	header = TRUE,
	na.strings = c("\"0,0\"", "0,0", "\"0\"", "0"),
	colClasses = "character",
	row.names = "Individ")

# genotypes <- read.table( "Eur/chr22data/d1_test.txt", header = TRUE, na.strings = c("\"0,0\"", "0,0", "\"0\"", "0"), colClasses = "character", row.names = "Individ")

message("Data reading finished.")

number.loci <- ncol(genotypes)
number.indivs <- nrow(genotypes)

size.by.percent <- function(n, pc) {
	if (pc < 0) {
		stop("Percent must be positive")
	}
	if (pc > 100) {
		stop("Percent must be at most 100. Given: ", pc)
	}
	return(round(n * pc / 100.0))
}

missing.loci <- sort(sample(
                   c(1:number.loci),
	           size = size.by.percent(number.loci, percent.missing.loci),
	           replace = FALSE))

missing.indivs <- sort(sample(
	               c(1:number.indivs),
	               size = size.by.percent(number.indivs, percent.missing.indivs),
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

# Apply the offset
missing.indivs = missing.indivs + indiv.offset

cat(
	paste("maskedindivs = ", paste(missing.indivs, collapse = " ")),
	paste("maskedloci = ", paste(missing.loci, collapse = " ")),
	sep = "\n",
	file = out.index.file)

########################################################################
#  Write the data in fastPHASE format
########################################################################

source("FastPhaseConverter.R")

loci = read.table(
	in.loci.file,
	header = TRUE,
	sep = "\t",
	row.names = "SNPid")

message("Saving fastPHASE file to ", out.fastphase.file, ".")
save.fastphase(genotypes, loci, out.fastphase.file)
message("Saving finished.")

