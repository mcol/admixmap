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
out.genotypes.file <- get.io.filename("_masked.txt")
out.index.file <- get.io.filename("_index.txt")
out.observed.genotypes.file <- get.io.filename("_observed.txt")
out.observed.genotypes.dput.file <- get.io.filename("_observed_dput.txt")

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
	"maskedindivs = ", missing.indivs, "\n",
	"maskedloci = ", missing.loci, "\n",
	file = out.index.file)
