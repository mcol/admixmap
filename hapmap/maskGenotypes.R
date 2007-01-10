##script to mask genotypes (set them to missing) in order to assess prediction of missing genotypes in HAPMIXMAP

# Base name for all the input and output files.

io.file.path <- "Eur/chr22data"
io.file.basename <- "genotypes5000"

# Returns a file name derived from the path and the base name
get.io.filename <- function(extension) {
	return(paste(
		io.file.path,
		"/",
		io.file.basename,
		extension,
		sep = ""))
}

# Saves the table with our custom settings
genepi.write.table <- function(obj, out.file.name) {
	write.table(
		obj,
		file = out.file.name,
		row.names = FALSE,
		col.names = TRUE,
		sep = "\t",
		quote = TRUE,
		na = missing.genotype)
}

in.genotypes.file <- get.io.filename(".txt")
out.genotypes.file <- get.io.filename("_masked.txt")
out.index.file <- get.io.filename("_index.txt")
out.original.genotypes.file <- get.io.filename("_original.txt")
out.changed.genotypes.file <- get.io.filename("_changed.txt")

percent.missing.indivs <- 10
percent.missing.loci <- 10

##note:assuming genotypes to be masked are all diploid 
missing.genotype <- "\"0,0\""
genotypes <- read.table(
	in.genotypes.file,
	header = TRUE,
	na.strings = c("\"0,0\"", "0,0", "\"0\"", "0"),
	colClasses = "character")
genotypes <- data.frame(genotypes)

nLoci <- ncol(genotypes)
nIndivs <- nrow(genotypes)

missing.loci <- sample(
	c(1:nLoci),
	size = round((nLoci * percent.missing.loci / 100.0)),
	replace = FALSE)

missing.indivs <- sample(
	c(1:nIndivs),
	size = round((nIndivs * percent.missing.indivs / 100.0)),
	replace = FALSE)

# Save original genotypes
genepi.write.table(genotypes, out.original.genotypes.file)

# Save changed genotypes, in their original version.
# This file should be small.
genepi.write.table(genotypes[missing.indivs, ], out.changed.genotypes.file)

# Erase the appropriate data by inserting NA values
genotypes[missing.indivs, missing.loci] <- NA

# Save the (big) changed file.
genepi.write.table(genotypes, out.genotypes.file)

cat(
	"maskedindivs = ", sort(missing.indivs), "\n",
	"maskedloci = ", sort(missing.loci), "\n",
	file = out.index.file)
