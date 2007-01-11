## Functions file for maskGenotypes.R

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
# It assumes that rows have labels stored not as a simple column, but as
# a row names in the manner that read.table handles row.names
# parameter. Hence row.names = TRUE.
genepi.write.table <- function(obj, out.file.name) {
	write.table(
		obj,
		file = out.file.name,
		row.names = TRUE,
		col.names = NA,
		sep = "\t",
		quote = TRUE,
		na = missing.genotype)
}


