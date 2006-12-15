##script to mask genotypes (set them to missing) in order to assess prediction of missing genotypes in HAPMIXMAP

in.genotypes.file <- "Eur/chr22data/genotypes_5000.txt"
out.genotypes.file <- "Eur/chr22data/genotypes_5000_masked.txt"
out.index.file <- "Eur/chr22data/masks.txt"

percent.missing.indivs <- 10
percent.missing.loci <- 10

##note:assuming genotypes to be masked are all diploid 
missing.genotype <- "\"0,0\""
genotypes <- read.table(in.genotypes.file, header=T, na.strings=c("\"0,0\"", "0,0", "\"0\"", "0"), colClasses="character")
genotypes <- data.frame(genotypes)

nLoci <- ncol(genotypes)
nIndivs <- nrow(genotypes)

missing.loci <- sample(c(1:nLoci), size=round((nLoci*percent.missing.loci / 100)), replace=F)
missing.indivs <- sample(c(1:nIndivs), size = round((nIndivs*percent.missing.indivs / 100)), replace=F)

genotypes[missing.indivs, missing.loci] <- missing.genotype

write.table(genotypes, file=out.genotypes.file, row.names=F, col.names=T)

cat("maskedindivs = ", sort(missing.indivs), "\nmaskedloci = ", sort(missing.loci), "\n", file=out.index.file)
