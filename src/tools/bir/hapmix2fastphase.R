## author: Maciej Blizinski

## Converts data in hapmapmix format into fastPHASE format.
##
## Data shoud be given as two objects:
##
## 1. genotypes, a matrix where rows are individuals and columns are
## diploid genotypes defined by values "1,1", "1,2", "2,2"
## 
## 2. loci, a matrix where rows are loci and the second column stores
## distance between current and previous loci.
##
## Typical call: save.fastphase(genotypes, loci, "fastphase.inp")

## Decode diploid data, two functions
## Return first genotype
##
## It could be also done with substr(x, 1, 1), but this seemingly stupid
## way is actually good for returning NA values. It also assures that
## only "1,1", "1,2", "2,1" and "2,2" values are recognized and warnings
## are provided if a different value is spotted.

get.fastphase.diploid.genotype <- function(x, n) {
  if (!(n == 1 || n == 2)) {
    stop("Accepted second arguments are 1 and 2.")
  }
  if (is.na(x)) { return(NA) }
  else if (x == "1,1") { return("0") }
  else if (x == "1,2") { return(format(c(1,2)[n]-1)) }
  else if (x == "2,1") { return(format(c(2,1)[n]-1)) }
  else if (x == "2,2") { return("1") }
  else if (x == "0,0") { return("?") }
  else {
    stop("Unrecognized value: ", x)
    return(NA)
	}
}

write.fastphase.haplotypefile <- function(haploid.genotypes, fastphase.file){
  ## takes HAPMIXMAP haploid genotypes (haplotypes) and locusfile
  ## and write a fastPHASE format haplotypefile, supplied with -b option
  ## haploid.genotypes is an N*(1+L) array of 1s and 2s, with Ids in the first column
  
  haploid.out <- t(cbind(haploid.genotypes[,1], matrix(apply(haploid.genotypes[,-1], 1, paste, collapse=" "), nrow=nrow(haploid.genotypes), byrow=T)))
  cat(nrow(haploid.genotypes), haploid.out, file=fastphase.file, sep="\n")
  
}

fastphase.diploid <- function(diploid.genotypes, loci, fastphase.file) {
  ## takes HAPMIXMAP-format diploid genotypes and a locusfile and writes a fastPHASE format input file
  ## there may be more loci in the locusfile than the genotypesfile
  
  number.diploid.indivs <- nrow(diploid.genotypes)
  number.loci <- nrow(loci)
  number.diploid.loci <- ncol(diploid.genotypes)-1
  diploid.loci <- dimnames(diploid.genotypes)[[2]][-1]
  
  ## Three-dimensional array, dimensions:
  ## 1. Individual
  ## 2. Locus
  ## 3. Genotype, for two fastPHASE data file lines
  diploid.out <- array("?", dim = c(number.diploid.indivs, number.loci, 2), dimnames=list(character(0), loci[,1], character(0)))
  
  ## Map the diploid data 
  diploid.out[,diploid.loci , 1] <- sapply(as.matrix(diploid.genotypes[,-1]), get.fastphase.diploid.genotype, 1)
  diploid.out[,diploid.loci , 2] <- sapply(as.matrix(diploid.genotypes[,-1]), get.fastphase.diploid.genotype, 2)


  ##calculate positions in bp
  pos <- cumsum(c(0,loci[-1,3]))*1e6
 
  cat(number.diploid.indivs, number.loci, sep="\n", file = fastphase.file)
  cat("P", pos, "\n", file=fastphase.file, append=T)
  for (indiv in 1:number.diploid.indivs) {
    ## message("Individual: ", indiv)
    write(paste("# ", diploid.genotypes[indiv,1]), file = fastphase.file, append=T)
    write(diploid.out[indiv,,1], file = fastphase.file, append=T, ncol=number.loci)
    write(diploid.out[indiv,,2], file = fastphase.file, append=T, ncol=number.loci)
  }

}

## example
##locusfile <- "chr22data/train_loci.txt"
##diploid.genotypesfile <- "chr22data/test_genotypes.txt"
##haploid.genotypesfile <- "chr22data/train_genotypes.txt"

##loci <- read.table(locusfile, na.strings=c("NA, #"), comment.char="%", header=T, colClasses="character")
##haploid.genotypes <- read.table(haploid.genotypesfile, header=T, colClasses=c("character", rep("integer", nrow(loci))))
##stopifnot(nrow(loci) == ncol(haploid.genotypes)-1)

##write.fastphase.haplotypefile(haploid.genotypes, "chr22data/Rfastphase5khaps.inp")
##rm(haploid.genotypes)

##diploid.genotypes <- read.table(diploid.genotypesfile, header=T, colClasses="character")
##rm(diploid.genotypes)
##rm(loci)
