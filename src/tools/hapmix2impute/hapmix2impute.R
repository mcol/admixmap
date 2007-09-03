## R script to convert HAPMIXMAP data to IMPUTE format
##
## Input:
## hapmixmap genotypesfile with haploid (training) genotypes
## hapmixmap test genotypesfile with diploid individuals under test
## legend file from HapMap
## recombination rate map file from HapMap (optional)
##
## Output:
## haplotype, genotype, legend and map(if required) file for IMPUTE
##
## WARNING: this script may not be usable for large data sets
## consider also the hapmix2impute program

##TODO: allow filenames as arguments to the script. For now, just edit them below
hapmixmap.train.genotypes.file <- "hapmixdata/genotypes.txt"
hapmixmap.test.genotypes.file <- "hapmixdata/obs_genotypes.txt"
hapmap.legend.file <- "hapmixdata/chr4_legend.txt"
hapmap.map.file <- "hapmixdata/chr4_map.txt"

outputdir <- "."
impute.haplotype.file <- "haplo.txt"
impute.genotype.file <- "geno.txt"
impute.legend.file <- "legend.txt"
impute.map.file <- "map.txt"

num.training.gametes <- 120#number of rows in haploid genotypes file

## convert genotypesfile to haplotype file
##read genotypes file, remove header and first col and change 1 to 0 and 2 to 1
haplotypes <- read.table(hapmixmap.train.genotypes.file, header=T, nrow=num.training.gametes)
first.locus <- dimnames(haplotypes)[[2]][2]
last.locus <- dimnames(haplotypes)[[2]][ncol(haplotypes)]
haplotypes <- t(as.matrix(haplotypes[,-1])-1)
##write transposes the matrix
write(haplotypes, paste(outputdir, impute.haplotype.file, sep="/"), ncol=ncol(haplotypes))
rm(haplotypes)

## process legend file
legend <- read.table(hapmap.legend.file, header=T)
first.index <- match(first.locus, legend[,1])
last.index <- match(last.locus, legend[,1])

write.table(legend[first.index:last.index,],paste(outputdir, impute.legend.file, sep="/"), row.names=F, col.names=T, quote=F)

## process map file
map <- read.table(hapmap.map.file, header=T)
write.table(map[first.index:last.index,], paste(outputdir, impute.map.file, sep="/"), row.names=F, col.names=T, quote=F)
rm(map)

## convert test(diploid) genotypes
legend.data <- data.frame(legend[,-1])
dimnames(legend.data) <- list(legend[,1], dimnames(legend)[[2]][-1])
rm(legend)

genotypes <- read.table(hapmixmap.test.genotypes.file, header=T, colClasses="character")
genotypes[genotypes=="1,1"] <- "1 0 0"
genotypes[genotypes=="1,2"] <- "0 1 0"
genotypes[genotypes=="2,1"] <- "0 1 0"
genotypes[genotypes=="2,2"] <- "0 0 1"
genotypes[genotypes=="0,0"] <- "0 0 0"

typed.loci <- dimnames(genotypes)[[2]][-1]
genotype.data <- data.frame(c(1:(ncol(genotypes)-1)), typed.loci, legend.data[typed.loci,], t(genotypes[,-1]))
write.table(genotype.data, paste(outputdir, impute.genotype.file, sep="/"), col.names=F, row.names=F, quote=F)
rm(genotypes, genotype.data, legend.data)



