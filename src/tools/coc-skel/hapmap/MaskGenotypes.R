## script to obtain correct genotype freqs and observed genotypes for masked individuals
## 

data.dir <- "/ichec/work/ndlif006b/maciej/chr22-data/Eur/chr22data"
num.loci <- 34026
num.gametes <- 120

##determine how many individuals were masked and at which loci 
results.dir <- "/ichec/work/ndlif006b/maciej/irw-2-arp-12-0.5-2-afpp-2-8-s1/hapmap/Eur/Chr22results6States2"
PPG <- dget(paste(results.dir, "PPGenotypeProbs.txt", sep="/"))
masked.loci <- dimnames(PPG)[[2]]
num.masked.indivs <- dim(PPG)[3]
num.masked.loci <- length(masked.loci)
rm(PPG)

## read in the original phased genotypes
all.geno <- read.table(paste(data.dir,"phased_genotypes.txt", sep="/"), header=T, nrows=num.gametes, colClasses=c("character", rep("integer", num.loci)))

##separate into masked and unmasked individuals and exclude unmasked loci
num.unmasked.gametes <- num.gametes - (num.masked.indivs*2)
masked.geno <- all.geno[(num.unmasked.gametes + 1):num.gametes, masked.loci]
unmasked.geno <- all.geno[1:num.unmasked.gametes,masked.loci]
rm(all.geno)


## recode observed diploid masked genotypes as integers 1, 2,3
## and save to file
masked.diploid.geno <- NULL
for(i in 1:num.nmasked.indivs){
  masked.diploid.geno <- rbind(masked.diploid.geno,masked.geno[(i*2 -1),]+masked.geno[(i*2),]-1)
}

dput(masked.diploid.geno, paste(data.dir,"obs_masked_genotypes.txt", sep="/"))


##compute genotype frequencies over unmasked individuals at masked loci
##first, calculate frequency of allele 1, p
p <- colSums(unmasked.geno==1)/num.unmasked.gametes
gfreqs <- cbind(p*p, 2*p*(1-p), (1-p)*(1-p) )

dimnames(gfreqs) <- list(dimnames(unmasked.geno)[[2]], c("1", "2", "3"))
dput(gfreqs, paste(data.dir,"genotype_freqs.txt", sep="/"))


