## script to mask genotypes and produce these outputs:
## 1. haploid training data for HAPMIXMAP with monomorphic loci removed
## 2. diploid data at a subset of the training loci (ie omitting the masked genotypes)
## 3. locus file with the same monomorphic loci removed
## 4. observed values of the masked genotypes, coded as 1, 2, 3
## 5. allele counts over the unmasked individuals at the masked loci
##

MaskGenotypes <- function(results.dir, data.dir, num.loci, num.gametes, fastphase.file){
  ##determine how many individuals were masked and at which loci 
  ## normally we determine which individuals' genotypes are to be masked at at which loci some other way (eg random selection).
  ## For now, determine this from the program output.
  PPG <- dget(paste(results.dir, "PPGenotypeProbs.txt", sep="/"))
  masked.loci <- dimnames(PPG)[[2]]
  num.masked.indivs <- dim(PPG)[3]
  num.masked.loci <- length(masked.loci)
  num.unmasked.gametes <- num.gametes - (num.masked.indivs*2)
  rm(PPG)
  
  ## read in the original phased genotypes
  all.geno <- read.table(paste(data.dir,"phased_5000_genotypes.txt", sep="/"), header=T, nrows=num.gametes, colClasses=c("character", rep("integer", num.loci)))
  
  ##separate into masked and unmasked individuals and exclude unmasked loci
  
  ##unmasked gametes at masked loci
  unmasked.geno <- all.geno[1:num.unmasked.gametes,masked.loci]
  ##unmasked gametes at all loci
  unmasked.geno.all.loci <- all.geno[1:num.unmasked.gametes, -1]
  
  all.loci <- dimnames(all.geno)[[2]][-1]

  ##determine which loci are now monomorphic in the unmasked individuals
  ## and exclude monomorphic loci from separated tables
  col.sums <- colSums(unmasked.geno)
  monomorphic.masked.loci.indices <- which((( col.sums == num.unmasked.gametes) | (col.sums == 2*num.unmasked.gametes)), masked.loci)
  dimorphic.masked.loci <- masked.loci[-monomorphic.masked.loci.indices]
  num.masked.loci <- length(dimorphic.masked.loci)
  rm(col.sums)
  
  col.sums <- colSums(unmasked.geno.all.loci)
  monomorphic.loci.indices <- which((( col.sums == num.unmasked.gametes) | (col.sums == 2*num.unmasked.gametes)), all.loci)
  all.dimorphic.loci <- all.loci[-monomorphic.loci.indices]
  dimorphic.unmasked.loci <- all.dimorphic.loci[-match(dimorphic.masked.loci, all.dimorphic.loci)]
  rm(col.sums)
  
  ##masked indivs at dimorphic masked loci
  masked.geno <- all.geno[(num.unmasked.gametes + 1):num.gametes, dimorphic.masked.loci]
  ##masked gametes at unmasked loci (include ID column)
  masked.haploid.geno <- all.geno[(num.unmasked.gametes + 1):num.gametes, dimorphic.unmasked.loci]
  
  ##revise unmasked indivs at dimorphic masked loci
  unmasked.geno <- unmasked.geno[,-monomorphic.masked.loci.indices]
  
  ## revise unmasked indivs at all loci
  unmasked.geno.all.loci <- data.frame(ID=all.geno[1:num.unmasked.gametes,1], all.geno[1:num.unmasked.gametes, all.dimorphic.loci])
  ## write genotypes for unmasked individuals at all dimorphic loci back to file
  write.table(unmasked.geno.all.loci, paste(data.dir, "train_genotypes.txt", sep="/"), row.names=F, col.names=T)
  ##rm(all.geno)
  
  
  ## recode observed diploid masked genotypes as integers 1, 2,3
  ## and save to file
  obs.masked.geno <- NULL
  for(i in 1:num.masked.indivs){
    obs.masked.geno <- rbind(obs.masked.geno, masked.geno[(i*2 -1),]+masked.geno[(i*2),]-1)
  }
  dput(obs.masked.geno, paste(data.dir,"obs_masked_genotypes.txt", sep="/"))
  rm(obs.masked.geno)
  
  ## write unmasked diploid genotypes of masked individuals to file in HAPMIXMAP format, coded as pairs of integers
  masked.diploid.geno <- NULL
  for(i in 1:num.masked.indivs){
    masked.diploid.geno <- rbind(masked.diploid.geno, paste(masked.haploid.geno[(i*2-1),], masked.haploid.geno[(i*2),], sep=","))
  }
  masked.diploid.geno <- data.frame(ID=all.geno[seq(num.unmasked.gametes + 1, num.gametes, 2),1], masked.diploid.geno)
  dimnames(masked.diploid.geno)[[2]] <- c("ID", dimorphic.unmasked.loci)
  write.table(masked.diploid.geno, paste(data.dir, "test_genotypes.txt", sep="/"), row.names=F, col.names=T)
  rm(masked.diploid.geno)
  
  ##compute genotype frequencies over unmasked individuals at masked loci
  ##first, calculate frequency of allele 1, p
  count.allele1 <- colSums(unmasked.geno==1)
  count.allele2 <- colSums(unmasked.geno==2)
  
  ##p <- count.allele1/num.unmasked.gametes
  ##gfreqs <- cbind(p*p, 2*p*(1-p), (1-p)*(1-p) )
  ##dimnames(gfreqs) <- list(dimnames(unmasked.geno)[[2]], c("1", "2", "3"))
  ##dput(gfreqs, paste(data.dir,"genotype_freqs.txt", sep="/"))
  
  gcounts <- rbind(count.allele1, count.allele2)
  dimnames(gcounts) <- list(c("allele1", "allele2"), dimnames(unmasked.geno)[[2]])
  dput(gcounts, paste(data.dir,"obs_allele_counts.txt", sep="/"))
  
  ## read locusfile and exclude monomorphic loci
  locusfile <- read.table(paste(data.dir, "phased_5000_loci.txt", sep="/"), header=T, colClasses="character", nrow=num.loci)
  pos <- cumsum(locusfile[,3])[-monomorphic.loci.indices]
  pos2 <- c(0, pos[1:(length(pos)-1)])
  
  locusfile.no.mm <- locusfile[-monomorphic.loci.indices,]
  ##replace distances
  locusfile.no.mm[,3] <- pos-pos2
  rm(pos2)
  write.table(locusfile.no.mm, paste(data.dir, "train_loci.txt", sep="/"), row.names=F, col.names=T, quote=F)
  
  ## write fastPHASE input file
  write(num.gametes/2, file=fastphase.file)
  write(length(all.dimorphic.loci), file=fastphase.file, append=T)
  cat("P ", as.numeric(pos)*1e6, file=fastphase.file, append=T)
  write("\nSSSSS", file=fastphase.file, append=T)
  
  ## write unmasked gametes
  num.dimorphic.loci <- length(all.dimorphic.loci)
  for(i in 1:(num.unmasked.gametes/2)){
    ## id line
    cat( "# ", all.geno[(i*2-1),1], "\n", file=fastphase.file, append=T)
    ## 2 lines per individual
    cat(as.numeric(c(unmasked.geno.all.loci[(i*2-1),-1]))-1, "\n", file=fastphase.file, append=T)
    cat(as.numeric(c(unmasked.geno.all.loci[(i*2),-1]))-1, "\n", file=fastphase.file, append=T)
  }
  
  ## write masked gametes
  masked.geno.all.loci <- format(all.geno[(num.unmasked.gametes+1):num.gametes, all.dimorphic.loci]-1)
  masked.geno.all.loci[,dimorphic.masked.loci] <- "?"
  
  for(i in 1:num.masked.indivs){
    ## id line
    cat( "# ", all.geno[(num.unmasked.gametes+i*2-1),1], "\n", file=fastphase.file, append=T)
    ## 2 lines per individual
    cat(unlist(masked.geno.all.loci[(i*2-1),]), "\n", file=fastphase.file, append=T)
    cat(unlist(masked.geno.all.loci[(i*2),]), "\n", file=fastphase.file, append=T)
  }

}


## loop over populations
Pops <- c("Eur", "Afr", "Asian")
Panels <- c("CEU", "YRI", "JPTCHB")
num.loci <- c(34026,38224, 32713)
num.gametes <- c(120, 120, 180)

for(pop in 1:3){
  results.dir <- paste("/ichec/work/ndlif006b/maciej/irw-2-arp-12-0.5-2-afpp-2-8-s1/hapmap", Pops[pop], sep="/")
  results.dir <- paste(results.dir, "Chr22Results6States2", sep="/")
  
  data.dir <- paste("/ichec/work/ndlif006b/genepi/hapmap/data/chr22/hapmixmap", Panels[pop], sep="/")

  ## file to put fastphase data
  fastphase.file <- paste("/ichec/work/ndlif006b/genepi/hapmap/data/chr22/fastphase", Panels[pop], sep="/")

  MaskGenotypes(results.dir, data.dir, num.loci[pop], num.gametes[pop], fastphase.file)
}
