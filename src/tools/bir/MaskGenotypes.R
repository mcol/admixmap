## function to mask genotypes and produce these outputs:
## 1. haploid training data for HAPMIXMAP with monomorphic loci removed
## 2. diploid data at a subset of the training loci (ie omitting the masked genotypes)
## 3. locus file with the same monomorphic loci removed
## 4. observed values of the masked genotypes, coded as 1, 2, 3
## 5. allele counts over the unmasked individuals at the masked loci
##

MaskGenotypes <- function(percent.masked.loci, num.masked.indivs, prefix,
                          data.dir, num.loci, num.gametes, fastphase.diploid.file, fastphase.haploid.file){
  ## percent.masked.loci
  ## num.masked.indivs = number of individuals to mask
  ## prefix = prefix of original genotypesfile and locusfile
  ## data.dir = directory to write data
  ## num.loci = number of loci in data (for convenience)
  ## num.gametes = number of HapMap gametes (for convenience)
  ## fastphase.diploid.file, fastphase/haploid.file = names of fastphase data files to write
  
  ##determine how many individuals are to be masked and at which loci 
  num.masked.loci <- num.loci*percent.masked.loci/100
  num.unmasked.gametes <- num.gametes - (num.masked.indivs*2)
  masked.loci.indices <- sort(sample(num.loci, size=num.masked.loci, replace=F))

  ## write a list of the masked loci to file
  cat(masked.loci.indices, file=paste(data.dir, "masked_loci.txt", sep="/"))
  
  ## read in the original phased genotypes
  all.geno <- read.table(paste(prefix,"_genotypes.txt", sep=""),
                         header=T, nrows=num.gametes, colClasses=c("character", rep("integer", num.loci)))
  masked.loci <- dimnames(all.geno)[[2]][masked.loci.indices+1]##+1 to skip ID column  
  
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
  if(length(monomorphic.masked.loci.indices)>0){
    dimorphic.masked.loci <- masked.loci[-monomorphic.masked.loci.indices]

    ##revise unmasked indivs at dimorphic masked loci
    unmasked.geno <- unmasked.geno[,-monomorphic.masked.loci.indices]
  
  }else{##no monomorphic masked loci
    dimorphic.masked.loci <- masked.loci
  }
  num.dimorphic.masked.loci <- length(dimorphic.masked.loci)
  
  col.sums <- colSums(unmasked.geno.all.loci)
  monomorphic.loci.indices <- which((( col.sums == num.unmasked.gametes) | (col.sums == 2*num.unmasked.gametes)), all.loci)
  if(length(monomorphic.loci.indices)>0){
    all.dimorphic.loci <- all.loci[-monomorphic.loci.indices]

  }else{#no monomorphic loci
    all.dimorphic.loci <- all.loci
  }
  dimorphic.unmasked.loci <- all.dimorphic.loci[-match(dimorphic.masked.loci, all.dimorphic.loci)]
  rm(col.sums)

  ##TODO: should stop if all loci are monomorphic
  
  ##masked indivs at dimorphic masked loci
  masked.geno <- all.geno[(num.unmasked.gametes + 1):num.gametes, dimorphic.masked.loci]
  ##masked gametes at unmasked loci (include ID column)
  masked.haploid.geno <- all.geno[(num.unmasked.gametes + 1):num.gametes, dimorphic.unmasked.loci]
  
  #########################################
  ## write genotypes for unmasked individuals at all dimorphic loci back to file
  #########################################
  unmasked.geno.all.loci <- data.frame(ID=all.geno[1:num.unmasked.gametes,1], all.geno[1:num.unmasked.gametes, all.dimorphic.loci])
  write.table(unmasked.geno.all.loci, file=paste(data.dir, "train_genotypes.txt", sep="/"), row.names=F, col.names=T)

  ######################################
  ## recode observed diploid masked genotypes as integers 1, 2,3
  ## and save to file
  ######################################
  obs.masked.geno <- NULL
  for(i in 1:num.masked.indivs){
    obs.masked.geno <- rbind(obs.masked.geno, masked.geno[(i*2 -1),]+masked.geno[(i*2),]-1)
  }
  dput(obs.masked.geno, paste(data.dir,"obs_masked_genotypes.txt", sep="/"))
  rm(obs.masked.geno)
  
  ##################################
  ## write unmasked diploid genotypes of masked individuals to file in HAPMIXMAP format,
  ## coded as pairs of integers
  ##################################
  masked.diploid.geno <- NULL
  for(i in 1:num.masked.indivs){
    masked.diploid.geno <- rbind(masked.diploid.geno, paste(masked.haploid.geno[(i*2-1),], masked.haploid.geno[(i*2),], sep=","))
  }
  masked.diploid.geno <- data.frame(ID=all.geno[seq(num.unmasked.gametes + 1, num.gametes, 2),1], masked.diploid.geno)
  dimnames(masked.diploid.geno)[[2]] <- c("ID", dimorphic.unmasked.loci)
  write.table(masked.diploid.geno, paste(data.dir, "test_genotypes.txt", sep="/"), row.names=F, col.names=T)
  rm(masked.diploid.geno)
  
  ###################################
  ## compute genotype frequencies over unmasked individuals at masked loci
  ####################################
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
  
  #################################
  ## read locusfile and exclude monomorphic loci
  #################################
  locusfile <- read.table(paste(prefix, "_loci.txt", sep=""), na.strings=c("NA, #"),
                          comment.char="%", header=T, colClasses="character", nrow=num.loci)
  locusfile[1,3] <- "0"
  if(length(monomorphic.loci.indices)>0){
    locusfile.no.mm <- locusfile[-monomorphic.loci.indices,]
    pos <- cumsum(locusfile[,3])[-monomorphic.loci.indices]
  }else{
    locusfile.no.mm <- locusfile
    pos <- cumsum(locusfile[,3])
  }
  pos2 <- c(0, pos[1:(length(pos)-1)])
  

  ##replace distances
  locusfile.no.mm[,3] <- format(pos-pos2)
  rm(pos2)
  locusfile.no.mm[1,3] <- "#"
  write.table(locusfile.no.mm, paste(data.dir, "train_loci.txt", sep="/"), row.names=F, col.names=T, quote=F)

  ##############################
  ## write fastPHASE input files
  ###############################

  ## write masked diploid genotypes
  if(!is.null(fastphase.diploid.file)){
    write(num.masked.indivs, file=fastphase.diploid.file)
    write(length(all.dimorphic.loci), file=fastphase.diploid.file, append=T)
    cat("P ", as.numeric(pos)*1e6, "\n", file=fastphase.diploid.file, append=T)
    ##  write("\nSSSSS", file=fastphase.diploid.file, append=T)
    ## NOTE: S line seems to confuse fastPHASE when P line is there ??

    ##convert to character an dreplace missing values
    masked.geno.all.loci <- format(all.geno[(num.unmasked.gametes+1):num.gametes, all.dimorphic.loci]-1)
    masked.geno.all.loci[,dimorphic.masked.loci] <- "?"##missing genotypes are denoted by '?'
    
    for(i in 1:num.masked.indivs){
      ## id line
      cat( "# ", all.geno[(num.unmasked.gametes+i*2-1),1], "\n", file=fastphase.diploid.file, append=T)
      ## 2 lines per individual
      cat(unlist(masked.geno.all.loci[(i*2-1),]), "\n", file=fastphase.diploid.file, append=T)
      cat(unlist(masked.geno.all.loci[(i*2),]), "\n", file=fastphase.diploid.file, append=T)
    }
  }
  if(!is.null(fastphase.haploid.file)){
    ## write unmasked gametes to haplotype file
    haploid.out <- t(cbind(all.geno[1:num.unmasked.gametes,1], matrix(apply(unmasked.geno.all.loci[,-1], 1, paste, collapse=" "),
                                                         nrow=num.unmasked.gametes, byrow=T)))
    cat(num.unmasked.gametes, haploid.out, file=fastphase.haploid.file, sep="\n")    
  }
}

