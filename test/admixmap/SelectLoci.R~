indexCompoundLoci <- function(v.Chr, v.x, c2) {
  ## indexes compound loci defined by minimum separation distance from other compound loci
  ## positions must be in ascending order within each chromosome
  chr <- as.factor(v.Chr)
  chr.levels <- levels(chr)
  num.chr <- length(chr.levels)
  clocus <- integer(0)
  for(lev in 1:num.chr) {
    x <- v.x[chr==chr.levels[lev]]
    new.clocus <- c(T, x[-1] > x[-length(x)] + c2)
    clocus <- c(clocus, as.integer(new.clocus))
  }
  return(cumsum(clocus))  
}

lociToDrop <- function(clocus, x, c2) {
  ## within each compound locus that spans >= c2, drop loci until
  ## an interval between flanking loci > c2 
  drop <- rep(F, length(clocus))
  num.cloci <- length(levels(as.factor(clocus)))
  seq.cloci <- 1:num.cloci
  first <- match(seq.cloci, clocus)  # indexes first snp in each clocus
  last <- 1 + length(clocus) - match(seq.cloci, rev(clocus)) # indexes last snp in each clocus
  clocus.span <- x[last] - x[first]  # vector of distances spanned by each compound locus
  for(j in seq.cloci) {
    if(clocus.span[j] > c2) {
      a <- first[j]
      b <- last[j]
      m <- length(a:b)
      kept <- rep(T, m)
      z <- x[a:b][kept]
      intervals <- z[-c(1,2)] - z[-c(m-1, m)]
      while(max(intervals) <= c2) {
        z <- x[a:b][kept]
        m <- length(z)
        intervals <- z[-c(1,2)] - z[-c(m-1, m)]
        kept[kept][1 + which.max(intervals)] <- F # drop locus in largest interval
      }
      drop[a:b] <- !kept
    }
  }
  return(drop)
}     

lociToKeep <- function(clocus, x, nS) {
  ## keeps up to n sloci in each compound locus
  num.cloci <- length(levels(as.factor(clocus)))
  seq.cloci <- 1:num.cloci
  first <- match(seq.cloci, clocus)  # indexes first snp in each clocus
  last <- 1 + length(clocus) - match(seq.cloci, rev(clocus)) # indexes last snp in each clocus
  keep <- rep(T, length(clocus))
  clocus.levels <- levels(as.factor(clocus))
  new.clocus <- match(clocus.levels, clocus)
  for(j in seq.cloci) {
    ## get number of sloci
    a <- first[j]
    b <- last[j]
    m <- b - a + 1
    while(m > nS) {
      d.first <- x[a+1] - x[a]
      d.last <- x[b] - x[b-1]
      if(d.first > d.last) {
        keep[a] <- F
        a <- a + 1
      } else {
        keep[b] <- F
        b <- b - 1
      }
      m <- b - a + 1
    }
  }
  return(keep)  
}

distanceFromLast <- function(v.Chr, v.x) {
## given vector of chr labels and positions x 
## this function returns vector of distances from last locus
  d <- numeric(0)
  chr <- as.factor(v.Chr)
  chr.levels <- levels(chr)
  num.chr <- length(chr.levels)
  for(lev in 1:num.chr) {
   x <- v.x[chr==chr.levels[lev]]
   x <- x[order(x)]
   d <- c(d, NA, c(x[-1]) - x[-length(x)])
 }
  return(d)
}

# call with arguments such as 
# selectLoci("data/SNPPositionsAll.txt", 0.4, 4)
#
#
#
selectLoci <- function(snpfile, c2, nS) {
  SNPs <- read.table(file=snpfile, sep="\t", as.is=T,
                     colClasses=c("character", "factor", "numeric", "numeric"),
                     comment.char="", na.strings="NA", header=F) 
  dimnames(SNPs)[[2]] <- c("rsnumber", "chr", "bp", "cM")
  SNPs$rsnumber <- paste("rs", SNPs$rsnumber, sep="")
  
  ## drop rows with missing position from SNPs
  SNPs <- SNPs[!is.na(SNPs$cM), ]
  ## drop duplicates from SNPs
  SNPs <- SNPs[!duplicated(SNPs$rsnumber), ] 
  
  ## sort remaining loci by map position 
  SNPs.kept <- SNPs[order(1000*as.numeric(SNPs$chr)+SNPs$cM), ]
  ## index compound loci
  clocus <- indexCompoundLoci(SNPs.kept$chr, SNPs.kept$cM, c2)
  cat("\n", dim(SNPs.kept)[1], "SNPs grouped into", length(levels(as.factor(clocus))),
      "compound loci\n")
  
  i <- 1
  while(i < 4) {## drop SNPs so as to split any compound loci that span > c2
    drop <- lociToDrop(clocus, SNPs.kept$cM, c2)
    SNPs.kept <- SNPs.kept[!drop, ]
    clocus <- indexCompoundLoci(SNPs.kept$chr, SNPs.kept$cM, c2)
    cat("\n", dim(SNPs.kept)[1], "SNPs grouped into", length(levels(as.factor(clocus))),
        "compound loci after iteration", i, "of dropping SNPs to split compound loci\n")
    i <- i + 1
  }
  
  ## keep nS SNPs in each compound locus
  ## chop off snps from either end until num snps < nS
  keep <- lociToKeep(clocus, SNPs.kept$cM, nS)
  SNPs.kept <- SNPs.kept[keep, ]
  clocus <- indexCompoundLoci(SNPs.kept$chr, SNPs.kept$cM, c2)
  cat("\n", dim(SNPs.kept)[1], "SNPs grouped into", length(levels(as.factor(clocus))),
      "compound loci after keeping only", nS, "SNPs in each compound locus\n")
  
  ## recalculate distances
  cMFromLast <- distanceFromLast(SNPs.kept$chr, SNPs.kept$cM)
  cMFromLast[cMFromLast < c2] <- 0
  bpFromLast <- distanceFromLast(SNPs.kept$chr, SNPs.kept$bp)
  LocusNames <- SNPs.kept$rsnumber
  
  ## write locus table
  d <- 0.01 * cMFromLast
  loci <- data.frame(LocusNames, rep(2, length(LocusNames)), d, SNPs.kept$chr)
  dimnames(loci)[[2]] <- c("LocusName", "NumAlleles", "DistanceFromLast", "Chr")
  write.table(loci, file=paste("loci_", snpfile, sep=""), sep="\t", row.names=F)
}

