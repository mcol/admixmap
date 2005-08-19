simulateHaploidAlleles <- function(M,ro,d,L) {
#  T <- t(matrix(data=c(M + (1-M)*exp(-ro*x),  1-M - (1-M)*exp(-ro*x),
#                M - M*exp(-ro*x),      1-M + M*exp(-ro*x)    ),
#              nrow=2, ncol=2))
  gameteAncestry <- numeric(L)
  randAnc <- runif(L)
  x <- numeric(L)
  gameteAncestry[1] <- ifelse(randAnc[1] > M, 2, 1) # M is prob pop 1.
  print('new gamete')
  print(gameteAncestry[1])
  x[1] <- ifelse(runif(1) > alleleFreqs[1,gameteAncestry[1]], 2, 1)
  for(locus in 2:L) {
    T <- t(matrix(data=c(M + (1-M)*exp(-ro*d[locus-1]), 1-M - (1-M)*exp(-ro*d[locus-1]),
                    M - M*exp(-ro*d[locus-1]), 1-M + M*exp(-ro*d[locus-1])    ),
                  nrow=2, ncol=2))
    gameteAncestry[locus] <- ifelse(randAnc[locus] > T[gameteAncestry[locus-1], 1],
                                    2, 1)
    x[locus] <- ifelse(runif(1) > alleleFreqs[1,gameteAncestry[locus]], 2, 1)
    if( d[locus-1] > 99.0 ) print(gameteAncestry[locus])
  }
  return(x)
}

simulateDiploidGenotypes <- function(Mpat,Mmat,Rhopat,Rhomat,x,L) {
  paternalGamete <- simulateHaploidAlleles(Mpat,Rhopat,x,L)
  maternalGamete <- simulateHaploidAlleles(Mmat,Rhomat,x,L)
  simulateGenotypes <- paste(paternalGamete, ",", maternalGamete, sep="")
}

simulateHaploidGenotypes <- function(M,ro,x,L) {
  paternalGamete <- simulateHaploidAlleles(M,ro,x,L)
  simulateGenotypes <- paste(paternalGamete, ",", paternalGamete, sep="")
}

#ro <- 5
#L <- 30
#x <- runif(L-1)*0.1
#x <- rep(0.01,L-1)
#x <- c(rep(0.1,L-1),rep(c(100,rep(0.1,L)),21))
#xx<-read.table('~/admix/minstead/test/data/loci_frudakis.txt',header=T)
#x<-xx[2:length(xx[,1]),3]

chr.L<-c(292,272,233,212,197,201,184,166,166,181,156,169,117,128,110,130,128,123,109,96,59,58)
x<-c(rep(0.05,chr.L[1]/5))
yy<-seq(from=1,to=length(x)+1,by=3)
last<-length(x)+1
for( i in 2:22){
  chrn<-c(100,rep(0.05,chr.L[i]/5))
  x<-c(x,chrn)
  y<-seq(from=1,to=length(chrn),by=3)
  yy<-c(yy,y+last)
  last<-last+length(chrn)
}
L <- length(x)
L1 <- length(x)+1

I <- 9
alleleFreqs <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow=2, ncol=2)
#alleleFreqs <- matrix(c(1.0, 0.0, 0.0, 1.0), nrow=2, ncol=2)
#freqs<-dget('~/admix/minstead/frudakis/data/allelefreqs_RObj.txt')
#alleleFreqs <- matrix(nrow=length(freqs[2,]), ncol=2)
#alleleFreqs[,1] <- as.numeric(freqs[2,])
#alleleFreqs[,2] <- as.numeric(freqs[3,])
genotypes <- character(L1)
for(individual in 1:I){
  Mpat <- 0.5
  Mmat <- 0.8
  Rhopat <- 1
  Rhomat <- 5
#  obs <- simulateHaploidGenotypes(M, ro, x, L)
  obs <- simulateDiploidGenotypes(Mpat, Mmat, Rhopat, Rhomat, x, L1)
  genotypes <- rbind(genotypes, obs)
}
genotypes <- genotypes[-1,]
id=seq(1:I)
genotypes <- data.frame(as.character(id), rep("0",I), genotypes, row.names=NULL)
write.table(genotypes, file="simdata/genotypes_sim.txt", row.names=FALSE)
write.table(genotypes[,c(1,2,yy+2)],
            file="simdata/genotypes_sim15.txt", row.names=FALSE)

loci <- data.frame(as.vector(dimnames(genotypes)[[2]][-c(1,2)]),
                   rep(2,L1),
#                   c(100, rep(x,L-1)),
                   c(100, x),
                   rep("1",L1),
                   row.names=NULL)
dimnames(loci)[[2]] <- c("Locus", "NumAlleles", "Distance", "Chromosome")
write.table(loci, file="simdata/loci_new.txt", row.names=FALSE)
write.table(loci[yy,], file="simdata/loci_sim15.txt", row.names=FALSE)

allelefreqs <- data.frame(as.vector(loci[,1]),
                          rep(alleleFreqs[1,1],L1),
                          rep(alleleFreqs[1,2],L1))
dimnames(allelefreqs)[[2]] <- c("Locus", "Pop1", "Pop2")
write.table(allelefreqs, file="simdata/allelefreqfile.txt", row.names=FALSE)

priorallelefreqs <- numeric(2)
locusnames <- character(0)
for(locus in 1:L1) {
  priorallelefreqs <- rbind(priorallelefreqs, 500*alleleFreqs + 0.5)
  locusnames <- c(locusnames, rep(as.vector(loci[locus,1]), 2))
}
priorallelefreqs <- data.frame(locusnames, priorallelefreqs[-1,], row.names=NULL)
write.table(priorallelefreqs, file="simdata/priorallelefreqs_sim.txt", row.names=FALSE)
write.table(priorallelefreqs[sort(c(2*yy,2*yy-1)),], file="simdata/priorallelefreqs_sim15.txt", row.names=FALSE)
                               
