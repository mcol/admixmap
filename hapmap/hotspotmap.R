## Script to plot posterior means of lambda parameters (expected numbers of arrivals between loci) in hapmix model along with recombination rates
## from HapMap
##currently only does one chromosome
##TODO: make a collection of plots for the whole genome; also plot subsets of chromosomes
## NOTE: the loci used in ADMIXMAP and the map may be different. This is accounted for in the script

lambdafile <- "Eur/Results6States/lambdaPosteriorMeans.txt"##file with lambda posterior means
locusfile  <- "Eur/chr22data/loci.txt"                     ##locus file as used with ADMIXMAP, for distances
ratefile   <- "Eur/chr22data/genetic_map_chr22.txt"        ##recombination rates from HapMap
outfile <- "HotspotMap.ps"

##note: there may be different loci in the locusfile (and therefore used in the program) than in the ratesfile
##so we must use care when combining the two plots

##read posterior means
lambda<-read.table(lambdafile, header=F)
lambda<-as.vector(lambda[,1])
##plot arrival numbers as time series
#plot(1:22985, lambda, type='l')

rates<-read.table(ratefile, header=T)[,1:2]
##rescale positions to start at 0 and measure in Mb (HapMap files measure in bp)
pos0<-rates[1,1]
rates[,1]<-rates[,1]-pos0
rates[,1]<-rates[,1]/1e06

##read distances from locusfile
distances<-read.table(locusfile, header=T, comment.char="", na.strings="#")[,3]
##convert distances in locusfile into positions in Mb
positions <- c(0, cumsum(distances[-1]))

##compute distances from positions in ratefile
#d2<-numeric(nrow(rates))
#d2[1]<-0
#for(i in 2:nrow(rates)){
#  d2[i]<-rates[i,2]-rates[i-1,2]
#}

max.posn <- max(positions[length(positions)], rates[nrow(rates),2])

##plot recombination rates and expected numbers of arrivals above each other
##overlaying the plots makes them hard to read

postscript(outfile)
split.screen(c(2,1))##create a plotting page with 2 panels
screen(1)##plot recombination rates on upper panel
plot(rates[,1], rates[,2], xlab="position (Mb)", ylab="Recomb rate (cM/Mb)",
     main="Hotspot map of chromosome 22", type='l')

screen(2)##plot posterior means on lower panel
plot(positions[-1], lambda, col=2, type='l', xlab="", ylab="Exp Num Arrivals")
close.screen(all=T)
dev.off()

