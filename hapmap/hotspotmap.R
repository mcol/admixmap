## Script to plot posterior means of lambda parameters (expected numbers of arrivals between loci) in hapmix model along with recombination rates
## from HapMap
##currently only does one chromosome
##TODO: make a collection of plots for the whole genome; also plot subsets of chromosomes
## NOTE: the loci used in ADMIXMAP and the map may be different. This is accounted for in the script

################ filenames
Eur.lambdafile <- "Eur/Results6States/lambdaPosteriorMeans.txt"##file with lambda posterior means
Afr.lambdafile <- "Afr/Results6States/lambdaPosteriorMeans.txt"
Asian.lambdafile <- "Asian/Results6States/lambdaPosteriorMeans.txt"

##locus file as used with ADMIXMAP, for distances
Eur.locusfile <- "Eur/chr22data/loci5000_phased.txt"
Afr.locusfile <- "Afr/chr22data/loci5000_phased.txt"
Asian.locusfile <- "Asian/chr22data/loci5000_phased.txt"

##HapMap legend file used to generate locusfile, for positions
Eur.legendfile <- "Eur/chr22_legend.txt"                    # (only for phased data)
Afr.legendfile <- "Afr/chr22_legend.txt"
Asian.legendfile <- "Asian/chr22_legend.txt"

ratefile   <- "genetic_map_chr22.txt"        ##recombination rates from HapMap
outfile <- "HotspotMap.ps"

################## Script starts here

#read rates from ratefile
rates<-read.table(ratefile, header=T)[,1:2]
##rescale positions to start at 0 and measure in Mb (HapMap files measure in bp)
##rates[,1]<-rates[,1]-pos0
##rates[,1]<-rates[,1]/1e06
pos0<-as.numeric(rates[1,1])

##rates.cM <- rates[,2]/1e6
rates.data <- data.frame(rate=as.numeric(rates[,2]*1e6), index=1:nrow(rates))##rates converted to cM
dimnames(rates.data)[[1]] <- as.character(rates[,1])##use positions as labels
rm(rates)

##read posterior means
##Eur
##Read positions in the raw data file *_legend.txt
Eur.distances<-read.table(Eur.locusfile, header=T, comment.char="", na.strings="#")[-1,3]
Eur.positions <- read.table(Eur.legendfile, header=T, )[1:(length(Eur.distances)+1),2]
Eur.positions <- as.character(Eur.positions)

Eur.lambda <-as.vector(read.table(Eur.lambdafile, header=F)[,1])
Eur.lambda <- data.frame(Eur.lambda / Eur.distances)
dimnames(Eur.lambda)[[1]] <- Eur.positions[-1]
rm(Eur.distances)

##Afr
Afr.distances<-read.table(Afr.locusfile, header=T, comment.char="", na.strings="#")[-1,3]
Afr.positions <- read.table(Afr.legendfile, header=T, )[1:(length(Afr.distances)+1),2]
Afr.positions <- as.character(Afr.positions)

Afr.lambda <-as.vector(read.table(Afr.lambdafile, header=F)[,1])
Afr.lambda <- data.frame(Afr.lambda / Afr.distances)
dimnames(Afr.lambda)[[1]] <- Afr.positions[-1]
rm(Afr.distances)

##Asians
Asian.distances<-read.table(Asian.locusfile, header=T, comment.char="", na.strings="#")[-1,3]
Asian.positions <- read.table(Asian.legendfile, header=T, )[1:(length(Asian.distances)+1),2]
Asian.positions <- as.character(Asian.positions)
all.positions <- c(all.positions, Asian.positions)
Asian.lambda <-as.vector(read.table(Asian.lambdafile, header=F)[,1])
Asian.lambda <- data.frame(Asian.lambda / Asian.distances)
dimnames(Asian.lambda)[[1]] <- Asian.positions[-1]
rm(Asian.distances)


##extract elements at positions common to all 3 pops and average
all.positions <- c(Eur.positions, Afr.positions)
positions <- sort(all.positions[duplicated(all.positions)])
all.positions <- c(positions, Asian.positions)
positions <- sort(all.positions[duplicated(all.positions)])
rm(Eur.positions)
rm(Afr.positions)
rm(Asian.positions)
rm(all.positions)
lambda <- cbind(Eur.lambda[positions,1],
                Afr.lambda[positions,1],
                Asian.lambda[positions,1])
avg.lambda <- apply(lambda, 1, mean)
rm(Eur.lambda)
rm(Afr.lambda)
rm(Asian.lambda)

##read distances from locusfile
##distances<-read.table(locusfile, header=T, comment.char="", na.strings="#")[,3]
##convert distances in locusfile into positions in Mb
##NB: assumes the first locus is the same in both files
##positions <- c(0, cumsum(distances[-1]))+pos0

##adjust rates by removing loci not used in analyses

##rates.modified <- sum(rates.data[ 1:(rates.data$as.character(positions[1])) ] )
rates.modified <- numeric(length(positions))
rates.modified[1] <- rates.data[1,1]##first element will be ignored anyway (should be 0)
from.index <- 1
for(d in 2:length(positions)){
  to.index <- rates.data[positions[d],2]
  rates.modified[d] <- sum( rates.data[from.index : to.index, 1]   )
  from.index <- to.index+1
}
rm(rates.data)
##convert back to cM/Mb
rates.modified <- rates.modified / 1e6
##shift positions back to start at 0 and rescale to Mb
positions <- (as.numeric(positions) - pos0)/1e6

##plot recombination rates and expected numbers of arrivals above each other
##overlaying the plots makes them hard to read
postscript(outfile)
split.screen(c(2,1))##create a plotting page with 2 panels
screen(1)##plot recombination rates on upper panel
plot(positions, rates.modified, xlab="position (Mb)", ylab="Recomb rate (cM/Mb)",
     main="Hotspot map of chromosome 22", type='l')

screen(2)##plot posterior means on lower panel
plot(positions, avg.lambda, type='l', col=2, xlab="", ylab="Exp Num Arrivals")
close.screen(all=T)

##scatter plot of lambda vs rates
##plot(log(avg.lambda), log(rates.modified[-1]), xlab="lambda", ylab="Recomb rate")
dev.off()

##max.posn <- max(positions[length(positions)], rates[nrow(rates),2])




