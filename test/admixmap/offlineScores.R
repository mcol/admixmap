
N <- 200
L <- 50
genotypes <- as.matrix(read.table("simdata/genotypes.txt", header=T, nrow=N, colClasses="character")[,-c(1:2)])
g <- matrix(0, nrow=N, ncol=L)
g[genotypes=="1,1"] <- 2
g[genotypes=="1,2"] <- 1
g[genotypes=="2,1"] <- 1
g[genotypes=="2,2"] <- 0
##allele1 => ancestry state 2

mean.counts.allele1 <- apply(g, 2, mean)
mean.counts.allele2 <- apply(2-g, 2, mean)

admix.post.means <- read.table("simResultsUnlinked3.6/IndAdmixPosteriorMeans.txt", header=T, nrow=N)
mean.admix <- colMeans(admix.post.means)

scores <- N*(mean.counts.allele1 - 2*mean.admix[1])
