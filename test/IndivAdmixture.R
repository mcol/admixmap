library(MASS)
setwd("C:/CVS/genepi/admix/test")
args <- read.table("Indresults//args.txt", sep="", header=FALSE)
user.options <- as.data.frame(matrix(data=NA, nrow=1, ncol=dim(args)[1]))
user.options[1,] <- args[,2]
dimnames(user.options) <- list("Value", args[,1])
## admixmap does not write the option randommatingmodel to file args.txt - should fix
randommatingmodel <- as.vector(1, mode="numeric")
user.options <- data.frame(user.options, randommatingmodel)
thinning <- as.numeric(user.options$every)
results <- "Indresults/"


## read posterior samples of individual vars
  if(!is.null(user.options$indadmixturefile)) {
  samples <- dget(user.options$indadmixturefile)
  K<-(dim(samples)[[1]]-3)/2
  populations.labels<-dimnames(samples)[[1]][seq(from=1,to=2*K,by=2)]

  n.iterations <- dim(samples)[3]
  n.individuals <- dim(samples)[2]
  ## rearrange as 4-way array if random-mating model
  if(!is.null(user.options$randommatingmodel) && user.options$randommatingmodel==1) {
    samples <- array(as.vector(samples), dim=c(2, K + 1, n.individuals, n.iterations))
# K if globalrho
    ## re-assign dimnames
#    if(!is.null(user.options$globalrho) && user.options$globalrho==0) {
      dimnames(samples) <- list(c("Parent1", "Parent2"),
                                c(population.labels, "sumIntensities"), NULL, NULL)
#    } else {
#    dimnames(samples) <- list(c("Parent1", "Parent2"), population.labels, NULL, NULL)
#  }
    ## identify parent1 as parent with more ancestry from population 1 
    moreAfr.parent <- 2 - as.numeric((samples[1,1,,] - samples[2,1,,]) > 0)
                                        # evaluates as 1 if parent1 more Afr, 2 otherwise
    samples.1moreAfr <- samples
    for(iteration in 1:n.iterations) {
      if(moreAfr.parent[iteration]==2) {
        samples.1moreAfr[1, , , iteration] <- samples[2, , , iteration]
        samples.1moreAfr[2, , , iteration] <- samples[1, , , iteration]
      }
    }
    s.quantiles <- array(data=NA, dim=c(K + 1, 4, 2), # K if globalrho
                         dimnames=list(c(population.labels, "sumIntensities"),
                           c("Mean", "Median", "2.5%", "97.5%"), c("Parent1", "Parent2")))
    s.quantiles[, 1, ] <- t(apply(samples.1moreAfr, 1:2, mean))
    for(parent in 1:2) {
      for(row in 1:(K + 1)) {
        s.quantiles[row, 2:4, parent]  <-
          quantile(samples.1moreAfr[parent,row,,], probs = c(0.5, 0.025, 0.975))
      }
    }
    print(round(s.quantiles, digits=2))
    write.table(round(s.quantiles[,,1], digits=2), file="Indresults/quantiles.txt", quote=numeric(0))
    write.table(round(s.quantiles[,,2], digits=2), file="Indresults/quantiles.txt", quote=numeric(0),
                append=TRUE)
    samples.meanparents <- apply(samples, 2:4, mean) # array of dim K, ? may need aperm
  } else {
    samples.meanparents <- samples  # array of dim K, 1, n.iterations
  }
  s.quantiles.meanparents <- array(data=NA, dim=c(K+1, 4),
                                   dimnames=list(c(population.labels, "sumIntensities"),
                                     c("Mean", "Median", "Pct2.5", "Pct97.5")))
  s.quantiles.meanparents[, 1] <- apply(samples.meanparents, 1:2, mean)
  for(pop in 1:K) {
    s.quantiles.meanparents[pop, 2:4]  <-
      quantile(samples.meanparents[pop,,], probs = c(0.5, 0.025, 0.975))
  }
  
  write.table(round(s.quantiles.meanparents, digits=2),
                file="Indresults/quantilesMeanParents.txt", quote=numeric(0))
  
  ## plot posterior densities of individual admixture proportions
  postscript("Indresults/IndivAdmixture.ps")
  popcols <- c("grey", "blue", "red", "yellow", "orange", "green")
  for(pop in 1:K) {
    if(!is.null(user.options$randommatingmodel) && user.options$randommatingmodel==1) {
      plot(density(0.5*(samples[1, pop, 1, ] + samples[2, pop, 1, ]),
                   adjust=0.5, from=0, to=1),
           main=paste(population.labels[pop], "admixture proportion: mean of both parents"),
           ylab="Posterior density", xlim=c(0,1))
      plot(density(abs(samples[1, pop, 1, ] - samples[2, pop, 1, ]),
                   adjust=0.5, from=0, to=1), 
           main=paste(population.labels[pop], "admixture proportion: difference between parents"),
           ylab="Posterior density", xlim=c(0,1))
      plot(density(samples.1moreAfr[1, pop, 1, ],
                   adjust=0.5, from=0, to=1),
           main=paste(population.labels[pop], "admixture proportion of parent with more Afr ancestry"),
           ylab="Posterior density", xlim=c(0,1))
      plot(density(samples.1moreAfr[2, pop, 1, ],
                   adjust=0.5, from=0, to=1), 
           main=paste(population.labels[pop], "admixture proportion of parent with less Afr ancestry"),
           ylab="Posterior density", xlim=c(0,1))
      parents.pop <- kde2d(samples[1, pop,,],samples[2,pop, ,], lims=c(0,1,0,1))
      contour(parents.pop$x, parents.pop$y, parents.pop$z,
            main=paste("Contour plot of posterior density of parental", population.labels[pop],
              "admixture proportions"),
              xlab="Parent 1", ylab="Parent 2")
      persp(parents.pop$x, parents.pop$y, parents.pop$z, col=popcols[pop],
            main=paste("Perspective plot of bivariate density of parental", population.labels[pop],
              "admixture proportions"),
            xlab="Parent 1", ylab="Parent 2", zlab="Posterior density")
    } else {
      plot(density(samples[pop, 1, ],
                   adjust=0.5, from=0, to=1),
           main=paste(population.labels[pop], "admixture proportion: assumed same for both parents"),
           ylab="Posterior density", xlim=c(0,1))
    }
  }
  dev.off()
  
  ## plot samples
  postscript("Indresults/PosteriorSamples.ps")
  for(pop in 1:K) {
    if(!is.null(user.options$randommatingmodel) && user.options$randommatingmodel==1) {
      plot(seq(1:n.iterations), samples[1, pop, 1, ], type="l", col="red",
           ylim=c(0,1),
           main=paste("Sample plot", population.labels[pop]),
           xlab="Iterations", ylab="Admixture proportion")
      lines(seq(1:n.iterations)/thinning, samples[2, pop, 1, ], col="blue")
    }
  }
  dev.off()
 
} # ends block conditional on indadmixturefile
