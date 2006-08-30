# rm(list = ls())  ## remove (almost) everything in the working environment.
library(MASS)
## script should be invoked from folder one level above subfolder specified by resultsdir
## to run this script from an R console session, set environment variable RESULTSDIR
## by typing 'Sys.putenv("RESULTSDIR" = "<path to directory containing results>")'
message <- "\n\nStarting R script\n";

if(nchar(Sys.getenv("RESULTSDIR")) > 0) {
  resultsdir <- Sys.getenv("RESULTSDIR")
} else {
  ## resultsdir set to default directory 
  resultsdir <- "results"
}

cbindIfNotNull <- function(table1, table2) {
## cbind two tables if both are not null
  if(!is.null(table1) && !is.null(table2)) {
    table3 <- cbind(table1, table2)
  } else if(!is.null(table1)) {
    table3 <- table1
  } else if(!is.null(table2)) {
    table3 <- table2
  } else {
    table3 <- NULL
  }
  return(table3)
}

getUserOptions <- function(argsfilename) {
  ## read table of user options
  if(!file.exists(argsfilename)){
    message <- c(message,"Error: cannot find argsfile\n")
    ##cat("Error: cannot find argsfile\n", file=outfile, append=T)
    quit(save="no", status=1, runLast=F)
  }else{
    args <- read.table(argsfilename, sep="=", header=FALSE, comment.char="")
    user.options <- as.data.frame(matrix(data=NA, nrow=1, ncol=dim(args)[1]))
    user.options[1,] <- args[,2]
    dimnames(user.options) <- list("Value", args[,1])
    user.options$every <- as.numeric(user.options$every)
    return(user.options)
  }
}

getAdmixturePrior <- function(K, user.options) {
  AdmixturePrior <- matrix(data=1, nrow=2, ncol=K)
  if(!is.null(user.options$admixtureprior)) {
    AdmixturePrior[1, ] <- as.numeric(unlist(strsplit(user.options$admixtureprior, ",")))
    AdmixturePrior[2, ] <- as.numeric(unlist(strsplit(user.options$admixtureprior, ",")))
  }
  if(!is.null(user.options$admixtureprior1)) {
    AdmixturePrior[2, ] <- as.numeric(unlist(strsplit(user.options$admixtureprior1, ",")))
  }
  return(AdmixturePrior)
}

getSumIntensitiesPrior <- function(user.options) {
  SumIntensitiesPrior <- numeric(2)
  if(!is.null(user.options$sumintensitiesprior)) {
    SumIntensitiesPrior <- as.numeric(unlist(strsplit(user.options$globalsumintensitiesprior, ",")))
   }
  return(SumIntensitiesPrior)
}

getParentsIdentified <- function(AdmixturePrior) {
  return(sum(AdmixturePrior[1, ] != AdmixturePrior[2, ]))
}

getIsAdmixed <- function(AdmixturePrior) {
  IsAdmixed <- logical(2)
  IsAdmixed[1] <- sum(AdmixturePrior[1, ]==0) < dim(AdmixturePrior)[2] - 1
  IsAdmixed[2] <- sum(AdmixturePrior[2, ]==0) < dim(AdmixturePrior)[2] - 1
  return(IsAdmixed)
}
  
readLoci <- function(){
  filename <- paste(resultsdir, "LocusTable.txt", sep="/")
  ##cols are locus name, number of allele/haplotypes, map distance in cM, chromosome number
  loci.compound <- read.table(file=filename, header=T, colClasses=c("character", "integer", "numeric", "integer"))
  return(loci.compound)
}

getNumSubpopulations <- function(user.options) {
  if(!is.null(user.options$populations)) {
    K <- as.numeric(user.options$populations)
  }
  else {cat("Error: number of subpopulations not written to file\n", file=outfile, append=T);}
  return(K)
}

getPopulationLabels <- function(k, user.options) {
  if(k==1) {
    population.labels <- "SinglePop"
  } else {
    if(!is.null(user.options$paramfile) &&
       file.exists(paste(resultsdir,user.options$paramfile,sep="/")) &&  
       length(scan(paste(resultsdir,user.options$paramfile, sep="/"),
                   what='character', quiet=TRUE)) != 0) {
      population.labels <-
        dimnames(read.table(paste(resultsdir, user.options$paramfile, sep="/"),
                            header=TRUE))[[2]][1:k]
    } else {
      if(!is.null(user.options$priorallelefreqfile)) {
        population.labels <- dimnames(read.table(user.options$priorallelefreqfile,
                                                 header=TRUE))[[2]][-1]
      } else {
        if(!is.null(user.options$historicallelefreqfile)) {
          population.labels <- dimnames(read.table(user.options$historicallelefreqfile,
                                                   header=TRUE))[[2]][-1]
        }
        else {
          population.labels <- rep("Pop",K)
          for(i in 1:K)population.labels[i] <- paste(population.labels[i], format(i),sep="")
        }        
      }
    }
  }
  return(population.labels)
}    

getOutcomeType <- function(header) {
  ## determine whether continuous outcome by testing for column headed "precision"
  if(match("precision", header, nomatch=0) == 0) 
    outcome.continuous <- 0 else outcome.continuous <- 1
  return(outcome.continuous)
}

getNumCovariates <- function(user.options) {
  ## n.covariates is set to 0 if covariatesfile option not specified 
  n.covariates <- 0
  if(!is.null(user.options$covariatesfile)) {
    covariates <-  read.table(user.options$covariatesfile, sep="", na.strings=c("NA", "#"),
                              comment.char="*", header=TRUE)
    n.covariates <- dim(covariates)[2]
  }
  return(n.covariates)
}

getRegressionParamsForAdmixture <- function(user.options, k, n.covariates, population.labels) {
  ## number col.coeff1 of column containing regression coefficient for effect of subpopulation 1
  ## is determined as  n.covariates + 1
  col.coeff1 <- n.covariates + 1   #  + 1 for intercept
  col.coefflast <- col.coeff1 + k - 2  # only k- 1 coefficients
  beta <- data.frame(read.table(paste(resultsdir,user.options$regparamfile,sep="/"),
                                sep="",header=TRUE))[, col.coeff1:col.coefflast]
  dimnames(beta)[[2]] <- paste("slope", population.labels[-1], sep="." )
  return(beta)
}

getPrecision <- function(user.options) {
  beta <- read.table(paste(resultsdir,user.options$regparamfile,sep="/"),sep="",header=TRUE)
  beta <- data.frame(beta[,ncol(beta)])
  dimnames(beta)[[2]] <- "Precision"
  return(beta)
}

plotAutocorrelations <- function(table.samples, thinning) {
  ## plot autocorrelations
  ## drop columns with zero variance or NaN
  vcols <- apply(table.samples,2,var)
  pcols <- !is.nan(vcols) & vcols > 0
  table.samples <- table.samples[, pcols]
  if(is.null(dim(table.samples))) {
    table.samples <- data.frame(table.samples)
  }
  numcols <- dim(as.matrix(table.samples))[2]
  if(numcols > 0) {
    for(j in 1:dim(table.samples)[2]) {
      ## next line generates warning messages
      acfRes <- acf(table.samples[, j], plot=FALSE)
      if(length(acfRes$acf[is.nan(acfRes$acf)]) == 0){
        plot(thinning*acfRes$lag, abs(acfRes$acf), type="l", ylim = c(0, 1),
             main=paste("Autocorrelation plot ", dimnames(table.samples)[[2]][j]),
             xlab=paste("Lag"), ylab="Autocorrelation") 
      }
    }
  }
}

getiter <- function (link, iter) 
{
    result <- NULL
    idx <- is.element(dimnames(link)[[1]], iter)
    if (any(idx)) 
        result <- link[idx, , drop = FALSE]
    return(result)
}

gewekePwr <-function (link) 
{
    spec <- NULL
    pnames <- dimnames(link)[[2]]
    n <- nrow(link)
    nspans <- min(1 + sqrt(n)/0.3, n - 1)
    if (n > 2) {
        for (i in pnames) {
            spec <- c(spec, spec.pgram(link[, i], spans = nspans, 
                demean = TRUE, detrend = FALSE, plot = FALSE)$spec[1])
        }
        spec <- 10^(spec/10)
    }
    else {
        spec <- rep(NA, length(pnames))
    }
    return(structure(spec, names = pnames))
}

getgeweke<-function (link, p.first, p.last) 
{
    iter <- unique(as.numeric(dimnames(link)[[1]]))
    n <- length(iter)
    link.first <- getiter(link, iter[1:round(p.first * n)])
    link.last <- getiter(link, iter[(n - round(p.last * n) + 
        1):n])
    result <- (colMeans(link.first) - colMeans(link.last))/sqrt(gewekePwr(link.first)/nrow(link.first) + 
        gewekePwr(link.last)/nrow(link.last))
    result <- cbind(result, 2 * (1 - pnorm(abs(result))))
    dimnames(result)[[2]] <- c("Z-Score", "p-value")
    return(result)
}

checkConvergence <- function(table.samples, listname, outputfile) {
  ## test for adequate burn-in based on Geweke (1992)
  ## in: Bayesian Stats 4. ed. Bernardo JM et al, Oxford, OUP)
  table.samples <- table.samples[, apply(table.samples,2,var)>0]
  if(is.null(dim(table.samples))) {
    table.samples <- data.frame(table.samples)
  }
  numcols <- dim(as.matrix(table.samples))[2]
  if(numcols > 0) {
    table.geweke <- getgeweke(link=table.samples, p.first=0.1, p.last=0.5)
    table.geweke <- round(table.geweke, digits=3)
    write.table(table.geweke, file=outputfile, sep="\t")
  }
}

plotErgodicAverages <- function(ergodicaveragefile, thinning) {
  table.averages <- read.table(file=ergodicaveragefile, header=TRUE)
  if(dim(table.averages)[1] > 5) {
    iterations <- 10*thinning*seq(1:dim(table.averages)[1])	
    omit <- seq( 1:(floor(length(iterations)/5)) ) #exclude first 20% from plot
    ## plot ergodic averages
    start <- 1
    if( K > 1 && user.options$indadmixhiermodel==1) {
      start <- K+1
      dispersion <- apply(table.averages[-omit,1:K], 1, sum)
      plot(iterations[-omit],dispersion, type='l', main="Running posterior mean", 
           xlab="Iterations", ylab="PopAdmix Dispersion");
      for(j in 1:K){
        plot(iterations[-omit], table.averages[-omit,j]/dispersion, type="l", 
             main="Running posterior mean", 
             xlab="Iterations", ylab=dimnames(table.averages)[[2]][j]) # should fix dimnames
      }
    }
    for(j in start:dim(table.averages)[2]) {
      plot(iterations[-omit], table.averages[-omit,j], type="l", 
           main="Running posterior mean", 
           xlab="Iterations", ylab=dimnames(table.averages)[[2]][j]) # should fix dimnames
    }
  }
}

popAdmixProportions <- function(population.labels, params, k) {
  ## calculate population admixture proportions from Dirichlet parameters
  pop.admix.prop <- as.data.frame(matrix(data=NA, nrow=n, ncol=k))
  for(j in 1:k) {
    dimnames(pop.admix.prop)[[2]][j] <- paste("prop", population.labels[j], sep=".")
    for(i in 1:n) {
      pop.admix.prop[i,j] <- params[i,j]/sum(params[i,1:k])
    }
  }
  return(pop.admix.prop)
}

effectEstimates <- function(beta, pop.admix.prop, n, k) {
  ## effect estimate is difference of 2 expectations
  ## E(Y | M[i]=1) = intercept + beta[i]
  ## E(Y | M[i]=0) = intercept + E(Y | M[j]=1) averaged over all j != i
  ## with weights given by Dirichlet parameters or population admixture proportions
  ## E(Y | M[1]=1) = intercept (because 1st population is baseline category and beta[1] = 0) 
  ## intercept cancels in subtraction
  ## with k populations, matrix of effects given M=1 is defined by betas augmented with 1st column of 0s
  beta <- data.frame(rep(0, n), beta)
  w.beta      <- matrix(data=NA, nrow=n, ncol=k,
                        dimnames=list(character(0), population.labels))
  sum.w.beta  <- matrix(data=NA, nrow=n, ncol=k,
                        dimnames=list(character(0), population.labels))
  sum.w       <- matrix(data=NA, nrow=n, ncol=k,
                        dimnames=list(character(0), population.labels))
  effect.Meq0 <- matrix(data=NA, nrow=n, ncol=k,
                        dimnames=list(character(0), population.labels))
  effect.pop  <- matrix(data=NA, nrow=n, ncol=k,
                        dimnames=list(character(0), population.labels))
  for(pop in 1:k) {
    dimnames(effect.pop)[[2]][pop] <- paste("EffectVsAllOther.", population.labels[pop], sep="")
  }
  
  for(pop in 1:k) {
    ## calculate for each pop a vector w.beta as weight[pop]*beta[pop]
    w.beta[,pop] <- pop.admix.prop[,pop]*beta[,pop]
  }
  for(pop in 1:k) {
    ## calculate for each pop a vector of numerators by summing w.beta over all columns except pop
    sum.w.beta[,pop] <- apply(w.beta[,-pop], 1, sum)
    ## calculate for each pop a vector of denominators by summing weights over all columns except pop
    sum.w[,pop] <- apply(pop.admix.prop[,-pop], 1, sum)
  }
  for(pop in 1:k) {
    ## calculate for each pop a vector E(Y | M[pop]=0) by dividing numerators by denominators
    effect.Meq0[,pop] <- sum.w.beta[,pop]/sum.w[,pop]
    effect.pop[,pop] <- beta[,pop] - effect.Meq0[,pop]
  }
  return(effect.pop)
}

calculateAndPlotQuantiles <- function(param.samples, nvars) {
  ## calculate means and 95% credible intervals
  post.quantiles <- matrix(data=NA, nrow=nvars, ncol=4, 
                           dimnames=list(dimnames(param.samples)[[2]], 
                             c("Median", "Mean", "Pct2.5", "Pct97.5")))
  outputfile <- paste(resultsdir, "ParameterPosteriorDensities.ps", sep="/" )
  postscript(outputfile)
  for(j in 1:nvars) {
    post.quantiles[j, c(1,3,4)] <- quantile(param.samples[,j], probs = c(0.5, 0.025, 0.975), na.rm = T)
    post.quantiles[j, 2] <- mean(param.samples[,j], na.rm=T) 
    ## plots kernel density of each variable
    if(is.finite(post.quantiles[j, 2])) {
      plot(density(param.samples[,j], adjust=0.5), 
           main="Posterior kernel density", xlab=dimnames(post.quantiles)[[1]][j])
    }
  }
  dev.off()
  outputfile <- paste(resultsdir, "PosteriorQuantiles.txt", sep="/" )
  write.table(signif(post.quantiles, digits=3), file=outputfile, quote=FALSE, 
              row.names=TRUE, col.names=NA, sep="\t")
  return(post.quantiles)
}

plotlogpvalues <- function(psfilename, log10pvalues, table.every, title, hist ) {
  ## function to plot cumulative p-values
  ## arguments: psfilename, matrix of -log10pvalues, table.every, plot title
  postscript(psfilename)
  ntests <- dim(log10pvalues)[1]
  evaluations <- dim(log10pvalues)[2]
  ## generate plot for test 1
  plot(table.every*seq(1:evaluations), log10pvalues[1,],
       type="l", ylim=c(0, 4), ## c( 0, max(log10pvalues*is.finite(log10pvalues), na.rm=T) ),
       main=title, xlab="Iterations", ylab="-log10(p-value)")
  ## add lines to plot for all other tests 
  for(test in 2:ntests) {
    lines(table.every*seq(1:evaluations), log10pvalues[test,], type="l")
  }
  ## label lines for which final p<0.01
  log10pvalues.final <- log10pvalues[, dim(log10pvalues)[2]]
  siglabels <- dimnames(log10pvalues)[[1]][log10pvalues.final > 2 | is.na(log10pvalues.final)]
  siglabels.x <- 0.9*table.every*dim(log10pvalues)[2]
  siglabels.y <- as.vector(log10pvalues.final[log10pvalues.final > 2])
  siglabels.y[is.na(siglabels.y)] <- -1
  if(length(siglabels) > 0) {
    text(siglabels.x, siglabels.y, labels=siglabels, pos=3, offset=0.5) 
  }
  if(hist){
    ## histogram of p-values
    hist(10^(-log10pvalues[, evaluations]), xlim=c(0,1), xlab = "p-value", 
         main="Histogram of one-sided p-values")
  }
  dev.off()
}

plotPValuesKPopulations <- function(outfile, pvalues, thinning) {
  ## stdNormDev is a 3-way array: k populations, loci, draws 
  log10pvalues <- -log10(pvalues)
  outputfile <- paste(resultsdir, outfile, sep="/")
  outputfile <- paste(outputfile, ".ps", sep="")
  postscript(outputfile)
  colours <- c("black", "blue", "red", "green")
  header <- paste("Running computation of p-values for ",outfile,sep="")
  header <- paste(header," score tests",sep="");
  ## set up plot for first population and first locus
  plot(10*thinning*seq(1:dim(log10pvalues)[3]), log10pvalues[1,1,], type="l", ylim=c(0,5),
       main=header, 
       xlab="Iterations", ylab="-log10(p-value)")
  ## loop over populations and loci to add lines
  for(pop in 1:dim(log10pvalues)[1]) {
    for(locus in 1:dim(log10pvalues)[2]) {
      lines(10*thinning*seq(1:dim(log10pvalues)[3]), log10pvalues[pop, locus, ],
            type="l", col=colours[1 + ((pop-1) %% dim(log10pvalues)[1])])
      ## print(colours[1 + (test %% k)])
    }
  }
  ## label loci at which tests are significant
  siglabels.x <- 0.9*10*thinning*dim(log10pvalues)[3] ## x coordinate at which to write labels
  log10pvalues.final <-  log10pvalues[,,dim(log10pvalues)[3]]
  siglabels.y <- as.vector(log10pvalues.final[log10pvalues.final > 2]) # label only where p < 0.01
  siglabels.y[is.na(siglabels.y)] <- -1
  siglabels <- dimnames(log10pvalues.final)[[2]][log10pvalues.final > 2]
  if(length(siglabels.y) > 0) {
    text(siglabels.x, siglabels.y, labels=siglabels, pos=3, offset=0.5) 
  }
  dev.off()
}

##used for allelic and haplotype association score tests
plotScoreTest <- function(scorefile, haplotypes, outputfilePlot, thinning) {
  scoretest <- dget(paste(resultsdir,scorefile,sep="/"))
  
  ## rows are: locus,(haplotype), -log10pvalue
  ## extract testnames and drop 1st row
  if (!haplotypes) {
    testnames <- as.vector(scoretest[1,,1])
    ## convert pvalues to numeric
    log10pvalues <- matrix(as.numeric(scoretest[-1, ,]), nrow=dim(scoretest)[[2]], ncol=dim(scoretest)[[3]])
  } else {
    testnames <- as.vector(paste(scoretest[1,,1], scoretest[2,,1]))
    log10pvalues <- matrix(as.numeric(scoretest[-c(1:2), ,]), nrow=dim(scoretest)[[2]], ncol=dim(scoretest)[[3]])
  }

  dimnames(log10pvalues) <- list(testnames, NULL)
  log10pvalues[is.nan(log10pvalues)] <- NA

  plotlogpvalues(outputfilePlot, log10pvalues,
              10*thinning, "Running computation of p-values for allelic association", F)
}

#used to plot output of score test for heterozygosity
plotHWScoreTest <- function(scorefile, k) {
  scoretest <- read.table(paste(resultsdir,scorefile,sep="/"),header=TRUE, row.names="Locus")
  
  #qq plot of scores
  point.list <- scoretest[,6]
  if(length(point.list[!is.na(point.list)]) > 0){
    outputfile <- paste(resultsdir, "QQPlotHWTest.ps", sep="/" )
    postscript(outputfile)
    title <- "QQ plot of H-W Test z-scores"
    point.list <- qqnorm(scoretest[,6], main = title)
    lines(x = c(min(point.list$x,na.rm=T), max(point.list$x,na.rm=T)), y = c(min(point.list$x,na.rm=T), max(point.list$x,na.rm=T)))
    dev.off()
  }
}
  
## used to plot output of Rao-Blackwellized score tests for ancestry association and affectedsonly
plotAncestryScoreTest <- function(scorefile, testname, Pops, population.labels, thinning) {
  KK <- Pops
  poplabels <- population.labels
  if(Pops == 2 ) {
    KK <- 1
    poplabels <- population.labels[2]
  }
  scoretests <- dget(paste(resultsdir,scorefile,sep="/"))
  ## extract first row containing locus names
  locusnames <- scoretests[1, seq(1, dim(scoretests)[2], by=KK), 1]
  testnames <- paste(scoretests[1,,1], scoretests[2,,1])
  ## drop first two rows and reformat as 4-way array
  scoretests.n <- array(data=scoretests[-c(1:2),,],dim=c(dim(scoretests)[1]-2,dim(scoretests)[2],dim(scoretests)[3]))
  dim3way <- dim(scoretests.n)
  dim4way <- c(dim3way[1], KK, dim3way[2]/KK, dim3way[3])
  dimnames4way <- list(dimnames(scoretests)[[1]][-c(1:2)], poplabels, locusnames, NULL)
  scoretests4way <- array(as.numeric(as.vector(scoretests.n)), dim=dim4way, dimnames=dimnames4way)
  scoretests4way[is.nan(scoretests4way)] <- NA
  
  ## plot cumulative p-values in K colours
  pvalues <- array(data=scoretests4way[8,,,],dim=c(dim(scoretests4way)[2:4]),dimnames=c(dimnames(scoretests4way)[2:4]))
  if(length(!is.na(pvalues)) > 0) {
    plotPValuesKPopulations(testname, pvalues, thinning)
  }
  ## extract final table as 3-way array: statistic, locus, population
  
  scoretest.final <- array(data=scoretests4way[,,,dim(scoretests4way)[4]],dim=c(dim(scoretests4way)[1:3]),
                           dimnames=c(dimnames(scoretests4way)[1:3]))
  
  ## set test statistic to missing if obs info < 1
  if(getOutcomeType(dimnames(param.samples)[[2]]) == 1){#continuous outcome
    scoretest.final[7,,][scoretest.final[3,,] < 1] <- NA
  }
  ##plot z-scores across genome
  zscores <- array(data=scoretest.final[7,,], dim=c(dim(scoretest.final)[2:3]),dimnames=c(dimnames(scoretest.final)[2:3]))
  plotScoreMap(loci.compound,zscores, KK, testname) 
  ## plot information content
  info.content <- array(data=scoretest.final[4, , ],dim=c(dim(scoretest.final)[2:3]),dimnames=c(dimnames(scoretest.final)[2:3]))
  plotInfoMap(loci.compound, info.content, KK, testname)
  
  ## calculate high and low cutoffs of population risk ratio r that can be excluded at
  ## a likelihood ratio of 0.01
  ## based on quadratic approximation to the log-likelihood as a function of log r
  u <- as.vector(scoretest.final[1,1,])
  v <- as.vector(scoretest.final[3,1,])
  ## set to missing if obs info < 1
  v[v < 1] <- NA
  r.exclude.hi <- exp(u/v + sqrt(u^2 + 2*v*log(100))/v)
  r.exclude.lo <- exp(u/v - sqrt(u^2 + 2*v*log(100))/v)
  ## plotExclusionMap not implemented at present

  ##qq plot of scores
  outputfile <- paste(resultsdir, "QQPlot", sep="/" )
  outputfile <- paste(outputfile, testname, sep="")
  outputfile <- paste(outputfile, ".ps", sep="")
  
  postscript(outputfile)
  title <- paste("QQ plot of z-scores,", testname,sep="" )
  for(k in 1:(nrow(zscores))){
    if( length(zscores[k, ][!is.na(zscores[k, ])]) > 0 ) {
      point.list <- qqnorm(zscores[k,], main = title, sub=poplabels[k])
      lines(x = c(min(point.list$x,na.rm=T), max(point.list$x,na.rm=T)), y = c(min(point.list$x,na.rm=T),
                                                                           max(point.list$x,na.rm=T)))
    }
  }
  dev.off()
}

plotScoreMap <- function(loci.compound, zscores, K, testname){
  outputfile <- paste(resultsdir, testname, sep="/")
  outputfile <- paste(outputfile, "ScoreMap.ps", sep="")
  postscript(outputfile)
  for(chr in 1:n.chr) {
    for(pop in 1:K) {
      plot(loci.compound$MapPosition[loci.compound$Chromosome==chr],
           zscores[pop, loci.compound$Chromosome==chr], 
           type="l", ylim=c(-5,5),
           xlab="Map position (cM) from first locus", ylab="z-score",
           main=paste("z-scores for chromosome", chr, "- pop", pop)
           )
      lines(c(0,max(loci.compound$MapPosition[loci.compound$Chromosome==chr])), c(qnorm(0.975),qnorm(0.975)))
      lines(c(0,max(loci.compound$MapPosition[loci.compound$Chromosome==chr])), c(qnorm(0.995),qnorm(0.995)), lty=2)
      lines(c(0,max(loci.compound$MapPosition[loci.compound$Chromosome==chr])), c(qnorm(0.025),qnorm(0.025)))
      lines(c(0,max(loci.compound$MapPosition[loci.compound$Chromosome==chr])), c(qnorm(0.005),qnorm(0.005)), lty=2)
    }
  }
  dev.off()
}

plotExclusionMap <- function(loci.compound, info.content, cutoffs.lo, cutoffs.hi) {
  for(chr in 1:n.chr) {
    ## cutoffs plotted for risk ratio associated with 1st population only
    plot(loci.compound$MapPosition[loci.compound$Chromosome==chr],
         cutoffs.lo[loci.compound$Chromosome==chr], 
         type="l", ylim=c(0,10),
         xlab="Map position (cM) from first locus", ylab="Population risk ratio at LOD = -2",
         main=paste("Exclusion map for chromosome", chr))
    text(loci.compound$MapPosition[loci.compound$Chromosome==chr],
         cutoffs.lo[loci.compound$Chromosome==chr]
         + 0.5 + 2*(seq(1:length(loci.compound$Chromosome[loci.compound$Chromosome==chr])) %% 3),
         labels=as.vector(loci.compound[,1][seq(1, dim(scoretest.affectedsonly.final)[1]/k, by=k)]
           [loci.compound$Chromosome==chr]),adj=c(0,0), # pos=3,
         srt=90, offset=0.5)
  }
}

plotInfoMap <- function(loci.compound, info.content, K, testname) {
  ## info.content is matrix in which rows index populations, cols index loci
  outputfile <- paste(resultsdir, testname, sep="/")
  outputfile <- paste(outputfile, "InformationContentMap.ps", sep="")
  postscript(outputfile)
  for(chr in 1:n.chr) {
    for(pop in 1:K) {
      plot(loci.compound$MapPosition[loci.compound$Chromosome==chr],
           info.content[pop, loci.compound$Chromosome==chr], 
           type="l", ylim=c(-10,100),
           xlab="Map position (cM) from first locus", ylab="Percent info extracted",
           main=paste("Map information content for chromosome", chr, "- pop", pop)
           )
      ## text offset vertically in sequences of 3 for legibility
      text(loci.compound$MapPosition[loci.compound$Chromosome==chr],
           info.content[pop, loci.compound$Chromosome==chr] +  
           1 + 20*(seq(1:length(loci.compound$Chromosome[loci.compound$Chromosome==chr])) %% 3),
           labels=as.vector(loci.compound[,1][loci.compound$Chromosome==chr]),
           adj=c(0,0), # pos=3,
           srt=90, offset=0.5)
    }
  }
  dev.off()
}

plotResidualAllelicAssocScoreTest <- function(scorefile, outputfile, thinning){
  scoretest <- dget(paste(resultsdir,scorefile,sep="/"));
  ## dimensions are cols (2), pairs of loci, number of evaluations
  evaluations <- dim(scoretest)[3]
  ntests <- dim(scoretest)[2]
  locusnames <- scoretest[1,,evaluations]

  minuslog10pvalues <- as.numeric(scoretest[2, , ])
  minuslog10pvalues[is.nan(minuslog10pvalues)] <- NA
  minuslog10pvalues <- data.frame(matrix(data=minuslog10pvalues, nrow=ntests, ncol=evaluations))
  dimnames(minuslog10pvalues)[[1]] <- locusnames
  plotlogpvalues(outputfile, minuslog10pvalues, 10*thinning,
                 "Running computation of p-values for residual allelic association", T)
}

## offline score tests for genotypes at loci that have not been included in the model (because we can't model large haplotypes)
ExtraScoreTests <- function(testgenotypesfile, outcomevarfile, expectedoutcomefile){
  p <- dget(file=expectedoutcomefile)
  p <- p[, 1, ] # first outcome var
  testvars <- read.table(testgenotypesfile, header=T, as.is=T, sep="\t")
  y <- read.table(outcomevarfile, header=T, as.is=T, sep="\t")[, 1]
  
  scoretable <- matrix(nrow=0, ncol=4)
  numvars <- dim(testvars)[2]
  numdraws <- dim(p)[2]
  scoredraws <- numeric(numdraws)
  infodraws <- numeric(numdraws)
  for(j in 1:numvars) {
    x <- testvars[, j] - mean(testvars[, j], na.rm=T)
    for(draw in 1:numdraws) { # evaluate score and info for each realization of complete data 
      scoredraws[draw] <- sum((y - p[, draw])*x, na.rm=T)
      infodraws[draw] <- sum(p[, draw]*(1 - p[, draw])*x^2, na.rm=T)
    }
    ## evaluate expectations over posterior distribution 
    score <- mean(scoredraws)
    completeinfo <- mean(infodraws)
    missinfo <- var(scoredraws)
    obsinfo <- completeinfo - missinfo
    ## calculate test statistic
    z <- score / sqrt(obsinfo)
    result <- c(score, completeinfo, obsinfo, z)
    scoretable <- rbind(scoretable, result)
  }
  
  dimnames(scoretable)[[1]] <- dimnames(testvars)[[2]]
  dimnames(scoretable)[[2]] <- c("score", "complete.info", "observed.info", "z")
  write.table(scoretable, file=paste(resultsdir, "ExtraScoreTests.txt", sep="/"), sep="\t")
}

convertAlleleFreqs <- function(allelefreq.samples) {
  ## drop first row containing locus names
  converted <- allelefreq.samples[-1,,]
  ## if this is a 2-way array (single population), rearrange as 3-way array
  if(length(dim(converted))==2) {
    converted <- array(as.numeric(converted), dim=c(1, dim(converted)))
  }
  ## permute so that dimension 1 is alleles x loci, dimension 2 is pops, dimension 3 is draws
  converted <- aperm(converted, c(2,1,3))
  ## convert array to numeric
  converted <- array(as.numeric(converted), dim=dim(converted),
                              dimnames=dimnames(converted))
  return(converted)
} 

listAlleleFreqs <- function(allelefreq.samples) {
  ## argument is a 3-way array of allele frequencies (locusname + pops, alleles x loci, draws)
  ## converts to a list of 3-way arrays each holding alleles x pops x draws for one locus
  ## dimension 1 of each array is number of alleles minus 1 
  ## extract vector of locus names
  converted.samples <- convertAlleleFreqs(allelefreq.samples)
  draws <- dim(allelefreq.samples)[3]
  locusnames <- allelefreq.samples[1, , 1]
  ## generate vector with locusnames coded as levels of a factor
  row.locusnumber <- match(locusnames, unique(locusnames))
  
  ## generate list holding 3-way arrays (alleles-1 x pops x draws) of allele freqs at each locus
  allelefreq.samples.list <- list()
  ## loop over loci
  for(j in 1:nrow(loci.compound)) {
    ## firstrow is position in dimension 1 of first allele at j th locus 
    firstrow <- match(j, row.locusnumber)
    ## lastrow is position in dimension 1 of last allele (but one) at j th locus
    lastrow <- length(row.locusnumber) + 1 - match(j, rev(row.locusnumber))
    ## freqarray should always be 3-way array even if diallelic locus or single population 
    freqarray <- converted.samples[firstrow:lastrow, , ,drop=FALSE]
    ##if(is.vector(freqarray)) { 
      ##freqarray <- array(freqarray, dim=c(1, 1, draws))
    ##} else {
      ##if(length(dim(freqarray))==2) {
        ##freqarray <- array(freqarray, dim=c(1, dim(freqarray)[1], draws))
      ##}
    ##}
    allelefreq.samples.list[[j]] <- freqarray
  }
  return(allelefreq.samples.list)
}

calculateLocusfValues <- function(allelefreq.samples.list) {
  ## calculates locus information content for ancestry
  ## correct only for two populations
  ## otherwise uses first two populations
  num.loci <- length(allelefreq.samples.list)
  num.draws <- dim(allelefreq.samples.list[[1]])[3]
  locusfValues.samples <- matrix(data=NA, nrow=num.loci, ncol=1) # should have K(K+1)/2 cols
  f.draws <- numeric(num.draws)
  f.means <- numeric(num.loci)
  ## should loop over all pairs of pops
  pop1 <- 1
  pop2 <- 2
  ## 3-way array (alleles-1 x pops x draws) of allele freqs at each locus
  for(locus in 1:num.loci) {
    aminus1 <- dim(allelefreq.samples.list[[locus]])[1]
    for(draw in 1:num.draws) {
      ## a alleles, 2 pops
      p <- matrix(allelefreq.samples.list[[locus]][, c(pop1, pop2), draw], nrow=aminus1, ncol=2)
      p <- rbind(p, 1 - apply(p, 2, sum))
      ratios <- p[, 1]*p[, 2]/(p[, 1] + p[, 2]) ## vector of length a
      f.draws[draw] <- 1 - 2*sum(ratios)
    }
    ## locusfValues.samples.list[[locus]] <- f.draws
    f.means[locus] <- mean(f.draws)
  }
  return(f.means) 
  ## return(locusfValues.samples.list)  
}

entropy <- function(p) {
  return(-sum(p*log(p)))
}
  
calculateLocusKLInfo <- function(allelefreq.samples.list) {
  ## calculates locus KL information content for ancestry
  num.loci <- length(allelefreq.samples.list)
  num.draws <- dim(allelefreq.samples.list[[1]])[3]
  KLInfo.samples <- matrix(data=NA, nrow=num.draws, ncol=num.loci) 
  ## 3-way array (alleles-1 x pops x draws) of allele freqs at each locus
  for(locus in 1:num.loci) {
    for(draw in 1:num.draws) {
      ## a alleles, K pops
      p <- matrix(allelefreq.samples.list[[locus]],
                  nrow=dim(allelefreq.samples.list[[locus]])[1], ncol=K)
      p <- rbind(p, 1 - apply(p, 2, sum))
      pbar <- apply(p, 2, mean) # mean freq of allele over K populations
      KLInfo.samples[draw, locus]  <- entropy(pbar) - mean(apply(p, 2, entropy))
    }
    KLInfo.means <- apply(KLInfo.samples, 2, mean)
  }
  return(KLInfo.means) 
}

listFreqMeansCovs <- function(allelefreq.samples.list) {
  ## generate lists to hold allele freq means and covariances
  allelefreq.means.list <- list()
  allelefreq.covs.list  <- list()
  K <- dim(allelefreq.samples.list[[1]])[2]
  ## loop over loci to compute means and covariances
  for(locus in 1:length(allelefreq.samples.list)) {
    ## loop over populations 
    for(pop in 1:K) {
      ## add row for freq last allele to matrix of draws of allele freqs
      f.exceptlast <- allelefreq.samples.list[[locus]][,pop,]
      if(is.vector(f.exceptlast)) {
        f.exceptlast <- matrix(data=f.exceptlast, nrow=1)
      }
      f.last <- 1 - apply(f.exceptlast, 2, sum)
      f.all <- rbind(f.exceptlast, f.last) # matrix in which rows index alleles, cols draws
      ## calculate means and covariances over draws
      allelefreqs.mean <- apply(f.all,1 , mean) # vector of length a
      allelefreqs.cov <- cov(data.frame(t(f.all)))/( (ncol(f.all)+1)/ncol(f.all) )
                                        # covariance matrix of order a
                                        #to get biased covariances
      ## assign means and covariances to list elements
      if(pop==1) { 
        allelefreq.means.list[[locus]] <- list()
        allelefreq.covs.list[[locus]]  <- list()
      }
      allelefreq.means.list[[locus]][[pop]] <- allelefreqs.mean
      allelefreq.covs.list[[locus]][[pop]]  <- allelefreqs.cov
    }
  }
  return(list(allelefreq.means.list, allelefreq.covs.list))
}

fitDirichletParams <- function(allelefreq.means.list, allelefreq.covs.list) {
  ## fit Dirichlet parameters by equating posterior means and variances
  k <- length(allelefreq.means.list[[1]])
  allelefreq.params.list  <- list()
  allelefreq.sumalphas <- matrix(data=NA, nrow=length(allelefreq.means.list), ncol=k)
  for(locus in 1:length(allelefreq.means.list)) {
    allelefreq.params.list[[locus]] <- list()
    for(pop in 1:k){
        p <- allelefreq.means.list[[locus]][[pop]]#[-1]
        v <- allelefreq.covs.list[[locus]][[pop]] #[-1,-1]

      if(sum(v)>0){
        if(length(p)==2) {#diallelic locus
          factor <- p[1]*(1 - p[1])/v[1,1]
        } else {
          covar.predicted <- matrix(data=NA, nrow=length(p), ncol=length(p))
          for(i in 1:length(p)) {
            for(j in 1:length(p)) {
              ## predicted covariance if sum.alphas = 0
              covar.predicted[i,j] <- ifelse(i==j, p[i]*(1-p[i]), p[i]*p[j])
            }
          }
          d.predicted <- sum(diag(covar.predicted))##det(covar.predicted)
          d.observed <- sum(diag(v))##det(v)
          factor <- (d.predicted/d.observed)##^(1/length(p))
        }
##        if(pop==1) { # create matrix of allele freqs: rows index alleles, cols index populations  
##          allelefreq.params.list[[locus]] <- matrix(data=NA, nrow=length(p), ncol=k,
##                                                    dimnames=dimnames(allelefreq.means.list[[locus]][-1]))
##                                        #list(character(0), population.labels))
##        }
        allelefreq.sumalphas[locus,pop] <- factor - 1
        allelefreq.params.list[[locus]][[pop]] <-
          allelefreq.sumalphas[locus,pop]*allelefreq.means.list[[locus]][[pop]] 
      }else{
        allelefreq.params.list[[locus]][[pop]] <-allelefreq.means.list[[locus]][[pop]]
      }
    }##end pop loop
  }## end locus loop
  return(allelefreq.params.list)
}


writeAlleleFreqs <- function(allelefreq.params.list, k, loci.compound, population.labels,
                                  outputfile) {
  ## write parameters of allele freq distributions to file in priorallelefile format
  allelefreq.params <- numeric(0)
  allelefreq.params.names <- character(0)
  for(locus in 1:length(allelefreq.params.list)) {
    allelefreq.params.names <- c(allelefreq.params.names,
                                 rep(as.vector(loci.compound[locus, 1]), loci.compound[locus, 2]))

    freqs.atlocus <- numeric(0)
    for(j in 1: length(allelefreq.params.list[[locus]]))
      freqs.atlocus <- cbind(freqs.atlocus, allelefreq.params.list[[locus]][[j]])
    allelefreq.params <- rbind(allelefreq.params, freqs.atlocus)
  }
  allelefreq.params <- data.frame(allelefreq.params.names,
                                  round(allelefreq.params,digits=2),
                                  check.rows=FALSE)
  dimnames(allelefreq.params)[[2]] <- c("Locus", population.labels)
  write.table(allelefreq.params, file=outputfile, row.names=FALSE, sep="\t")
}

plotAdmixtureDistribution <- function(alphas, samples, k) {
  ## plots histogram of posterior means of individual admixture proportions for comparison 
  ## with population distribution given by Dirichlet parameter vector alphas
  xvalues <- seq(0.005,0.995,0.01)
  outputfile <- paste(resultsdir, "DistributionIndividualAdmixture.ps", sep="/" )
  postscript(outputfile)
  par(cex=1.5, las=1) 
  for(pop in 1:k) {
    popM <- dbeta(xvalues, alphas[pop], sum(alphas[-pop])) # beta density of pop at xvalues
    indivadmixture.estimates <- apply(samples[pop,,], 1, mean) # proportionate admixture from pop 
    hist(indivadmixture.estimates, xlim=c(0,1), # ylim=c(0, 2*max(popM)), 
         ## breaks=seq(0,1, 0.05), 
         freq=FALSE,  
         xaxs="i", yaxs="i", 
         xlab=paste("Proportion ", population.labels[pop], " admixture"), ylab="Frequency", 
         main="Distribution of individual admixture estimates")
    points(xvalues, popM, type="l", col="blue")
  }
  dev.off()
}

writePosteriorMeansIndivAdmixture <- function(samples, K) {
  ## compute posterior means of individual admixture and ancestry diversity
  ## dim 1 is populations, dim2 is individuals, dim3 is samples
  n.individuals <- dim(samples)[2]
  n.samples <- dim(samples)[3]
  
  M.squared <- samples.meanparents[1:K,1:n.individuals,1:n.samples, drop=F]^2
  ancestry.diversity <- 1 - apply(M.squared, 2:3, sum)
  AncestryDiversity <- apply(ancestry.diversity, 1, mean)
  IndividualAdmixture <- t(apply(samples.meanparents, 1:2, mean))
  admixture.table <- data.frame(IndividualAdmixture, AncestryDiversity)
  outputfile <- paste(resultsdir, "IndivAdmixPosteriorMeans.txt", sep="/" )
  write.table(signif(admixture.table, digits=3), file=outputfile, sep="\t", 
              quote=FALSE, row.names=FALSE)
}

plotPosteriorDensityIndivParameters <- function(samples.admixture, samples.sumIntensities, user.options,
                                                population.labels, AdmixturePrior, IsAdmixed, ParentsIdentified) {
  ## samples.admixture: 3-way array with dimensions that index subpopulations, parents, draws
  outputfile <- paste(resultsdir, "IndivParameters.ps", sep="/")
  postscript(outputfile)
  popcols <- c("grey", "blue", "red", "yellow", "orange", "green")
  if(IsAdmixed[1] & IsAdmixed[2]) { # both parents admixed
    if(ParentsIdentified) { # bivariate plots if both parents have ancestry from pop
      for(pop in 1:K) {
        if(AdmixturePrior[1, pop] > 0 & AdmixturePrior[2, pop] > 0) { # bivariate plot
          parents.pop <- kde2d(samples.admixture[pop, 1, ],samples.admixture[pop, 2, ], lims=c(0,1,0,1))
          contour(parents.pop$x, parents.pop$y, parents.pop$z,
                  main=paste("Contour plot of posterior density of parental", population.labels[pop],
                    "admixture proportions"),
                  xlab="Parent 1", ylab="Parent 2")
          persp(parents.pop$x, parents.pop$y, parents.pop$z, col=popcols[pop],
                main=paste("Perspective plot of bivariate density of parental", population.labels[pop],
                  "admixture proportions"),
                xlab="Parent 1", ylab="Parent 2", zlab="Posterior density")
        }
      }
      parents.pop <- kde2d(log(samples.sumIntensities[1, ]),
                           log(samples.sumIntensities[2, ]), lims=c(-1,3,-1,3))
      contour(parents.pop$x, parents.pop$y, parents.pop$z,
              main="Contour plot of posterior density of parental sum-intensities",
              xlab="Parent 1", ylab="Parent 2")
      persp(parents.pop$x, parents.pop$y, parents.pop$z, col="blue",
            main="Perspective plot of bivariate density of parental sum-intensities",
            xlab="Parent 1", ylab="Parent 2", zlab="Posterior density")
    } else { # parents unidentified - concatenate array and do bivariate plots 
      samples.admixture <- rbind(t(samples.admixture[,1,]), t(samples.admixture[,2,]))
      samples.sumIntensities <- c(samples.sumIntensities[1, ], samples.sumIntensities[2, ])
      for(pop in 1:K) {
        if(AdmixturePrior[1, pop] > 0 & AdmixturePrior[2, pop] > 0) { # bivariate plot
          parents.pop <- kde2d(samples.admixture[, pop], log(samples.sumIntensities), lims=c(0,1,-1,3))
          contour(parents.pop$x, parents.pop$y, parents.pop$z,
                  main=paste("Contour plot of posterior density of parental", population.labels[pop],
                    "admixture proportions"), xlab="Admixture proportion", ylab="Sum-intensities")
          persp(parents.pop$x, parents.pop$y, parents.pop$z, col=popcols[pop],
                main=paste("Perspective plot of bivariate density of parental", population.labels[pop],
                  "admixture proportions"), xlab="Admixture propotion", ylab="log sum-intensities", zlab="Posterior density")
        }
      }
      if(K > 2) {
        for(pop1 in 1:K) {
          for(pop2 in 2:K) {
            if((AdmixturePrior[1, pop1] > 0 & AdmixturePrior[2, pop2] > 0) & (pop1 < pop2)) { # bivariate plot
              parents.pop <- kde2d(samples.admixture[, pop], samples.admixture[, pop], lims=c(0,1,0,1))
              contour(parents.pop$x, parents.pop$y, parents.pop$z,
                      main="Contour plot of posterior density of parental admixture proportions",
                      xlab=paste(population.labels[pop1], "admixture proportion"),
                      ylab=paste(population.labels[pop2], "admixture proportion"))
            }
          }
        }
      }
    }
  } else { # only one parent admixed - univariate plots
    admixedParent <- 0
    if(IsAdmixed[1]) admixedParent <- 1
    if(IsAdmixed[2]) admixedParent <- 2
    if(admixedParent > 0) {
      for(pop in 1:K) {
        if(AdmixturePrior[admixedParent, pop] > 0) {
          plot(density(samples.admixture[pop, admixedParent, ], adjust=0.5, from=0, to=1),
               main=paste(population.labels[pop], "admixture proportion: parent ", admixedParent),
               ylab="Posterior density", xlim=c(0,1))
        }
      }
      plot(density(samples.sumIntensities[admixedParent, ], adjust=0.5, from=1, to=20),
           main=paste("Sum-intensities for parent ", admixedParent),
           ylab="Posterior density", xlim=c(1,20))
      x <- seq(1,20, by=0.1)
      prior <- getSumIntensitiesPrior(user.options)
      ##print(prior)
      points(x, dgamma(x, prior[1], rate=prior[2]), type="l", col="blue")
    }
  }
  dev.off()
}


###################################################################################
## start of script
###################################################################################
options(echo=TRUE)
#par(cex=2,las=1)  # ? this line causes X connection to break
graphics.off()
ps.options(pointsize=16)

##cat("Starting R script\n", file=outfile, append=F)

## read table of user options
message <- c(message, "reading user options...")
##cat("reading user options...", file=outfile, append=T)
user.options <- getUserOptions(paste(resultsdir, "args.txt", sep="/"))
##cat(" done\n", file=outfile, append=T)
message <- c(message, " done\n")

outfile <- paste(resultsdir, user.options$logfile, sep="/")
cat(message, file=outfile, append=T, sep="")
rm(message)

## read table of loci and calculate map positions
loci.compound <- readLoci()
n.chr <- nlevels(factor(loci.compound$Chromosome))

K <- getNumSubpopulations(user.options)
population.labels <- getPopulationLabels(K, user.options)
AdmixturePrior <- getAdmixturePrior(K, user.options)
ParentsIdentified <- getParentsIdentified(AdmixturePrior)
IsAdmixed <- getIsAdmixed(AdmixturePrior)
                                            
param.samples <- NULL
effect.pop <- NULL
pop.admix.prop <- NULL

## read population parameter samples
if(is.null(user.options$paramfile)) {
  cat("no paramfile\n", file=outfile, append=T)
} else {
  paramfile <- paste(resultsdir,user.options$paramfile, sep="/")
  if(!file.exists(paramfile)) {
    cat("paramfile specified but file does not exist\n", file=outfile, append=T)
  } else {
    if(length(scan(paramfile,  what='character', quiet=TRUE)) == 0) {
      cat("paramfile empty\n", file=outfile, append=T)
    } else {
      cat("reading paramfile...", file=outfile, append=T)
      ## param.samples columns contain:    # K Dirichlet parameters 
                                        # global sum of intensities or gamma shape param if hierarchical
      param.samples <- read.table(paramfile, header=TRUE)
      n <- dim(param.samples)[1]
      if(user.options$hapmixmodel ==1){
        checkConvergence(param.samples, "Population sumintensities parameters",
                         paste(resultsdir, "SumIntensitiesConvergenceDiags.txt", sep="/"));
        postscript( paste(resultsdir, "PopSumIntensitiesAutocorrelations.ps", sep="/" ))     
        plotAutocorrelations(param.samples, user.options$every)
        dev.off()
      }else{
        for(pop in 1:K) {##prepend "Dirichlet." to labels of Dirichlet params 
          dimnames(param.samples)[[2]][pop] <- paste("Dirichlet.",
                                                     population.labels[pop], sep="")
        }
        ## Geweke convergence diagnostics,  autocorrelations and ergodic average plots
        checkConvergence(param.samples, "Population admixture parameters",
                         paste(resultsdir, "PopAdmixParamConvergenceDiags.txt", sep="/"));
        postscript( paste(resultsdir, "PopAdmixParamAutocorrelations.ps", sep="/" ))     
        plotAutocorrelations(param.samples, user.options$every)
        dev.off()
        
        if(K > 1) {
          ## extract Dirichlet admixture parameters
          admixparams <- param.samples[, 1:K,drop=FALSE]
          ## calculate population admixture proportions from Dirichlet parameters
          pop.admix.prop <- popAdmixProportions(population.labels, admixparams, K)
        } else {
          pop.admix.prop <- NULL
        }
      }##end nonhapmixmodel block
      cat(" done\n", file=outfile, append=T)
    }
  }
}
  
## read regression parameter samples
if(is.null(user.options$regparamfile) ||
           length(scan(paste(resultsdir, user.options$regparamfile, sep="/"),
                       what='character',quiet=TRUE)) == 0)  {
  cat("No regression paramfile\n", file=outfile, append=T);
  regparam.samples <- NULL
  beta.admixture<-NULL
} else {
  cat("reading regression parameters...", file=outfile, append=T)
  regparam.samples <- read.table(paste(resultsdir, user.options$regparamfile, sep="/"), header=TRUE)
  n.covariates <- getNumCovariates(user.options)
  
  ## Geweke convergence diagnostics, autocorrelation and ergodic average plots
  checkConvergence(regparam.samples, "Regression parameters",
                   paste(resultsdir, "RegressionParamConvergenceDiags.txt", sep="/"))
  postscript(paste(resultsdir, "RegressionParamAutocorrelations.ps", sep="/" ))     
  plotAutocorrelations(regparam.samples, user.options$every)
  dev.off()
  
  #beta.admixture <- getRegressionParamsForAdmixture(user.options, K, n.covariates, population.labels)
  #if(K > 2 && !is.null(pop.admix.prop)) {
    ## calculate estimate of effect of each pop vs all others if there are >2 populations
    #effect.pop <- effectEstimates(beta.admixture, pop.admix.prop, n, K)
  #} else {
  #  effect.pop <- NULL
  #}
  outcome.continuous <- getOutcomeType(dimnames(param.samples)[[2]])  
  ## calculate residual standard deviation
  if(outcome.continuous == 1) {
    residual.SD <- getPrecision(user.options)^-0.5
  }
  cat(" done\n", file=outfile, append=T)
}

## read allele freq dispersion parameter samples
if(is.null(user.options$dispparamfile)||
           length(scan(paste(resultsdir, user.options$dispparamfile, sep="/"),
                       what='character',quiet=TRUE)) == 0)  {
  eta.samples <- NULL
} else {
  cat("reading allele frequency dispersion parameter...", file=outfile, append=T)
  eta.samples<-read.table(paste(resultsdir, user.options$dispparamfile,sep="/"), header=TRUE)
  ## label dispersion parameters
  if(!is.null(user.options$historicallelefreqfile)) {
    dimnames(eta.samples)[[2]] <- paste("eta", population.labels, sep="." )
  }
  checkConvergence(eta.samples, "Dispersion parameters",
                   paste(resultsdir, "DispParamConvergenceDiags.txt", sep="/"))
  postscript(paste(resultsdir, "DispParamAutocorrelations.ps", sep="/" ))     
  plotAutocorrelations(eta.samples, user.options$every)
  dev.off()
  cat(" done\n", file=outfile, append=T)
}   

## combine samples of Dirichlet params, admixture proportions, dispersion params, regression params
param.samples.all <- cbindIfNotNull(param.samples, pop.admix.prop)
param.samples.all <- cbindIfNotNull(param.samples.all, eta.samples)
param.samples.all <- cbindIfNotNull(param.samples.all, regparam.samples)
param.samples.all <- cbindIfNotNull(param.samples.all, effect.pop)
## calculate posterior quantiles
if(!is.null(param.samples.all) && (dim(param.samples.all)[2] > 0)) {
  nvars <- dim(param.samples.all)[2]
  post.quantiles <- calculateAndPlotQuantiles(param.samples.all, nvars)
}

## get population admixture Dirichlet parameters: either posterior means, or values specified in model
if(K == 1) {
  alphas <- c(1)
} else {
  if(user.options$hapmixmodel==1){
    alphas <- rep(1/K, K)
  }else{
    if(!is.null(param.samples)) {
      alphas <- post.quantiles[1:K, 1]
    } else {
      if(!is.null(user.options$admixtureprior)) {
        alphas <- AdmixturePrior[1,]
      } else alphas <- rep(1, K)
    }
  }
}

##plot ergodic averages
if(is.null(user.options$ergodicaveragefile) ) {
  cat("no ergodicaveragefile\n", file=outfile, append=T)
} else {
  if(length(scan(paste(resultsdir,user.options$ergodicaveragefile, sep="/"),  what='character', quiet=TRUE)) == 0) {
    cat("ergodicaveragefile is empty\n", file=outfile, append=T)
  } else {
    cat("plotting ergodic averages...", file=outfile, append=T)
    postscript( paste(resultsdir, "ErgodicAverages.ps", sep="/" ))
    plotErgodicAverages(paste(resultsdir, user.options$ergodicaveragefile, sep="/"), user.options$every)
    dev.off()
    cat(" done\n", file=outfile, append=T)
  }
}

if(!is.null(user.options$thermo) && user.options$thermo == 1){
  cat("plotting energy against coolness...", file=outfile, append=T)
  anneal.table <- read.table(paste(resultsdir, "annealmon.txt", sep="/"), header=TRUE, row.names = NULL)
  postscript(paste(resultsdir, "annealplot.ps", sep="/"))
  ##plot raw points
  plot(anneal.table[,1],anneal.table[,2], xlab="Coolness", ylab="Mean energy", type = 'p')
##  ##fit smoothed spline and overlay on points
##  fit.spline <- smooth.spline(anneal.table[,1], anneal.table[,2], w=-1/anneal.table[,3],spar=0.6)
##  lines(fit.spline, lty=2)
  dev.off()
  cat(" done\n", file=outfile, append=T)
}

#read output of test for heterozygosity and plot
if(!is.null(user.options$hwtestfile) && file.exists(paste(resultsdir,user.options$hwtestfile, sep="/"))){
  cat("plotting scores in test for heterozygosity...", file=outfile, append=T)
  plotHWScoreTest(user.options$hwtestfile, K)
  cat(" done\n", file=outfile, append=T)
}
## read output of score test for allelic association, and plot cumulative results
if(!is.null(user.options$allelicassociationscorefile) && file.exists(paste(resultsdir,user.options$allelicassociationscorefile, sep="/"))) {
  cat("plotting scores in test for allelic association...", file=outfile, append=T)
  outputfilePlot <- paste(resultsdir, "TestsAllelicAssociation.ps", sep="/" )
  plotScoreTest(user.options$allelicassociationscorefile, FALSE, outputfilePlot, user.options$every)
  cat(" done\n", file=outfile, append=T)
}

## read output of score test for association with haplotypes, and plot cumulative results
if(!is.null(user.options$haplotypeassociationscorefile) && file.exists(paste(resultsdir,user.options$haplotypeassociationscorefile, sep="/"))) {
  cat("plotting scores in test for haplotype association...", file=outfile, append=T)
  outputfilePlot <- paste(resultsdir, "TestsHaplotypeAssociation.ps", sep="/" )
  plotScoreTest(user.options$haplotypeassociationscorefile, TRUE, outputfilePlot, user.options$every)
  cat(" done\n", file=outfile, append=T)
}

## read output of regression model score test for ancestry, and plot cumulative results
if(!is.null(user.options$ancestryassociationscorefile) && file.exists(paste(resultsdir,user.options$ancestryassociationscorefile, sep="/"))) {
  cat("plotting scores in test for ancestry association...", file=outfile, append=T)
  ## produces warning
  plotAncestryScoreTest(user.options$ancestryassociationscorefile, "TestsAncestryAssoc",K, population.labels, user.options$every)
  cat(" done\n", file=outfile, append=T)
}

## read output of affecteds-only score test for ancestry, and plot cumulative results
if(!is.null(user.options$affectedsonlyscorefile) && file.exists(paste(resultsdir,user.options$affectedsonlyscorefile, sep="/"))) {
  cat("plotting scores in affecteds-only test...", file=outfile, append=T)
  plotAncestryScoreTest(user.options$affectedsonlyscorefile, "TestsAffectedsOnly",K, population.labels, user.options$every)
  cat(" done\n", file=outfile, append=T)
}

## read output of score test for residual allelic association, and plot cumulative results
if(!is.null(user.options$residualallelicassocscorefile) && file.exists(paste(resultsdir,user.options$residualallelicassocscorefile, sep="/"))) {
  cat("plotting scores in test for residual allelic association", file=outfile, append=T)
  psfile <- paste(resultsdir, "TestsResidualAllelicAssoc.ps", sep="/")
  plotResidualAllelicAssocScoreTest(user.options$residualallelicassocscorefile, psfile, user.options$every)
  cat(" done\n", file=outfile, append=T)
}

if(!is.null(user.options$outcomevarfile) && !is.null(user.options$testgenotypesfile) && file.exists(user.options$testgenotypesfile)){
  cat("performing extra score tests...", file=outfile, append=T)
  ExtraScoreTests(user.options$testgenotypesfile, user.options$outcomevarfile, paste(resultsdir, "ExpectedOutcomes.txt", sep="/"))
  cat(" done\n", file=outfile, append=T)
}

if(is.null(user.options$allelefreqoutputfile) || user.options$fixedallelefreqs==1 || !file.exists( paste(resultsdir,user.options$allelefreqoutputfile, sep="/"))) {
  cat("no allelefreqoutputfile\n", file=outfile, append=T)
} else {
  allelefreq.samples <- dget(paste(resultsdir,user.options$allelefreqoutputfile,sep="/"))
  ## prevent script crashing when an allelefreqoutputfile has been specified with fixed allele frequencies
  if(dim(allelefreq.samples)[2]==0) {
    cat("allelefreqoutputfile is empty\n", file=outfile, append=T)
  } else {
    ## read posterior samples of allele frequencies as 3-way array (pops, alleles within loci,draws)
    cat("Converting allele frequency samples array to list...", file=outfile, append=T)
    allelefreq.samples.list <- listAlleleFreqs(allelefreq.samples)
    ##allelefreq.samples.list <- convertAlleleFreqs(allelefreq.samples)
    cat(" done\n", file=outfile, append=T)
    
    ## calculate posterior means of sampled fvalues at each locus
    if(K > 1) {
      cat("calculating posterior means of sampled f-values...", file=outfile, append=T)
      fValues.means <- calculateLocusfValues(allelefreq.samples.list)
      fValues.means <- data.frame(as.vector(loci.compound[,1]), round(fValues.means, digits=4))
      dimnames(fValues.means)[[2]] <- c("LocusName",
                                        paste(population.labels[1], population.labels[2], sep="."))
      write.table(fValues.means, file=paste(resultsdir,"LocusfValues.txt", sep="/"),
                  row.names=TRUE, col.names=TRUE)
      write.table(fValues.means[order(fValues.means[, 2], decreasing=TRUE), ],
                  file=paste(resultsdir,"LocusfValuesSorted.txt", sep="/"),
                  row.names=FALSE, col.names=TRUE)
      cat(" done\n", file=outfile, append=T)
    }

    ## calculate posterior means of KL info for ancestry at each locus
    if(K > 1) {
      cat("calculating posterior means ok KL info...", file=outfile, append=T)
      KLInfo.means <- calculateLocusKLInfo(allelefreq.samples.list)
      KLInfo.means <- data.frame(as.vector(loci.compound[,1]), round(KLInfo.means, digits=4))
      dimnames(KLInfo.means)[[2]] <- c("LocusName", "KLInfo")
      write.table(KLInfo.means, file=paste(resultsdir,"LocusKLInfo.txt", sep="/"),
                  row.names=TRUE, col.names=TRUE)
      write.table(KLInfo.means[order(KLInfo.means[, 2], decreasing=TRUE), ],
                  file=paste(resultsdir,"LocusKLInfoSorted.txt", sep="/"),
                  row.names=FALSE, col.names=TRUE)
      cat(" done\n", file=outfile, append=T)
    }

    cat("writing posterior means of allele frequencies...", file=outfile, append=T)
    ## generate lists to hold allele freq means and covariances
    freqMeansCovs <- listFreqMeansCovs(allelefreq.samples.list)
    allelefreq.means.list <- freqMeansCovs[[1]]
    allelefreq.covs.list <- freqMeansCovs[[2]]

    ## write posterior means of allele freqs to file
    writeAlleleFreqs(allelefreq.means.list, K, loci.compound, population.labels,
                     paste(resultsdir, "AlleleFreqPosteriorMeans.txt", sep="/" ))
    
    ## fit Dirichlet parameters by equating posterior means and variances
    allelefreq.params.list <- fitDirichletParams(allelefreq.means.list, allelefreq.covs.list) 
    
    ## write Dirichlet parameters of allele freq distributions to file as R object
    dput(allelefreq.params.list, file=paste(resultsdir,"allelefreqparamsAsRObject.txt", sep="/"))
    
    ## write Dirichlet parameters of allele freq distributions to file
    writeAlleleFreqs(allelefreq.params.list, K, loci.compound, population.labels,
                     paste(resultsdir, "AlleleFreqPosteriorParams.txt", sep="/" ))
    cat(" done\n", file=outfile, append=T)

  }
}

  
if(!is.null(user.options$indadmixturefile) && K >1 && file.exists(paste(resultsdir, user.options$indadmixturefile, sep="/"))) {
  ## read posterior samples of individual admixture
  cat("reading posterior samples of individual params...", file=outfile, append=T)
  samples <- dget(paste(resultsdir,user.options$indadmixturefile,sep="/"))
  n.iterations <- dim(samples)[3]
  n.individuals <- dim(samples)[2]
  ## reformat as 4-way array if random mating model
  if(!is.null(user.options$randommatingmodel) && user.options$randommatingmodel==1) {
    ## drop any extra vars in the array
    samples.adm <- samples[1:(2*K),1:n.individuals ,1:n.iterations, drop=F]
    samples.sumintensities <- samples[-c(1:(2*K)),1:n.individuals ,1:n.iterations, drop=F]
    samples4way <- array(samples.adm, dim=c(K, 2, dim(samples)[2:3]))
    dimnames(samples4way) <- list(population.labels,c("Parent1", "Parent2"),
                                  character(0), character(0))
    samples.meanparents <- apply(samples4way, 2:4, mean)
    samples.bothparents <- array(samples4way,
                                 dim=c(K, 2*n.individuals, n.iterations))
    dimnames(samples.bothparents) <- list(population.labels, character(0), character(0))
  } else {
    samples.meanparents <- samples[1:K, 1:n.individuals ,1:n.iterations, drop=F]
    samples.bothparents <- samples[1:K, 1:n.individuals ,1:n.iterations, drop=F]
  }
  cat(" done\n", file=outfile, append=T)
  ## should plot only if subpopulations are identifiable in model
  if(n.individuals > 1) { # dim(samples.meanparents)[2] > 1) {
    cat("plotting posterior distribution of admixture...", file=outfile, append=T)
    plotAdmixtureDistribution(alphas, samples.bothparents, K)
    cat(" done\n", file=outfile, append=T)
    ##if(K > 1) {
    ##writePosteriorMeansIndivAdmixture(t(samples), K)
    ##}
    cat("writing posterior means of individual admixture...", file=outfile, append=T)
    sample.means <- apply(samples, 1:2, mean)
    write.table(format(round(t(sample.means),3), nsmall=3),
                paste(resultsdir, "IndAdmixPosteriorMeans.txt", sep="/"), quote=F, row.names=F, sep="\t")
    cat(" done\n", file=outfile, append=T)
    
  } else { #if(dim(samples.meanparents)[2]==1)
    if(IsAdmixed[1] | IsAdmixed[2]) {
      cat("plotting posterior distribution of individualadmixture...", file=outfile, append=T)
      ## convert samples4way to 3way array
      samplesOneIndiv <- samples4way[, , 1, ]
      samples.sumintensities <- samples.sumintensities[, 1, ]
      plotPosteriorDensityIndivParameters(samplesOneIndiv, samples.sumintensities, user.options,
                                          population.labels, AdmixturePrior, IsAdmixed, ParentsIdentified)
      cat(" done\n", file=outfile, append=T)
    }
  }
}else{
cat("No indadmixturefile\n", file=outfile, append=T)
}
cat("Script completed", file=outfile, append=T)
