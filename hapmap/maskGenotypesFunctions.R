## Functions file for maskGenotypes.R and other R scripts.
# Returns a file name derived from the path and the base name
get.io.filename <- function(extension) {
	return(paste(
		io.file.path,
		"/",
		io.file.basename,
		extension,
		sep = ""))
}

# Saves the table with our custom settings
# It assumes that rows have labels stored not as a simple column, but as
# a row names in the manner that read.table handles row.names
# parameter. Hence row.names = TRUE.
genepi.write.table <- function(obj, out.file.name) {
	write.table(
		obj,
		file = out.file.name,
		row.names = TRUE,
		col.names = NA,
		sep = "\t",
		quote = TRUE,
		na = missing.genotype)
}

# For getting options from the command-line. Not exactly like getopt,
# though...
get_option <- function(optname, s, int = FALSE) {
        opt <- unlist(strsplit(s, "="))
        if (opt[1] != paste("--", optname, sep = "")) {
                stop("Expected --", optname, " option name, got ", opt[1])
        }
        print(optname);
        print(opt);
        if (int) {
                return(as.integer(opt[2]))
        } else {
                return(opt[2])
        }
}

# Mutual Information functions come here as well.

# Evaluate entropy of a vector
entropy <- function(vec, na.rm = FALSE) {
	if (!na.rm) {
		na.fail(vec)
	}
	# This check is disabled because it false-positives this data:
	#
	#      t1,1 t1,2         t2,2
	#      p1,1    0    0 1.488159e-19
	#      p1,2    0    0 7.661409e-10
	#      p2,2    0    0 1.000000e+00
	#
	# if (max(vec, na.rm = TRUE) > 1) {
	# 	print(vec)
	# 	stop("argument elements must be at most 1.")
	# }
	if (min(vec, na.rm = TRUE) < 0) {
		stop("argument elements must be greater than 0.")
	}
	vecsum <- sum(vec, na.rm = na.rm)
	# Since it will never sum to exactly 1, there is certain precision
	if (abs(vecsum - 1.0) > 1e-5) {
		print(vec)
		stop("argument needs to sum to 1. Current sum is ", vecsum)
	}
	# Zeros do not contribute to entropy
	# H(a_1, a_2, ..., 0) = H(a_1, a_2, ...)
	vec <- vec[which(vec != 0)]
	return(sum(-vec * log(vec), na.rm = na.rm))
}

# Returns a function that adds a prefix to a string
# NB: It doesn't return a string, but a function, which can be later
# called with an argument.
#
# Usage:
# > add.something <- add.prefix("something ")
# > add.something("strange")
# [1] "something strange"
# > add.something("nice")
# [1] "something nice"
#
add.prefix <- function(s) {
	return(function(x) (paste(s, x, sep = "")))
}

get.clean.cm <- function() {
	m = matrix(numeric(9), nrow = 3, ncol = 3)
	# t: true
	# c("t1,1", "t1,2", "t2,2")
	colnames(m) <- lapply(genotypes, add.prefix("t"))
	# p: predicted
	# c("p1,1", "p1,2", "p2,2")
	rownames(m) <- lapply(genotypes, add.prefix("p"))
	# It should look like this:
	# > m
	#           t1,1 t1,2 t2,2
	#      p1,1    0    0    0
	#      p1,2    0    0    0
	#      p2,2    0    0    0
	return(m)
}

# Calculates a product of marginal distributions
marginal.dist.prod <- function(my.matrix) {
	if (!is.matrix(my.matrix)) {
		stop("Argument should be a matrix.")
	}
	mcols = colSums(my.matrix, na.rm = TRUE)
	# mcols = mcols / sum(mcols)
	mrows = rowSums(my.matrix, na.rm = TRUE)
	# mrows = mrows / sum(mrows)
	# marginal distributions product
	md.prod = matrix(mrows, ncol = 1) %*% matrix(mcols, nrow = 1)
	# md.prod[which(md.prod == 0)] <- NA
	dimnames(md.prod) <- dimnames(my.matrix)
	# print(md.prod)
	return(md.prod)
}

# Calculate the mutual information as given in the paper, page 18.
# \sum_{x,y} (p_{x,y} log p_{x,y} - p_x p_y log(p_{x.} p_{.y}))
# Uses the entropy() function.

sharpen.vector <- function(vec) {
	# print(vec)
	m = max(vec)
	maxs_idx = which(vec == m)
	if (length(maxs_idx) != 1) {
		# There is more than one genotype with
		# probability equal to the max probability.
		# warning("Can't find max element.", ", ", indiv, ", ", paste(vec, collapse = " "))
		warning("Can't find max element. ", paste(vec, collapse = " "))
		# In such a case, we will throw away all but the
		# first probability.
		maxs_idx <- maxs_idx[1]
	}
	vec[] = 0
	vec[maxs_idx] = 1
	return(vec)
}

# Finds the most probable genotypes and assigns p = 1 to them, setting
# all the other genotypes to zero. Needed for fair comparison with
# fastPHASE which doesn't output any uncertainity.
sharpen.probs<- function(a) {
	# print("dim(a)")
	# print(dim(a))
	if (dim(a)[2] == 0) {
		# Don't do anything if it's empty. Trying to do
		# something with it would trigger an error.
		return(a)
	}
	for (indiv in 1:dim(a)[2]) {
		# print(c("Indiv: ", indiv))
		# print( a[, indiv] )
		a[, indiv] <- sharpen.vector(a[, indiv])
	}       
	return(a)
}

joint.prob <- function(GP, locus.no, genotypes, throw_away_uncertainity = FALSE) {
	joint.pd <- get.clean.cm()
	indiv.total <- 0
	for (genotype in genotypes) {
		# message(c("Genotype start:", genotype))
		indiv.idx <- which(orig[, locus.no] == genotype)
		indiv.no <- length(indiv.idx)
		# Adding to total only individuals that have data
		indiv.total <- indiv.total + indiv.no
		# print(rowMeans(matrix(GP[, indiv.idx, locus.no], nrow = 3)))
		colname <- add.prefix("t")(genotype)
		# message(c("colname:", colname))
		# print(indiv.idx)
		# print(locus.no)
		# print(matrix(GP[, locus.no, indiv.idx], nrow = 3))
		indiv.probs <- matrix(GP[, locus.no, indiv.idx], nrow = 3)
		#     in1 in2 in3 in4 ...
		# 1,1 0.1 0.2
		# 1,2 0.2 0.3
		# 2,2 0.7 0.5
		# print(indiv.probs)
		if (throw_away_uncertainity) {
			indiv.probs = sharpen.probs(indiv.probs)
		}
		if (length(indiv.probs) != 0) {
			# Sum over the individuals
			joint.pd[, colname] <- (rowSums(indiv.probs))
		}
		# print(joint.pd)
		# message(c("Genotype end:", genotype))
	}
	# Make it sum to one alltogether.
	joint.pd <- (joint.pd / indiv.total)
	return(joint.pd)
}

matrix.debug <- function(m, msg = "Matrix debug") {
	print(msg)
	print(dim(m))
	print(m)
	print(colSums(m))
	print(sum(m))
	print(entropy(m))
}

# Mutual information, as defined in the Wikipedia article.
# Source:
# http://en.wikipedia.org/wiki/Mutual_information#Relation_to_other_quantities
mutual.information <- function(joint.pd) {
	# marginal distributions
  if (!is.matrix(joint.pd)) {
    stop("Argument should be a matrix.")
  }
	# TODO: Check if the matrix is square
  if (dim(joint.pd)[1] != dim(joint.pd)[2]) {
    stop("Matrix should be square.")
  }
  x.entropy = entropy(rowSums(joint.pd))
  y.entropy = entropy(colSums(joint.pd))
  joint.pd.entropy = entropy(joint.pd)
  mi = x.entropy + y.entropy - joint.pd.entropy
	# Probably because of http://actin.ucd.ie/trac/genepi/ticket/2
	# If the result is very close to zero, round it to zero.
  if (abs(mi) < 1e-7) {
    return(0)
  } else {
    return(mi)
  }
}

# After Wikipedia:
# http://en.wikipedia.org/wiki/Mutual_information#Normalized_variants
coefficient.of.constraint <- function(joint.pd) {
	mi <- mutual.information(joint.pd)
	Hy = entropy(colSums(joint.pd))
	if (Hy == 0) {
		# print(joint.pd)
		# stop("Zero entropy!")
		#
		# This usually happens when all individuals from the
		# case-control set are monomorphic.
		return(NA)
	} else {
		return(mi / Hy)
	}
}

##source: Dawy et al ()
##a classification measure
DCL <- function(joint.pd){
## marginal distributions
  if (!is.matrix(joint.pd)) {
    stop("Argument should be a matrix.")
  }
## TODO: Check if the matrix is square
  if (dim(joint.pd)[1] != dim(joint.pd)[2]) {
    stop("Matrix should be square.")
  }
  x.entropy = entropy(rowSums(joint.pd))
  y.entropy = entropy(colSums(joint.pd))
  joint.pd.entropy = entropy(joint.pd)
  mi = x.entropy + y.entropy - joint.pd.entropy
  
  dcl <- 1 - mi / max(x.entropy, y.entropy)
  return(dcl)  
}
entropy<-function(x){
  return(-sum(x[x>0]*log(x[x>0])))
}

##calculates the average information score according to
##Kononenko & Bratko, Machine Learning 6, 67-80 (1991)
InfoScore <- function(prior, posterior, observed){
  ##prior = vector of K probabilities
  ##posterior = NxK matrix of posterior predictive probabilities
  ##observed = vector of N observed values
  N <- nrow(posterior)
  K <- ncol(posterior)

  if(length(prior)!=K || length(observed)!=N){
    stop("mismatched dimensions")
  }
  sum.info.score <- 0
  for( i in 1:N){
    if(posterior[i, observed[i]] >= prior[observed[i]]){
      sum.info.score <- sum.info.score + log2(posterior[i, observed[i]]) - log2(prior[observed[i]])
    }else{
      sum.info.score <- sum.info.score + log2(1 - prior[observed[i]]) - log2(1 - posterior[i, observed[i]])
    }
  }
  return(sum.info.score / K)
  
}

info.score.denom <- function(prior, posterior, observed){
  N <- nrow(posterior)
  K <- ncol(posterior)

  if(length(prior)!=K || length(observed)!=N){
    stop("mismatched dimensions")
  }
  sum <- 0
  for( i in 1:N){
    ##if(posterior[i, observed[i]] >= prior[observed[i]]){
      sum <- sum - log2(prior[observed[i]])
    ##}
  }
  return(sum)
}

Relative.Info.Score <- function(prior, posterior, observed){
  ##joint.entropy <- sum(apply(posterior, 1, entropy))
 denom <- info.score.denom(prior, posterior, observed)
  return (100*InfoScore(prior, posterior, observed) / denom)
}

RIS.locus <- function(locus, prior, posterior, observed){
  RIS.average <- Relative.Info.Score(prior, t(posterior[5,,locus,]), observed[,locus])
  RIS.seeds <- c(Relative.Info.Score(prior, t(posterior[1,,locus,]), observed[,locus]),
                 Relative.Info.Score(prior, t(posterior[2,,locus,]), observed[,locus]),
                 Relative.Info.Score(prior, t(posterior[3,,locus,]), observed[,locus]),
                 Relative.Info.Score(prior, t(posterior[4,,locus,]), observed[,locus])
                 )
  return(list(RIS.average, mean(RIS.seeds)))
}
