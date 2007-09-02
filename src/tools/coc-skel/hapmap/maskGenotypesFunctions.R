# Functions file for maskGenotypes.R and other R scripts.
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
get.option <- function(optname, s, int = FALSE) {
        opt <- unlist(strsplit(s, "="))
        if (opt[1] != optname && opt[1] != paste("--", optname, sep = "")) {
                stop("Expected --", optname, " option name, got ", opt[1])
        }
##        print(optname);
##        print(opt);
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
	# Wikipedia,
	# http://en.wikipedia.org/wiki/Information_entropy
	# cites use of log2 instead of log (natural logarithm).
	# return(sum(-vec * log(vec), na.rm = na.rm))
	return(sum(-vec * log2(vec), na.rm = na.rm))
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

sharpen.vector <- function(vec) {
	# print(vec)
	# # When vec has NA values
	# vec = na.omit(vec)
	# if (length(vec) == 0) {
	# 	return(c(NA, NA, NA))
	# }
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

# Calculate the mutual information as given in the paper, page 18.
# \sum_{x,y} (p_{x,y} log p_{x,y} - p_x p_y log(p_{x.} p_{.y}))
# Uses the entropy() function.

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

coc.locus <- function (locus.no, GP, orig, genotypes) {
	joint.pd <- joint.prob(GP, locus.no, genotypes)
	joint.pd.no.uncert <- joint.prob(GP, locus.no, genotypes, throw_away_uncertainity = TRUE)
	# mi[locus.no, 1] <- coefficient.of.constraint(joint.pd)
	# mi[locus.no, 2] <- coefficient.of.constraint(joint.pd.no.uncert)
	# mi[locus.no, 3] <- length(unique(na.omit(orig[, locus.no])))
	return(c(
		coefficient.of.constraint(joint.pd),
		coefficient.of.constraint(joint.pd.no.uncert),
		length(unique(na.omit(orig[, locus.no])))))
}

# Return number of the genotype, based by the genotype.
# There probably is a better way to implement it using a dictionary.
num.by.genotype <- function (genotypes, genotype) {
  for (idx in 1:length(genotypes)) {
    if (genotype == genotypes[idx]) {
      return(idx)
    }
  }
}

get.coc.table <- function (GP, orig, genotypes) {
# NaN problem: when GP has some NaN values, we need to skip them.
# Indices of loci with NaN
# All locus indices
# Iterate over `forbidden' indices and remove the data. I'm almost sure
# that there is a non-for-looped way to calculate a difference between
# two sets in R, but I couldn't find it.
	nan_idx <- as.vector(which(is.na(GP[1,,1])))
	# Array for mutual information.
	mi <- matrix(nrow = length(GP[1,,1]), ncol = 3)
	colnames(mi) <- c(
		"Mutual information",
		"Mutual information, no uncertainity",
		"Unique genotypes")

	for (locus.no in 1:length(GP[1,,1])) {
		if (length(nan_idx) > 0) {
			# FIXME: This is an ad-hoc fix that only works for 1 NA
			# locus. If there are more NA loci, it will fail.
			if (locus.no == nan_idx[1]) {
				next
			}
		}
		# print(c("Number of unique genotypes: ", length(unique(na.omit(orig[, locus.no])))))
		joint.pd <- joint.prob(GP, locus.no, genotypes)
		joint.pd.no.uncert <- joint.prob(GP, locus.no, genotypes, throw_away_uncertainity = TRUE)
		# mi[locus.no, ] <- coc_locus(locus.no, GP, orig, genotypes)
		mi.tmp <- coc_locus(locus.no, GP, orig, genotypes)
		mi[locus.no, 1] <- mi.tmp[1]
		mi[locus.no, 2] <- mi.tmp[2]
		mi[locus.no, 3] <- mi.tmp[3]
		# print(as.matrix(cbind(c("Locus no.", locus.no), as.vector(mi[locus.no, ]))))
		# print((cbind(t(c("Locus no.", locus.no)), t(as.vector(mi[locus.no, ])))))
	}
	return(mi)
}

coc.dbg <- function (locus.no, GP, orig, genotypes) {
	# for (indiv.no in 1:length(GP[1, locus.no, ])) {
		# print(indiv.no)
		# print(rbind(GP[, locus.no, indiv.no], c("**", "", "")))
		# print(orig[indiv.no, locus.no])
		# print(cbind(GP[, locus.no, indiv.no], orig[indiv.no, locus.no]))
	# }
	print(rbind(GP[, locus.no, ], orig[, locus.no]))
	joint.pd <- joint.prob(GP, locus.no, genotypes)
	print("Joint probabilities")
	print(joint.pd)
	print("Marginal distribution product")
	print(marginal.dist.prod(joint.pd))
	print("Coefficient of constraint")
	print(coc_locus(locus.no, GP, orig, genotypes))
}


info.reward2 <- function(prior, predictive, t) {
   ## prior is a probability vector of length K
   ## predictive is a probability vector of length K
   ## t is the correct class: integer between 1 and K

   # prior should be calculated from the allele frequencies in the 50
   # (or whatever) hapmap individuals whose genotypes were not masked.
   # If the observed proportions of allele 1 and 2 on these 100 gametes
   # are p and q, the prior should be a vector with values p^2, 2*p*q,
   # and q^2.

   K <- length(prior)
   ## normalize prior probability vector
   prior <- prior / sum(prior)
   if (sum(predictive) == 0) { ## no prediction, return zero
     i.reward <- 0
   } else {
     predictive <- predictive / sum(predictive)
     ## i.plus evaluates to -Inf if predictive[t] is zero
     i.plus <- log(predictive[t] / prior[t])
     ## i.minus evaluates to -Inf if any element of predictive[-t] is 1
     i.minus <- sum( log( (1 - predictive[-t]) / (1 - prior[-t]) ) )
     i.reward <- (i.plus + i.minus) / K
   }
   return(i.reward)
}

bir.locus <- function (locus.no, GP, orig, genotypes, prior) {
##  tmp.prior <- prior[locus.no, ]
  ## print(dim(GP))
  ## print(c("GP"))
  ## print(GP[,locus.no,])
  ## print(dim(orig))
  ## print(orig[,locus.no])
  ## print(genotypes)

  ## print(prior[locus.no,])
  indivs <- length(orig[,locus.no])
  ## print(indivs)
  tmp.bir <- rep(0, length(indivs))
  for (idx in 1:indivs) {
    ## print(c(idx,
    ## orig[idx,locus.no],
    ## num.by.genotype(genotypes, orig[idx,locus.no]),
    ## GP[,locus.no,idx], prior[locus.no, ]))
    tmp.predictive <- GP[, locus.no, idx]
    correct.class <- num.by.genotype(genotypes, orig[idx,locus.no])
    tmp.bir[idx] <- info.reward2(tmp.prior, tmp.predictive, correct.class)
  }
  return(c(
           mean(tmp.bir),
           0,
           length(unique(na.omit(orig[, locus.no])))))
}

bir.locus2 <- function (locus.name, GP, obs, prior) {
  ## locus.name is the name of the locus, used to index the arrays (eg "rs123456")
  ## if there are L loci, N individuals (with masked data) and K states, then
  ## predicted is an L*N*K array of posterior predictive probabilities
  ## obs is an L*K array of observed values
  ## prior is an L*K array of prior predictive probabilities
  
  ##indivs <- length(obs[,locus.name])
  indivs <- dim(GP)[3]

  tmp.bir <- rep(0, length(indivs))
  for (idx in 1:indivs) {
    ## print(c(idx,
    ## orig[idx,locus.no],
    ## num.by.genotype(genotypes, orig[idx,locus.no]),
    ## GP[,locus.name,idx], prior[locus.name, ]))

    tmp.bir[idx] <- info.reward2(prior[locus.name,], GP[, locus.name, idx], obs[idx,locus.name])
  }
  return(c(
           mean(tmp.bir),
           0,
           length(unique(na.omit(obs[, locus.name])))))
}
