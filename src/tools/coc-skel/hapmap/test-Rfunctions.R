 
bir2 <- info.reward2( prior=c(0.7, 0.2, 0.1), predictive=c(0.8, 0.1, 0.1), t=1)
print(bir2)

bir2 <- info.reward2( prior=c(0.7, 0.2, 0.1), predictive=c(1, 0, 0), t=1)
print(bir2)

# for each locus we have: allele freqs p and q in sample from which to calculate
## prior genotype probs as p^2, 2pq, q^2
## (strictly we should marginalize over the posterior distribution of allele freqs
## not sure how to do this on the assumption of HWE - it's not simply a multinomial-Dirichlet


## for each masked genotype we have: predictive probs output by imputation program
## true masked genotype

