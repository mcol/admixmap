
##read posterior conditional genotype probs
CGProbs<-dget("PPGenotypeProbs.txt")
##get expected number of copies of allele 2
Ex <- CGProbs[2,,] + 2*CGProbs[3,,]
Exsq <- CGProbs[2,,] + 4*CGProbs[3,,]
Vx <- Exsq - Ex*Ex

L <- dim(CGProbs)[3]
rm(CGProbs)

EY <- as.matrix(dget("ExpectedOutcomes.txt")[,1,])
EY <- apply(EY, 1 ,mean)
Y <- Y <- as.vector(as.numeric(read.table("helen_outcome.txt", header=T, colClasses="numeric")[,1]))
YMinusEY <- Y - EY
rm(Y)

N <- length(EY)

Score <- rep(0,L)
CompleteInfo <- rep(0,L)
ScoreSq <- rep(0,L)

#for( locus in 1:L){

for(indiv in 1:N){
  Score <- Score + YMinusEY[indiv]*Ex[indiv, ]
  CompleteInfo <- CompleteInfo + EY[indiv]*(1-EY[indiv])*Exsq[indiv,]
  ScoreSq <- ScoreSq + YMinusEY[indiv]*YMinusEY[indiv]*Exsq[indiv, ]
}
#}

MissingInfo <- ScoreSq - Score*Score
ObservedInfo = CompleteInfo - MissingInfo

