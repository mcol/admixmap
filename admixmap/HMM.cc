#include "rand.h"
#include "HMM.h"
#include "functions.h"

using namespace std;

HMM::HMM()
{
  alpha = 0;
  beta = 0;
  LambdaBeta = 0;
  p = 0;
  StateArrivalProbs = 0;
  ThetaThetaPrime = 0;
  colProb = 0;
  Expectation0 = 0;
  Expectation1 = 0;
  rowProb = 0;
  cov = 0;
  f = 0;
  Lambda = 0;
  betaAllocated = false;
  alphaIsBad = true;
  betaIsBad = true;
}

//not currently used
//HMM objects are instantiated in Chromosome using default constructor above
//and dimensions set by SetDimensions below
HMM::HMM( int inTransitions, int pops) {
  SetDimensions(inTransitions, pops);
}

HMM::~HMM()
{
  delete[] p;
  delete[] alpha;
  if(betaAllocated) {// delete beta array only if allocated with new
    delete[] beta;
    delete[] LambdaBeta;
  }
  delete[] StateArrivalProbs;
  delete[] ThetaThetaPrime;
  delete[] rowProb;
  delete[] colProb;
  delete[] Expectation0;
  delete[] Expectation1;
  delete[] cov; //free_matrix(cov, K);
}

void HMM::SetDimensions( int inTransitions, int pops)
{
  //inTransitions = #transitions +1 = #Loci 
  //pops = #populations
  K = pops;
  DStates = K*K;
  Transitions = inTransitions;
  alpha = new double[Transitions*K*K];
  sumfactor=0.0;
  p = new double[Transitions];
  StateArrivalProbs = new double[Transitions * K * 2];
  ThetaThetaPrime = new double[K*K];
  if(K>2){
    rowProb = new double[K];
    colProb = new double[K];
    Expectation0 = new double[K];
    Expectation1 = new double[K];
    cov = new double[K*K]; 
  }
}

void HMM::SetGenotypeProbs(const double* const lambdain, const bool* const missing){
  Lambda = lambdain;
  missingGenotypes = missing;
  alphaIsBad = true;//new input so reset
  betaIsBad = true;
}


void HMM::SetStateArrivalProbs(const double* const fin, const double* const Theta, const int Mcol, 
			       const bool isdiploid){
  f = fin;
  theta = Theta;
  alphaIsBad = true;//new input so reset
  betaIsBad = true;
  if(isdiploid){
    for(int t = 1; t < Transitions; t++ ){        
      for(int j = 0; j < K; ++j){
	StateArrivalProbs[t*K*2 + j*2]    = (1.0 - f[2*t]) * Theta[j];
	StateArrivalProbs[t*K*2 + j*2 +1] = (1.0 - f[2*t + 1]) * Theta[K*Mcol +j ];
      }
      p[t] = f[2*t] * f[2*t + 1];
    }
    for(int j0 = 0; j0 < K; ++j0) {
      for(int j1 = 0; j1 < K; ++j1) {
	ThetaThetaPrime[j0*K + j1] = Theta[j0]*Theta[j1 + K*Mcol];
      }
    }
  }
}

void HMM::SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, 
			       int *SumLocusAncestry, vector<unsigned> &SumNumArrivals, bool SampleArrivals, unsigned startlocus)const {
  //this does not require forward or backward probs, just state arrival probs
  if(!Lambda || !theta || !f)throw string("Error: Call to HMM::SampleJumpIndicators when StateArrivalProbs are not set!");
  bool xi = false;//jump indicator
  double ProbJump; // prob jump indicator is 1
  // first locus not included in loop below
  for( unsigned int g = 0; g < gametes; g++ ){
    SumLocusAncestry[ g*K + LocusAncestry[g*Transitions] ]++;
  }
  for( int t = 1; t < Transitions; t++ ) {
    for( unsigned int g = 0; g < gametes; g++ ){
      xi = true;
      if( LocusAncestry[g*Transitions + t-1] == LocusAncestry[g*Transitions + t] ){
	ProbJump = StateArrivalProbs[t*K*2 +LocusAncestry[t + g*Transitions]*2 + g];  
	xi = (bool)(ProbJump / (ProbJump + f[2*t+g]) > Rand::myrand());
      } 
      if( xi ){ // increment sumlocusancestry if jump indicator is 1
	SumLocusAncestry[ g*K + LocusAncestry[t+g*Transitions] ]++;
	if(SampleArrivals) { // sample number of arrivals where jump indicator is 1
	  double u = Rand::myrand();
	  // sample distance dlast back to last arrival, as dlast = -log[1 - u(1-f)] / rho
	  // then sample number of arrivals before last as Poisson( rho*(d - dlast) )
	  // algorithm does not require rho or d, only u and f
	  unsigned int sample = Rand::genpoi( log( (1 - u*( 1 - f[2*t+g])) / f[2*t+g] ) );
	  SumNumArrivals[2*(startlocus + t) + g] += sample + 1;
	}
      }//end if xi true
    }//end gamete loop
  } // ends loop over intervals
}

/*
  returns log-likelihood
  This is the sum over states of products of alpha and beta
  and is the same for all t so it is convenient to compute for
  t = T-1 since beta here is 1 so no backward recursions are required. 
*/
double HMM::getLogLikelihood(bool isdiploid) 
{ 
  if(alphaIsBad){
    if(isdiploid)UpdateForwardProbsDiploid();
    else UpdateForwardProbsHaploid();
    }
  double sum = 0;
  if(isdiploid) {
    for( int j = 0; j < DStates; j++ ) {
      sum += alpha[(Transitions - 1)*DStates + j];
    } 
  } else {//haploid
    for( int j = 0; j < K; j++ ) {
      sum += alpha[(Transitions - 1)*K + j];
    }
  }
  return( sumfactor+log( sum ) );
}
      
/*
  Samples Hidden States
  ---------------------
  SStates          - an int array to store the sampled states
  isdiploid  - indicator for diploidy
*/
void HMM::Sample(int *SStates, const bool isdiploid)
{
  if(alphaIsBad){
    if(isdiploid)UpdateForwardProbsDiploid();
    else UpdateForwardProbsHaploid();
  }
  int j1,j2;
  if(isdiploid) { 
    double* V = new double[DStates]; //probability vector for possible states (haploid or diploid)
    int C = 0; // sampled state (haploid or diploid) coded as integer
    // array Sstates: elements 0 to T-1 represent paternal gamete, elements T to 2T-1 represent maternal gamete 
    int State = 0;
    // sample rightmost locus  
    for( int j = 0; j < DStates; ++j ) V[State++] = alpha[(Transitions - 1)*DStates + j];
    
    C = Rand::SampleFromDiscrete( V, DStates ); 
    SStates[Transitions-1] = (int)(C/K);
    SStates[Transitions - 1 + Transitions] = (C % K);
    
    for( int t =  Transitions - 2; t >= 0; t-- ) { // loop from right to left
      j1 = SStates[t+1];                     // ancestry on gamete 0 at locus t+1
      j2 = SStates[Transitions + t + 1];     // ancestry on gamete 1 at locus t+1
      State = 0;
      for(int i1 = 0; i1 < K; ++i1)for(int i2 = 0; i2 < K; ++i2) {
	V[State] = 
	  ( (i1==j1)*f[2*t+2] + StateArrivalProbs[(t+1)*K*2 + j1*2] ) * 
	  ( (i2==j2)*f[2*t+3] + StateArrivalProbs[(t+1)*K*2 + j2*2 +1] );
	V[State] *= alpha[t*DStates + i1*K + i2];
	++State;
      }
      C = Rand::SampleFromDiscrete( V, DStates );
      SStates[t] = (int)(C/K);//paternal
      SStates[t + Transitions] = (C % K);//maternal
    }
    delete[] V;
  } else {//haploid
    double* V = new double[K]; //probability vector for possible states (haploid or diploid)
    int* C = new int[Transitions]; // sampled state (haploid or diploid) coded as integer
    for( int j = 0; j < K; j++ )V[j] = alpha[(Transitions - 1)*K + j ];
    C[ Transitions - 1 ] = Rand::SampleFromDiscrete( V, K );
    SStates[Transitions-1] = C[Transitions-1];
    for( int t =  Transitions - 2; t >= 0; t-- ){
      for(int j = 0; j < K; j++)V[j] = (j == C[t+1])*f[2*t+1]+theta[C[t+1]]*(1.0 - f[2*t]);
      for( int j = 0; j < K; j++ )	V[j] *= alpha[t*K + j];
      C[ t ] = Rand::SampleFromDiscrete( V, K );
      SStates[t] = C[t];
    }
    delete[] C;
    delete[] V;
  }
}
/// Returns a vector of conditional probabilities of each hidden state at 'time' t.
/// If diploid, the vector is really a matrix of probabilities of pairs of states.
const std::vector<double> HMM::Get3WayStateProbs( const bool isDiploid, int t){
  if(alphaIsBad){
    if(isDiploid)UpdateForwardProbsDiploid();
    else UpdateForwardProbsHaploid();
  }
  if(betaIsBad){
    if(isDiploid)UpdateBackwardProbsDiploid();
    else UpdateBackwardProbsHaploid();
  }

  double sum = 0.0;
  // int State = 0;
  int States = 0;
  if(isDiploid) {
    States=DStates;
  } else {
    States=K;
  }
  std::vector<double> probs(States);

  for( int j = 0; j < States; j++ ){
    probs[j] = alpha[t*States + j] * beta[t*States + j];
    sum += probs[j];
  }

  for( int j = 0; j < States; j++ ){ //vectorization successful
     probs[j] /= sum;
   }

  return probs;
}

// ****** End Public Interface *******

/*
  Updates Forward probabilities alpha, diploid case only
*/
void HMM::UpdateForwardProbsDiploid()
{
  if(!Lambda || !theta || !f)throw string("Error: Call to HMM when inputs are not set!");
  // if genotypes missing at locus, skip multiplication by lambda and scaling at next locus   
  sumfactor = 0.0; // accumulates log-likelihood
  double Sum = 0.0;
  double scaleFactor = 0.0;
  
  if(!missingGenotypes[0]) {
    for(int j = 0; j < DStates; ++j) {
      alpha[j] =  ThetaThetaPrime[j] * Lambda[j];
    } 
  } else {
    for(int j = 0; j < DStates; ++j) {
      alpha[j] = ThetaThetaPrime[j];
    }
  }
  
  for( int t = 1; t < Transitions; ++t ){
    if(!missingGenotypes[t-1]) {
      Sum = 0.0;
      //scale previous alpha to sum to 1
      for( int j = 0; j <  DStates; ++j ) {
	Sum += alpha[(t-1)*DStates +j];
      }
      scaleFactor = 1.0 / Sum;
      for( int j = 0; j <  DStates; ++j ) {
	alpha[(t-1)*DStates +j] *= scaleFactor;
      }
      sumfactor += log(Sum);
    }
    
    RecursionProbs(p[t], f + 2*t, StateArrivalProbs + t*K*2, alpha + (t-1)*DStates, alpha + t*DStates);
    
    for(int j = 0; j < DStates; ++j){
      if(!missingGenotypes[t]) {
	alpha[t*DStates +j] *=  Lambda[t*DStates + j]; //*lam++; 
	//++lam;
      } //else ++lam;
    }
  }
  alphaIsBad = false;
}

void HMM::UpdateBackwardProbsDiploid()
{
  if(!Lambda || !theta || !f)throw string("Error: Call to HMM when inputs are not set!");
  if(!betaAllocated) { // allocate beta array if not already done
    beta =  new double[Transitions*K*K];
    LambdaBeta = new double[K*K];
    betaAllocated = true;
  }
  vector<double> rec(DStates);
  double scaleFactor, Sum;
  
  for(int j = 0; j < DStates; ++j){
    //set beta(T) = 1
    beta[(Transitions - 1)*DStates + j] = 1.0;
    rec[j] = 1.0 / ThetaThetaPrime[j];
  }
  
  for( int t = Transitions-2; t >= 0; --t ) {
    double f2[2] = {f[2*t + 2], f[2*t + 3]};
    Sum = 0.0;
    for(int j = 0; j < DStates; ++j){
      LambdaBeta[j] = Lambda[(t+1)*K*K + j] * beta[(t+1)*DStates + j] * ThetaThetaPrime[j];
      Sum += LambdaBeta[j];
    }
    //scale LambdaBeta to sum to 1
    scaleFactor = 1.0 / Sum;
    for( int j = 0; j <  DStates; ++j ) {
      LambdaBeta[j] *= scaleFactor;
    }
    
    RecursionProbs(p[t+1], f2, StateArrivalProbs+(t+1)*K*2, LambdaBeta, beta+ t*DStates);
    for(int j = 0; j < DStates; ++j){ // vectorization successful
      beta[t*DStates + j] *= rec[j];
    }
  }
  betaIsBad = false;
}

/*
  Updates forward probs, haploid case only
  Here Admixture is a column matrix and the last dimensions of f and lambda are 1.
*/
void HMM::UpdateForwardProbsHaploid(){
  if(!Lambda || !theta || !f)throw string("Error: Call to HMM when inputs are not set!");
  sumfactor = 0.0;
  double Sum = 0.0;
  double scaleFactor = 0.0;
  for(int j = 0; j < K; ++j){
    alpha[j] = theta[j] * Lambda[j];
  }
  
  for( int t = 1; t < Transitions; t++ ) {
    if(!missingGenotypes[t-1]) {
      Sum = 0.0;
      //scale previous alpha to sum to 1
      for(int j = 0; j < K; ++j){
	Sum += alpha[(t-1)*K + j];
      }
      scaleFactor = 1.0 / Sum;
      for( int j = 0; j <  K; ++j ) {
	alpha[(t-1)*K +j] *= scaleFactor;
      }
      sumfactor += log(Sum);
    }

    for(int j = 0; j < K; ++j){
      alpha[t*K + j] = f[2*t]*alpha[(t-1)*K +j] + (1.0 - f[2*t]) * theta[j];
      alpha[t*K + j] *= Lambda[t*K + j];
    }
  }
  alphaIsBad =  false;
}

void HMM::UpdateBackwardProbsHaploid(){
  if(!Lambda || !theta || !f)throw string("Error: Call to HMM when inputs are not set!");
  if(!betaAllocated) { // allocate diploid-sized beta array if not already done
    beta =  new double[Transitions*K*K];
    LambdaBeta = new double[K*K];
    betaAllocated = true;
  }
  double Sum = 0.0;
  for(int j = 0; j < K; ++j){
    beta[(Transitions-1)*K + j] = 1.0;
  }
  
  for( int t = Transitions-2; t >=0; t-- ){
    Sum = 0.0;
    for(int j = 0; j < K; ++j){
      LambdaBeta[j] = Lambda[(t+1)*K + j]*beta[(t+1)*K + j];
      Sum += theta[j]*LambdaBeta[j];
    }
    for(int j = 0; j < K; ++j){
      beta[t*K + j] = f[2*(t+1)]*LambdaBeta[j] + (1.0 - f[2*(t+1)+1])*Sum;
    }
  }
  betaIsBad = false;
}

// argument oldProbs is square array of size K * K
// for forward recursions, pass alpha_t and multiply newProbs by emission probs lambda_t 
// for backward recursions, pass array of products lambda_t+1[jj] * beta_t+1[jj] 
void HMM::RecursionProbs(const double ff, const double f2[2], const double* const stateArrivalProbs, 
			 const double* const oldProbs, double *newProbs) {
  if(K==2) RecursionProbs2(ff, f2, stateArrivalProbs, oldProbs, newProbs);
  else {
    int j0K = 0;
    for( int j0 = 0; j0 <  K; ++j0 ) {
      rowProb[j0] = 0.0;
      colProb[j0] = 0.0;
      for( int j1 =0; j1 < K; ++j1 ) {
	rowProb[j0] += oldProbs[j0*K + j1];
	colProb[j0] += oldProbs[j1*K + j0];
      }
    }
    // calculate expectations of indicator variables for each ancestry state on each gamete
    for( int j0 = 0; j0 <  K; ++j0 ) {
      Expectation0[j0] = f2[0]*rowProb[j0] + stateArrivalProbs[j0*2];
      Expectation1[j0] = f2[1]*colProb[j0] + stateArrivalProbs[j0*2 + 1];
    }

    // calculate expectation of product as covariance plus product of expectations
    for(int j0 = 0; j0 < K; ++j0) {
      for(int j1 = 0; j1 < K; ++j1) {
        j0K = j0 * K;
	cov[j0K + j1] = ff * ( oldProbs[j0K + j1] - rowProb[j0] * colProb[j1] );
	newProbs[j0K + j1] = cov[j0K + j1] + Expectation0[j0] * Expectation1[j1];
	// newProbs[1] is prob(paternal=1, maternal=0)
      }
    }
  }//end else
}
 
void HMM::RecursionProbs2(const double ff, const double f2[2], const double* const stateArrivalProbs, 
			  const double* const oldProbs, double *newProbs) {
  // version for 2 subpopulations
  double row0Prob;
  double col0Prob;
  double Exp0;
  double Exp1;
  // sum row 0 and col 0  
  row0Prob = ( oldProbs[0] + oldProbs[2] );
  col0Prob = ( oldProbs[0] + oldProbs[1] );
  // calculate expectations of indicator variables for ancestry=0 on each gamete
  Exp0 = f2[0]*row0Prob + stateArrivalProbs[0]; // paternal gamete
  Exp1 = f2[1]*col0Prob + stateArrivalProbs[1]; // maternal gamete
  // calculate covariance of indicator variables as ff * deviation from product of row and col probs
  newProbs[0] = Exp0 * Exp1 + ff * ( oldProbs[0] - row0Prob * col0Prob );; 
  newProbs[1] = Exp0 - newProbs[0]; //prob paternal ancestry=1, maternal=0 
  newProbs[2] = Exp1 - newProbs[0]; //prob paternal ancestry=0, maternal=1 
  newProbs[3] = 1 - Exp0 - Exp1 + newProbs[0];
}




