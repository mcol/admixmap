#include "Chromosome.h"
#include "Individual.h"
#include <iostream>
#include "functions.h"

#define PR(x) cerr << #x << " = " << x << endl;

using namespace std;

Chromosome::Chromosome(int size, int start, int inpopulations) : Genome(size)
{
  _startLoci = start;
  populations = inpopulations;
  D = populations * populations;
  L = size;

  SampleStates.SetDimensions( size, populations, true );

  Lambda = new double**[size];

  for(int i = 0; i < size ;++i){
    Lambda[i] = alloc2D_d(populations, populations);
  }

  CodedStates = new int[size];
  for(int j = 0; j < 2; ++j) f[j] = new double[size];

  StateArrivalProbs = new double**[L];
  for(int t = 1; t < size; t++ ){        
    StateArrivalProbs[t] = alloc2D_d(populations,2);
  }

}

void Chromosome::ResetStuffForX()
{
  D = populations;
  SampleStates.SetDimensions( L, populations, false );
}

void Chromosome::SetLabel( string label )
{
   _Label = label;
}

string Chromosome::GetLabel( int )
{
   return _Label;
}

Chromosome::~Chromosome()
{
  //TODO: delete these properly
  delete CodedStates;
  delete[] Lambda;
  delete[] f[0];
  delete[] f[1];
  for(unsigned  int t = 1; t < L; t++ ){
    for(int j = 0; j < populations; ++j){
      delete[] StateArrivalProbs[t][j];
    }
    delete[] StateArrivalProbs[t];
  }
  delete[] StateArrivalProbs;
}


// Returns the number of the num'th compositelocus on this chromosome
// eg if chromosome 2 consists of compositeloci 5,6,7 and 8,
// GetLocus(i) returns 5 + i
int
Chromosome::GetLocus(int num){
  return _startLoci + num;
}

//returns number of composite loci in the chromosome
unsigned int Chromosome::GetSize(){
  return L;
}

//Initialises f for global rho. Necessary since individual-level parameters updated before global rho (in Latent)
void Chromosome::InitialiseLociCorr(const double rho){
  for(unsigned int j = 1; j < L; j++ )
    f[0][j] = f[1][j] = ( -GetDistance( j ) * rho > -700) ? exp( -GetDistance( j ) * rho ) : 0.0;
}

//sets f for global rho, called after Latent is updated
void Chromosome::SetLociCorr(const double rho){
    for(unsigned int jj = 1; jj < L; jj++ ){
      f[0][jj] = f[1][jj] = exp( -GetDistance( jj ) * rho );
    }
}

void Chromosome::UpdateParameters(Individual* ind, AlleleFreqs *A, Matrix_d& Admixture, AdmixOptions* options,  
				  std::vector< double > _rho, bool fixedallelefreqs, bool diploid ){
  // f0 and f1 are arrays of scalars of the form exp - rho*x, where x is distance between loci
  // required to calculate transition matrices 
  if( options->getRhoIndicator() ){

    for( unsigned int jj = 1; jj < L; jj++ ){
      f[0][jj] = exp( -GetDistance( jj ) * _rho[0] );
      if( options->isRandomMatingModel() ){
	  f[1][jj] = exp( -GetDistance( jj ) * _rho[1] );
      }
      else
	f[1][jj] = f[0][jj];
    }
  }
  //global rho case already dealt with

 int locus = GetLocus(0);
  bool test = (options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry());
  if(diploid){
    //construct Lambda
    for(unsigned int j = 0; j < L; j++ ){
      if( !(ind->IsMissing(locus)) ){
	A->GetGenotypeProbs(Lambda[j], locus, ind->getGenotype(locus), ind->getPossibleHapPairs(locus), true, fixedallelefreqs );
      }
      else{
	for( int k1 = 0; k1 < populations; k1++ )for(int k2 =0; k2<populations; ++k2) Lambda[j][k1][k2] = 1.0;
      }
      locus++;
    }
    //construct StateArrivalProbs
    for(unsigned int t = 1; t < L; t++ ){        
      for(int j = 0; j < populations; ++j){
	StateArrivalProbs[t][j][0] = (1.0 - f[0][t]) * Admixture(j,0);
	if( options->isRandomMatingModel())
	  StateArrivalProbs[t][j][1] = (1.0 - f[1][t]) * Admixture(j,1);
	else StateArrivalProbs[t][j][1] = (1.0 - f[1][t]) * Admixture(j,0);
      }
    }
 
    //Update Forward/Backward Probs in HMM
    SampleStates.UpdateForwardProbsDiploid(StateArrivalProbs,f, Admixture, Lambda, options->isRandomMatingModel());
    if(test){
      double **ThetaThetaPrime = alloc2D_d(populations, populations);
      for(int j0 = 0; j0 < populations; ++j0)for(int j1 = 0; j1 < populations; ++j1)
	ThetaThetaPrime[j0][j1] = Admixture(j0,0)*Admixture(j1, options->isRandomMatingModel());
      SampleStates.UpdateBackwardProbsDiploid(StateArrivalProbs,f, ThetaThetaPrime, Lambda);
      free_matrix(ThetaThetaPrime, populations);
    }

  }

  else{//haploid
    for(unsigned int j = 0; j < L; j++ ){
      if( !(ind->IsMissing(locus)) ){
	A->GetGenotypeProbs(Lambda[j], locus, ind->getGenotype(locus), ind->getPossibleHapPairs(locus), false, fixedallelefreqs );
      }
      else{
	for( int k1 = 0; k1 < populations; k1++ )for(int k2 = 0; k2 < populations; ++k2) Lambda[j][k1][k2] = 1.0;
      }
    }
    //Update Forward/Backward Probs in HMM
    SampleStates.UpdateProbsHaploid(f, Admixture, Lambda, test);
  }

}

//void Chromosome::SampleLocusAncestry(Matrix_i *OrderedStates, Matrix_d &Admixture, double *f[], bool isdiploid){
void Chromosome::SampleLocusAncestry(Matrix_i *OrderedStates, Matrix_d &Admixture, bool isdiploid){

  SampleStates.Sample(CodedStates, Admixture, f, StateArrivalProbs, isdiploid);
  for(unsigned int j = 0; j < L; j++ ){
    if(isdiploid){
      //     OrderedStates( 0, j ) = 0;
      (*OrderedStates)( 0, j ) = (int)(CodedStates[j] / populations);
      (*OrderedStates)( 1, j ) = (CodedStates[j] % populations);
    }
    else //haploid
      {
	(*OrderedStates)(0, j ) = CodedStates[j];
      }
  }
}

void Chromosome::getAncestryProbs( int j, double AncestryProbs[][3] ){
  //
  //sets conditional probabilities of ancestry at locus j
  //One row per population, Cols 0,1,2 are probs that 0,1,2 of the 2 gametes have ancestry from that population
  //i.e. (i,2) = p_{ii}
  //     (i,1) = \sum_j{p_{ij}} +   \sum_j{p_{ji}} - 2.0*p_{ii}
  //     (i,0) = 1.0 - (i,1) - (i,2)
  //where p's are probs in StateProbs from HMM

  double *StateProbs;
  StateProbs = new double[D];//possibly should keep this at class scope
  
  SampleStates.GetStateProbs(StateProbs, j);
 
  for( int k1 = 0; k1 < populations; k1++ ){
    AncestryProbs[k1][2] = StateProbs[ ( populations + 1 ) * k1 ];
    AncestryProbs[k1][1] = 0.0;
    for( int k2 = 0 ; k2 < populations; k2++ )
      AncestryProbs[k1][1] += StateProbs[k1*populations +k2] + StateProbs[k2*populations +k1];
    AncestryProbs[k1][1] -= 2.0*AncestryProbs[k1][2];
    AncestryProbs[k1][0] = 1.0 - AncestryProbs[k1][1] - AncestryProbs[k1][2];
  }
  delete[] StateProbs;
}

//accessor for HMM Likelihood
double Chromosome::getLogLikelihood()
{
  return SampleStates.getLikelihood();
}

//samples jump indicators xi for this chromosome, 
//updates sumxi and Sumrho0
void Chromosome::SampleJumpIndicators(const Matrix_i &LocusAncestry, 
				      const unsigned int gametes, std::vector< std::vector<bool> > *xi, int *sumxi, 
				      double *Sumrho0){

  int locus;
  double Prob;
  for( unsigned int jj = 1; jj < L; jj++ ){
    locus = GetLocus(jj);    
    for( unsigned int g = 0; g < gametes; g++ ){
      if( LocusAncestry(g,jj-1) == LocusAncestry(g,jj) ){

	Prob = StateArrivalProbs[jj][ LocusAncestry(g,jj)][g] / (StateArrivalProbs[jj][ LocusAncestry(g,jj)][g] + f[g][jj] );
	if( Prob > myrand() ){
	  (*xi)[g][locus] = true;
	  sumxi[locus]++;
	} else {
	  (*xi)[g][locus] = false;
	  *Sumrho0 += GetDistance( jj );
	}
      } else {
	(*xi)[g][locus] = true;
	sumxi[locus]++;
      }
    }
  }
}
