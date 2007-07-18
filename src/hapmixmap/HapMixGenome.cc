#include "HapMixGenome.hh"
#include "HapMixHMM.hh"
#include "bclib/misc.h"//for myexp

///allocate a HapMixChromosome. overrides base class function which allocates a (base) Chromosome
void HapMixGenome::CreateChromosome(unsigned i, unsigned size, bool isX, unsigned cstart, int NumHiddenStates ){
  vHapMixChr.push_back(new HapMixChromosome(i, size, cstart, NumHiddenStates, isX));
  C[i] = vHapMixChr[i];
}

///sets locus-specific f 
void HapMixGenome::SetLocusCorrelation(const std::vector<double>& lambda){
  //check lambda has length equal to number of intervals
  //  if(lambda.size()<NumberOfCompositeLoci-NumberOfChromosomes)throw string("Bad arguments passed to Chromosome::SetLocusCorr");

  vector<double>::const_iterator lambda_iter = lambda.begin();
  for( unsigned int j = 0 ; j < NumberOfChromosomes; j++ ) {
    vHapMixChr[j]->SetLocusCorrelation(lambda_iter);
    lambda_iter += C[j]->GetSize()-1;
  }
}

void HapMixGenome::SetHMMDimensions(int NumHiddenStates, bool diploid){
  for(vector<HapMixChromosome*>::iterator c = vHapMixChr.begin(); c != vHapMixChr.end(); ++c ){
    (*c)->SetHMMDimensions(NumHiddenStates, diploid); 
  }
}

HapMixChromosome* HapMixGenome::getHapMixChromosome(unsigned j){
  return vHapMixChr[j];
}

///constructor for HapMixChromosome. The same as base class but creates HapMixHMM object
HapMixChromosome::HapMixChromosome(int n, int size, int start, int inNumHiddenStates, bool isx):
  Chromosome(n, size, start, inNumHiddenStates, isx){

  HMM = new HapMixHMM( size, inNumHiddenStates, f);
}

///Sets locus ancestry correlations f for locus-specific arrival rate lambda (= rho*distance).
void HapMixChromosome::SetLocusCorrelation(const std::vector<double>::const_iterator lambda_iter){
  for(unsigned int j = 1; j < NumberOfCompositeLoci; j++ ){
    double lambda = *(lambda_iter + j -1);
    if(isX)lambda *= 0.5;
    f[2*j] = f[2*j + 1] = bclib::eh_exp( - lambda );
  }

}

void HapMixChromosome::SetHMMDimensions(int NumHiddenStates, bool diploid){
  ((HapMixHMM*)HMM)->SetDimensions(NumberOfCompositeLoci, NumHiddenStates, f, diploid);
}
