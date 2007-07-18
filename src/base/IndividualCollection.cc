/** 
 *   IndividualCollection.cc
 *   Class to hold an array of Individuals
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "IndividualCollection.h"
#include "bclib/Regression.h"

using namespace std;

// **** CONSTRUCTORS  ****
// IndividualCollection::IndividualCollection() {
//   SetNullValues();
// }
void IndividualCollection::SetNullValues(){
  OutcomeType = 0;
  NumOutcomes = 0;
  NumCovariates = 0;
  NumberOfInputCovariates = 0;
  ReportedAncestry = 0;
  SumDeviance = SumDevianceSq = 0.0;
  //SumEnergy = 0; SumEnergySq = 0;
  _child = 0;
  SumLogLikelihood = 0.0;
  ReportedAncestry = 0;
  size = 0;
  NumDiploidIndividuals = 0;
}

IndividualCollection::IndividualCollection(unsigned numIndividuals, unsigned numPopulations, unsigned numCompositeLoci) :
  NumInd(numIndividuals), Populations(numPopulations), NumCompLoci(numCompositeLoci){
  SetNullValues();
  //Populations = options->getPopulations();
  //NumInd = Data->getNumberOfIndividuals();
  size = NumInd;
  //NumCompLoci = Loci->GetNumberOfCompositeLoci();
}

int IndividualCollection::getNumDiploidIndividuals()const{
  return NumDiploidIndividuals;
}

// ************** DESTRUCTOR **************
IndividualCollection::~IndividualCollection() {
  //cout << "\n Deleting individual objects\n" << flush;
  //  AdmixedIndividual::DeleteStaticMembers();
  for(unsigned int i = 0; i < size; i++){
    delete _child[i];
  }

//  delete[] _child;

  delete[] OutcomeType;
  delete[] ReportedAncestry;
}
void IndividualCollection::DeleteGenotypes(bool setmissing=false){
  for (unsigned int i = 0; i < size; i++) {
    if(setmissing)_child[i]->SetMissingGenotypes();
    _child[i]->DeleteGenotypes();
  }
}

// ************** INITIALISATION AND LOADING OF DATA **************

void IndividualCollection::LoadData(const Options* const options, const InputData* const data_, bool admixtureAsCovariate){

  if ( options->getNumberOfOutcomes()>0){
    delete[] OutcomeType;
    OutcomeType = new DataType[ options->getNumberOfOutcomes() ];
    data_->getOutcomeTypes(OutcomeType);

    if(strlen( options->getOutcomeVarFilename() ) != 0)LoadOutcomeVar(data_);
    LoadCovariates(data_, options, admixtureAsCovariate);
  }
}

void IndividualCollection::LoadCovariates(const InputData* const data_, const Options* const options, bool admixtureAsCovariate){
  if ( strlen( options->getCovariatesFilename() ) > 0 ){
    bclib::DataMatrix& CovData = (bclib::DataMatrix&)data_->getCovariatesMatrix();
    NumberOfInputCovariates = CovData.nCols();
    unsigned NumInds = CovData.nRows();//should already have been checked to be the same as in outcomevarfile

    if( admixtureAsCovariate){
      Covariates.setDimensions(NumInds, CovData.nCols() + options->getPopulations());
      for(unsigned i = 0; i < NumInds; ++i)for(int j = 0; j < options->getPopulations()-1; ++j)
	Covariates.set(i, j + CovData.nCols() + 1,  1.0 / (double)options->getPopulations() );
    } else
      Covariates.setDimensions(NumInds, CovData.nCols()+1);//+1 for intercept
    for(unsigned i = 0; i < NumInds; ++i)for(unsigned j = 0; j < CovData.nCols(); ++j)
      {
	Covariates.set(i,0, 1.0); //set to 1 for intercept
	Covariates.set(i,j+1, CovData.get(i, j) );
	Covariates.isMissing(i,j+1, CovData.isMissing(i,j));
      }
    if ( Covariates.hasMissing() ) Covariates.SetMissingValuesToColumnMeans();

    //centre covariates about their means
    // (this should already be done by user)    
    //vector<double> mean = Covariates.columnMeans();
    //for(unsigned int i = 0; i < NumInds; i++ )
    // for( int j = 1; j <= NumberOfInputCovariates; j++ )
    //Covariates( i, j ) -= mean[j];
  }
  else {//no covariatesfile
    unsigned NumInds = Outcome.nRows();
    if(NumInds <= 0 )NumInds = NumInd;
    if( admixtureAsCovariate ){
      Covariates.setDimensions(NumInds, options->getPopulations());
      for(unsigned i = 0; i < NumInds; ++i)for(int j = 1; j < options->getPopulations(); ++j)
	Covariates.set(i, j, 1.0 / (double)options->getPopulations() );
    } else
      Covariates.setDimensions(NumInds, 1);//just an intercept
    for(unsigned i = 0; i < NumInds; ++i)
      Covariates.set(i,0, 1.0 );
  }
  if( !admixtureAsCovariate )
    NumCovariates = NumberOfInputCovariates + 1;
  else
    NumCovariates = NumberOfInputCovariates + options->getPopulations();
}

void IndividualCollection::LoadOutcomeVar(const InputData* const data_){
  Outcome = data_->getOutcomeVarMatrix();
  NumOutcomes = Outcome.nCols();
 
}

void IndividualCollection::HMMIsBad(bool b){
  for(unsigned i = 0; i < size; i++)
    _child[i]->HMMIsBad(b);
}

/**
   Samples Haplotype pairs and upates allele/haplotype counts if requested
*/
void IndividualCollection::SampleHapPairs(const Options& options, AlleleFreqs *A, const Genome* const Loci,
					  bool skipMissingGenotypes, bool anneal, bool UpdateAlleleCounts){
  unsigned nchr = Loci->GetNumberOfChromosomes();
  unsigned locus = 0;
  // if annealthermo, no need to sample hap pair: just update allele counts if diallelic
  bool annealthermo = anneal && options.getThermoIndicator() && !options.getTestOneIndivIndicator();
  for(unsigned j = 0; j < nchr; ++j){
    for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){
      
      for(unsigned int i = 0; i < size; i++ ){
	//Sample Haplotype Pair
	//also updates allele counts unless using hamiltonian sampler at locus with > 2 alleles 
	_child[i]->SampleHapPair(j, jj, locus, A, skipMissingGenotypes, annealthermo, UpdateAlleleCounts);
      }
      locus++;
    }
  }
}

void IndividualCollection::AccumulateAlleleCounts(const Options& options, AlleleFreqs *A, const Genome* const Loci,
                                                  bool anneal){
  unsigned nchr = Loci->GetNumberOfChromosomes();
  unsigned locus = 0;
  // if annealthermo, no need to sample hap pair: just update allele counts if diallelic
  bool annealthermo = anneal && options.getThermoIndicator() && !options.getTestOneIndivIndicator();
  for(unsigned j = 0; j < nchr; ++j){
    for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){

      for(unsigned int i = 0; i < size; i++ ){
        _child[i]->UpdateAlleleCounts(j, jj, locus, A, annealthermo);
      }
      locus++;
    }
  }
}

// ************** ACCESSORS **************
int IndividualCollection::getSize()const {
  return size;
}

const vector<double> IndividualCollection::getOutcome(int j)const{
  vector<double> col;
  if(j < (int)Outcome.nCols())
    col = Outcome.getCol(j);
  return col;
}

double IndividualCollection::getOutcome(int j, int ind)const{
    return Outcome.get(ind, j);
}
bool IndividualCollection::isMissingOutcome(int j, int i)const{
    return Outcome.isMissing(i, j);
}
int IndividualCollection::getNumberOfOutcomeVars()const{
  return NumOutcomes;
}

DataType IndividualCollection::getOutcomeType(int i)const{
  return OutcomeType[i];
}

Individual* IndividualCollection::getIndividual(int num)const
{
  if (num < (int)size){
    return _child[num];
  } else {
    return 0;
  }
}

int IndividualCollection::GetNumberOfInputCovariates()const{
  return NumberOfInputCovariates;
}
int IndividualCollection::GetNumCovariates() const{
  return NumCovariates;
}

const bclib::DataMatrix& IndividualCollection::getCovariatesMatrix()const{
  return Covariates;
}
const bclib::DataMatrix& IndividualCollection::getOutcomeMatrix()const{
  return Outcome;
}

/**
 * returns a count of the copies of allele a at a comp locus.
 * Only works for diallelic loci.
 * used in strat test to determine loci with given percentage of observed genotypes 
 */
unsigned IndividualCollection::GetSNPAlleleCounts(unsigned locus, int allele)const{
  int AlleleCounts = 0;
  for(unsigned i = 0; i < size; i++){
    if(!_child[i]->GenotypeIsMissing(locus)){
      const int* haps = _child[i]->getSampledHapPair(locus);
      if(haps[0] == allele-1){// -1 because alleles count from 1 and haps count from 0
	AlleleCounts++;
      }
      if(haps[1] == allele-1){
	AlleleCounts++;
      }
    }
  }
  return AlleleCounts;
}

const vector<int> IndividualCollection::getAlleleCounts(unsigned locus, int pop, unsigned NumStates)const{
  int ancestry[2];
  vector<int> counts(NumStates);
  fill(counts.begin(), counts.end(), 0);
  for(unsigned i = 0; i < size; i++)
    if( !_child[i]->GenotypeIsMissing(locus)){
      _child[i]->GetLocusAncestry(locus, ancestry);
      const int* happair = _child[i]->getSampledHapPair(locus);
      // happair[1] ==  -1 if haploid
      if(ancestry[0] == pop && (happair[0] >= 0) )++counts[happair[0]];
      if(ancestry[1] == pop&& (happair[1] >= 0) )++counts[happair[1]];
    }
  return counts;
}
///count number of missing genotypes at locus
int IndividualCollection::getNumberOfMissingGenotypes(unsigned locus)const{
  int count = 0;
  for(unsigned i = 0; i < size; i++){
    if(_child[i]->GenotypeIsMissing(locus)){
      ++count;
    }
  }
  return count;
}

// ************** OUTPUT **************



double IndividualCollection::getLogLikelihood(const Options& options, bool forceupdate){
  double LogLikelihood = 0.0;
  for(unsigned i = 0; i < size; i++) {
    LogLikelihood += _child[i]->getLogLikelihood(options, forceupdate, true); // store result if updated
    _child[i]->HMMIsBad(false); // HMM probs overwritten by next indiv, but stored loglikelihood still ok
  }
  return LogLikelihood;
}

double IndividualCollection::getEnergy(const Options& options, const vector<bclib::Regression*> &R, 
				       const bool & annealed) {
  // energy is minus the unnannealed log-likelihood summed over all individuals under study from both HMM and regression 
  // called every iteration after burnin, after update of genotype probs and before annealing
  // accumulates sums of deviance and squared deviance
  double LogLikHMM = 0.0;
  double LogLikRegression = 0.0;
  double Energy = 0.0;
  // assume that HMM probs and stored loglikelihoods are bad, as this function is called after update of allele freqs  
  for(unsigned i = 0; i < size; i++) {
    LogLikHMM += _child[i]->getLogLikelihood(options, false, !annealed); // store result if not an annealed run
    // don't have to force an HMM update here - on even-numbered iterations with globalrho, stored loglikelihood is still valid
    
    if(annealed)  _child[i]->HMMIsBad(true); // HMM probs bad, stored loglikelihood bad
    else _child[i]->HMMIsBad(false); 
  }
  // get regression log-likelihood 
  for(unsigned c = 0; c < R.size(); ++c) LogLikRegression += R[c]->getLogLikelihood(getOutcome(c));
  Energy = -(LogLikHMM + LogLikRegression);
  return Energy;
} 


