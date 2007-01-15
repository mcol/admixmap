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
#include "regression/Regression.h"
#include "Comms.h"
#ifdef PARALLEL
#include <mpe.h>//for MPI event logging
#endif
using namespace std;

// **** CONSTRUCTORS  ****
IndividualCollection::IndividualCollection() {
    SetNullValues();
}
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
}

IndividualCollection::IndividualCollection(const AdmixOptions* const options, const InputData* const Data, Genome* Loci) {
  SetNullValues();
  Populations = options->getPopulations();
  NumInd = Data->getNumberOfIndividuals();
  size = NumInd;
  NumCompLoci = Loci->GetNumberOfCompositeLoci();
  worker_rank = 0;
  NumWorkers = 1;
#ifdef PARALLEL
    int global_rank = MPI::COMM_WORLD.Get_rank();
    //create communicator for workers and find size of and rank within this group
    workers = MPI::COMM_WORLD.Split( Comms::isWorker(), global_rank);
    NumWorkers = workers.Get_size();
    if(global_rank >1)
      worker_rank = workers.Get_rank();
    else worker_rank = size;//so that non-workers will not loop over Individuals
#endif
}

int IndividualCollection::getNumDiploidIndividuals(){
  int numdiploid = 0;
    for (unsigned int i = worker_rank; i < size; i += NumWorkers) {
      if(!_child[i]->isHaploidIndividual())++numdiploid;
    }
#ifdef PARALLEL
    Comms::Reduce(&numdiploid);
#endif
    return numdiploid;
}

// ************** DESTRUCTOR **************
IndividualCollection::~IndividualCollection() {
  //if(worker_rank==0)
  //cout << "\n Deleting individual objects\n" << flush;
  //  AdmixedIndividual::DeleteStaticMembers();
  for(unsigned int i = worker_rank; i < size; i+=NumWorkers){
    delete _child[i];
  }

//  delete[] _child;

  delete[] OutcomeType;
  delete[] ReportedAncestry;
}
void IndividualCollection::DeleteGenotypes(bool setmissing=false){
  for (unsigned int i = worker_rank; i < size; i += NumWorkers) {
    if(setmissing)_child[i]->SetMissingGenotypes();
    _child[i]->DeleteGenotypes();
  }
}

// ************** INITIALISATION AND LOADING OF DATA **************

void IndividualCollection::LoadData(const AdmixOptions* const options, const InputData* const data_, bool admixtureAsCovariate){

  if ( options->getNumberOfOutcomes()>0){
    delete[] OutcomeType;
    OutcomeType = new DataType[ options->getNumberOfOutcomes() ];
    data_->getOutcomeTypes(OutcomeType);

    if(strlen( options->getOutcomeVarFilename() ) != 0)LoadOutcomeVar(data_);
    LoadCovariates(data_, options, admixtureAsCovariate);
  }
  if ( strlen( options->getReportedAncestryFilename() ) != 0 ){
    LoadRepAncestry(data_);
  }
}

void IndividualCollection::LoadCovariates(const InputData* const data_, const AdmixOptions* const options, bool admixtureAsCovariate){
//   NumCovariates = 0;
//   DataMatrix& CovData = (DataMatrix&)data_->getCovariatesMatrix();
//   if(!(options->getNumberOfOutcomes()==1 && OutcomeType[0]==CoxData))
//     ++NumCovariates;//for intercept. No intercept for Cox regression
//   //TODO: what to do if Cox regression and another regression type
//   if( !options->getTestForAdmixtureAssociation() && !options->getHapMixModelIndicator() ){
//     NumCovariates += options->getPopulations()-1;//admixture term for each population, except the first
//   }
//   NumberOfInputCovariates = CovData.nCols();
//   NumCovariates += NumberOfInputCovariates;//add number of covariates in file
//   Covariates.setDimensions(size, NumCovariates);


  if ( strlen( options->getCovariatesFilename() ) > 0 ){
    DataMatrix& CovData = (DataMatrix&)data_->getCovariatesMatrix();
    NumberOfInputCovariates = CovData.nCols();
    unsigned NumInds = CovData.nRows()-1;//should already have been checked to be the same as in outcomevarfile

    if( admixtureAsCovariate){
      Covariates.setDimensions(NumInds, CovData.nCols() + options->getPopulations());
      for(unsigned i = 0; i < NumInds; ++i)for(int j = 0; j < options->getPopulations()-1; ++j)
	Covariates.set(i, j + CovData.nCols() + 1,  1.0 / (double)options->getPopulations() );
    } else
      Covariates.setDimensions(NumInds, CovData.nCols()+1);//+1 for intercept
    for(unsigned i = 0; i < NumInds; ++i)for(unsigned j = 0; j < CovData.nCols(); ++j)
      {
	Covariates.set(i,0, 1.0); //set to 1 for intercept
	Covariates.set(i,j+1, CovData.get(i+1, j) );
	Covariates.isMissing(i,j+1, CovData.isMissing(i+1,j));
      }
    if ( Covariates.hasMissing() ) Covariates.SetMissingValuesToColumnMeans();
    
    //vector<double> mean(NumberOfInputCovariates);
    
    //centre covariates about their means
    // (this should already be done by user)
    //    for( int j = 0; j < NumberOfInputCovariates; j++ ){
    //    int count = 0;
    //  mean[j] = 0.0;
    //for(unsigned int i = 0; i < NumInds; i++ )
    //if(!Covariates.IsMissingValue(i,j+1)){
    // mean[j] += Covariates(i,j+1);
    //++count;
    //}
    //mean[j] /= (double)count;
    //}
    
    //for(unsigned int i = 0; i < NumInds; i++ )
    // for( int j = 0; j < NumberOfInputCovariates; j++ )
    //Covariates( i, j+1 ) -= mean[j];
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
  //if(size != Outcome.nRows() && size!= NumInd)throw string("ERROR in outcomevarfile: wrong number of rows\n");
  NumOutcomes = Outcome.nCols();
 
}

void IndividualCollection::LoadRepAncestry(const InputData* const data_){
  ReportedAncestry = new DataMatrix[NumInd];
  DataMatrix& temporary = (DataMatrix&)data_->getReportedAncestryMatrix();
  for( unsigned i = 0; i < temporary.nRows() / 2; i++ )
    ReportedAncestry[i] = temporary.SubMatrix( 2*i, 2*i + 1, 0, temporary.nCols() - 1 );
 
}

void IndividualCollection::HMMIsBad(bool b){
  for(unsigned i = worker_rank; i < size; i+= NumWorkers)
    _child[i]->HMMIsBad(b);
}

void IndividualCollection::setGenotypeProbs(const Genome* const Loci, const AlleleFreqs* const
#ifdef PARALLEL
					    A
#endif
					    ){
  unsigned nchr = Loci->GetNumberOfChromosomes();
  unsigned locus = 0;
  for(unsigned j = 0; j < nchr; ++j){
    for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){
#ifdef PARALLEL
      //get pointer to allele probs for this locus
      const double* AlleleProbs = A->GetAlleleFreqs(locus);//need to get alleleprobs from A as workers have no CompositeLocus objects
      for(unsigned int i = worker_rank; i < size; i+= NumWorkers ) {
	_child[i]->SetGenotypeProbs(j, jj, locus, AlleleProbs);
      }
#else
      for(unsigned int i = worker_rank; i < size; i+= NumWorkers ) {
	_child[i]->SetGenotypeProbs(j, jj, locus, false);
      }
#endif
      locus++;
    }
  }
}  

void IndividualCollection::annealGenotypeProbs(unsigned nchr, const double coolness, const double* ){
  for(unsigned j = 0; j < nchr; ++j){
      for(unsigned int i = worker_rank; i < size; i+= NumWorkers) {
	_child[i]->AnnealGenotypeProbs(j, coolness);
      }
  }
}

/**
   Samples Haplotype pairs and upates allele/haplotype counts
*/
void IndividualCollection::SampleHapPairs(const AdmixOptions* const options, AlleleFreqs *A, const Genome* const Loci,
					  bool skipMissingGenotypes, bool anneal){
  unsigned nchr = Loci->GetNumberOfChromosomes();
  unsigned locus = 0;
  // if annealthermo, no need to sample hap pair: just update allele counts if diallelic
  bool annealthermo = anneal && options->getThermoIndicator() && !options->getTestOneIndivIndicator();
  for(unsigned j = 0; j < nchr; ++j){
    for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){
#ifdef PARALLEL
      const double* AlleleProbs = A->GetAlleleFreqs(locus);
#endif
      
      for(unsigned int i = worker_rank; i < size; i+=NumWorkers ){
	//Sample Haplotype Pair
	//also updates allele counts unless using hamiltonian sampler at locus with > 2 alleles 
#ifdef PARALLEL
	_child[i]->SampleHapPair(j, jj, locus, A, skipMissingGenotypes, annealthermo, AlleleProbs);
#else
	_child[i]->SampleHapPair(j, jj, locus, A, skipMissingGenotypes, annealthermo);//also updates allele counts
#endif
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

const DataMatrix& IndividualCollection::getCovariatesMatrix()const{
  return Covariates;
}
const DataMatrix& IndividualCollection::getOutcomeMatrix()const{
  return Outcome;
}

/**
 * returns a count of the copies of allele a at a comp locus.
 * Only works for diallelic loci.
 * used in strat test to determine loci with given percentage of observed genotypes 
 */
unsigned IndividualCollection::GetSNPAlleleCounts(unsigned locus, int allele)const{
  int AlleleCounts = 0;
  for(unsigned i = worker_rank; i < size; i += NumWorkers){
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
#ifdef PARALLEL
  MPI::COMM_WORLD.Barrier();
  int totalCounts;
  MPI::COMM_WORLD.Reduce(&AlleleCounts, &totalCounts, 1, MPI::INT, MPI::SUM, 0);
  AlleleCounts = totalCounts;
#endif
  return AlleleCounts;
}

const vector<int> IndividualCollection::getAlleleCounts(unsigned locus, int pop, unsigned NumStates)const{
  int ancestry[2];
  vector<int> counts(NumStates);
  fill(counts.begin(), counts.end(), 0);
  for(unsigned i = worker_rank; i < size; i += NumWorkers)
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
  for(unsigned i = worker_rank; i < size; i += NumWorkers){
    if(_child[i]->GenotypeIsMissing(locus)){
      ++count;
    }
  }
#ifdef PARALLEL
  MPI::COMM_WORLD.Barrier();
  int totalCount;
  MPI::COMM_WORLD.Reduce(&count, &totalCount, 1, MPI::INT, MPI::SUM, 0);
  count = totalCount;
#endif
  return count;
}

// ************** OUTPUT **************

double IndividualCollection::getDevianceAtPosteriorMean(const AdmixOptions* const options, vector<Regression *> &R, Genome* Loci,
							LogWriter &Log, const vector<double>& SumLogRho, unsigned numChromosomes
							, AlleleFreqs* A){
  //TODO: broadcast SumLogRho to workers
  //SumRho = ergodic sum of global sumintensities
  int iterations = options->getTotalSamples()-options->getBurnIn();
  
  //update chromosomes using globalrho, for globalrho model
  if(options->getPopulations() >1 && (options->isGlobalRho() || options->getHapMixModelIndicator()) ){
    vector<double> RhoBar(Loci->GetNumberOfCompositeLoci());
    if(Comms::isMaster())//master only
      for(unsigned i = 0; i < Loci->GetNumberOfCompositeLoci(); ++i)RhoBar[i] = (exp(SumLogRho[i] / (double)iterations));
#ifdef PARALLEL
    if(!Comms::isFreqSampler()) 
      Comms::BroadcastVector(RhoBar);
#endif
    //set locus correlation
    if(Comms::isWorker()){//workers only
      Loci->SetLocusCorrelation(RhoBar);
      if(options->getHapMixModelIndicator())
	for( unsigned int j = 0; j < numChromosomes; j++ )
	  //set global state arrival probs in hapmixmodel
	  //TODO: can skip this if xonly analysis with no females
	  //NB: assumes always diploid in hapmixmodel
	  //KLUDGE: should use global theta as first arg here; Theta in Individual should be the same
	  Loci->getChromosome(j)->SetStateArrivalProbs(options->isRandomMatingModel(), true);
    }
  }
  
  //set haplotype pair probs to posterior means (in parallel version, sets AlleleProbs(Freqs) to posterior means
  if(Comms::isFreqSampler())
    for( unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ )
      (*Loci)(j)->SetHapPairProbsToPosteriorMeans(iterations);
  
#ifdef PARALLEL
  //broadcast allele freqs
  if(!Comms::isMaster())A->BroadcastAlleleFreqs();
#endif
  
  //set genotype probs using happair probs calculated at posterior means of allele freqs 
  if(Comms::isWorker())setGenotypeProbs(Loci, A);
  
  //accumulate deviance at posterior means for each individual
  // fix this to be test individual only if single individual
  double Lhat = 0.0; // Lhat = loglikelihood at estimates
  for(unsigned int i = worker_rank; i < size; i+= NumWorkers ){
    if(options->getHapMixModelIndicator())Lhat += _child[i]->getLogLikelihood(options, false, false);
    else Lhat += _child[i]->getLogLikelihoodAtPosteriorMeans(options);
  }
#ifdef PARALLEL
  if(!Comms::isFreqSampler()){
    Comms::Reduce(&Lhat);
  }
#endif

  if(Comms::isMaster()){
    Log << Quiet << "DevianceAtPosteriorMean(IndAdmixture)" << -2.0*Lhat << "\n";
    for(unsigned c = 0; c < R.size(); ++c){
      double RegressionLogL = R[c]->getLogLikelihoodAtPosteriorMeans(iterations, getOutcome(c));
      Lhat += RegressionLogL;
      Log << "DevianceAtPosteriorMean(Regression " << c+1 << ")"
	  << -2.0*RegressionLogL << "\n";
    }
  }
  return(-2.0*Lhat);
}

double IndividualCollection::getLogLikelihood(const AdmixOptions* const options, bool forceupdate){
  double LogLikelihood = 0.0;
  for(unsigned i = worker_rank; i < size; i+= NumWorkers) {
    LogLikelihood += _child[i]->getLogLikelihood(options, forceupdate, true); // store result if updated
    _child[i]->HMMIsBad(true); // HMM probs overwritten by next indiv, but stored loglikelihood still ok
  }
#ifdef PARALLEL
  //send total to master
  Comms::Reduce(&LogLikelihood);
#endif
  return LogLikelihood;
}

double IndividualCollection::getEnergy(const AdmixOptions* const options, const vector<Regression*> &R, 
				       const bool & annealed) {
  // energy is minus the unnannealed log-likelihood summed over all individuals under study from both HMM and regression 
  // called every iteration after burnin, after update of genotype probs and before annealing
  // accumulates sums of deviance and squared deviance
  double LogLikHMM = 0.0;
  double LogLikRegression = 0.0;
  double Energy = 0.0;
  // assume that HMM probs and stored loglikelihoods are bad, as this function is called after update of allele freqs  
  for(unsigned i = worker_rank; i < size; i+= NumWorkers) {
    LogLikHMM += _child[i]->getLogLikelihood(options, false, !annealed); // store result if not an annealed run
    // don't have to force an HMM update here - on even-numbered iterations with globalrho, stored loglikelihood is still valid
    
    if(annealed)  _child[i]->HMMIsBad(true); // HMM probs bad, stored loglikelihood bad
    else _child[i]->HMMIsBad(false); 
  }
#ifdef PARALLEL
  //send total to master
  Comms::Reduce(&LogLikHMM);
#endif
  // get regression log-likelihood 
  if(Comms::isMaster())
    for(unsigned c = 0; c < R.size(); ++c) LogLikRegression += R[c]->getLogLikelihood(getOutcome(c));
  Energy = -(LogLikHMM + LogLikRegression);
  return Energy;
} 


