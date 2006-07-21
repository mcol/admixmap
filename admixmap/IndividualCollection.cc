/** 
 *   ADMIXMAP
 *   IndividualCollection.cc 
 *   Class to hold an array of Individuals
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "IndividualCollection.h"
#include "StringSplitter.h"
#include "StringConvertor.h"
#include "Regression.h"
#include "admixmap.h"//for MPI stuff
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
  ExpectedY = 0;
  ReportedAncestry = 0;
  SumDeviance = SumDevianceSq = 0.0;
  SumEnergy = 0; SumEnergySq = 0; 
  SumAncestry = 0;
  _child = 0;
  TestInd = 0;
  sizeTestInd = 0;
  indadmixoutput = 0;
  SumLogLikelihood = 0.0;
  //SumResiduals = 0;
  SumLogTheta = 0;
  SumAncestry = 0;
  ReportedAncestry = 0;
  //thetahat = 0;
  //sigma.resize(2);
  //sigma[0] = sigma[1] = 1.0;
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
    GlobalSumAncestry = 0;
    //create communicator for workers and find size of and rank within this group
    workers = MPI::COMM_WORLD.Split( (global_rank>1), global_rank);
    NumWorkers = workers.Get_size();
    if(global_rank >1)
      worker_rank = workers.Get_rank();
    else worker_rank = size;//so that non-workers will not loop over Individuals

    //create communicator for messages between freqsampler and workers
    workers_and_freqs = MPI::COMM_WORLD.Split( (global_rank > 0), global_rank);
    rank_with_freqs = workers_and_freqs.Get_rank();

    workers_and_master = MPI::COMM_WORLD.Split((global_rank!=1), global_rank);
    Populations = options->getPopulations();
#endif

  Individual::SetStaticMembers(Loci, options);
  
  unsigned i0 = 0; // used to offset numbering of other individuals (not the one under test)
  // Fill separate individuals.
  if(options->getTestOneIndivIndicator()) {
    sizeTestInd = options->getNumAnnealedRuns()+1;
    //create array of copies of TestInd
    TestInd = new Individual*[sizeTestInd];
    SumEnergy = new double[sizeTestInd];
    fill(SumEnergy, SumEnergy+sizeTestInd, 0.0);
    SumEnergySq = new double[sizeTestInd];
    fill(SumEnergySq, SumEnergySq+sizeTestInd, 0.0);

    for(int i = 0; i < sizeTestInd; ++i){
      TestInd[i] = new Individual(1, options, Data, true); 
    }     
    ++i0;
    --size;
  }

  if(worker_rank < (int)size){
    _child = new Individual*[size];
    for (unsigned int i = worker_rank; i < size; i += NumWorkers) {
      _child[i] = new Individual(i+i0+1, options, Data, false);//NB: first arg sets Individual's number
    }
  }
}

// ************** DESTRUCTOR **************
IndividualCollection::~IndividualCollection() {
  if(worker_rank==0)
	cout << "\n Deleting individual objects\n" << flush;
  Individual::DeleteStaticMembers();
  for(unsigned int i = worker_rank; i < size; i+=NumWorkers){
    delete _child[i];
  }
  if(TestInd){
    for(int i = 0; i < sizeTestInd; ++i)delete TestInd[i];
  }

  delete[] _child;
  delete[] TestInd;
  delete indadmixoutput;

  delete[] OutcomeType;
  free_matrix(ExpectedY, NumOutcomes);
  //free_matrix(SumResiduals, NumOutcomes);
  //delete[] thetahat;
  delete[] SumLogTheta;
  delete[] SumAncestry;
  delete[] ReportedAncestry;
#ifdef PARALLEL
  delete[] GlobalSumAncestry;
#endif
}
///finish writing expected outcome as R object
void IndividualCollection::FinishWritingEYAsRObject(unsigned NumIterations, const Vector_s Labels){
  //dimensions are NumIndividuals, NumOutcomes, NumIterations
  if(EYStream.is_open()){
    EYStream << ")," << endl << ".Dim = c(" << size << "," << NumOutcomes << "," << NumIterations << ")," << endl
	     << ".Dimnames=list(character(0),c(";
    //write outcome var labels
    for(unsigned j = 0; j < Labels.size(); ++j){
      EYStream << "\"" << Labels[j] << "\"";
      if(j < Labels.size()-1) EYStream << ",";
    }
    EYStream << ") , character(0)))" << endl;
    EYStream.close();  
  }
}
void IndividualCollection::DeleteGenotypes(bool setmissing=false){
  for (unsigned int i = worker_rank; i < size; i += NumWorkers) {
    if(setmissing)_child[i]->SetMissingGenotypes();
    _child[i]->DeleteGenotypes();
  }
}

// ************** INITIALISATION AND LOADING OF DATA **************

void IndividualCollection::Initialise(const AdmixOptions* const options, const Genome* const Loci, const string* const PopulationLabels,
				      const std::vector<std::vector<double> > &alpha, //double rhoalpha, double rhobeta, 
				      LogWriter &Log){
  Log.setDisplayMode(Quiet);
  //Open indadmixture file  
  if ( strlen( options->getIndAdmixtureFilename() ) ){
    Log << "Writing individual-level parameters to " << options->getIndAdmixtureFilename() <<"\n";
    indadmixoutput = new IndAdmixOutputter(options, Loci, PopulationLabels);
  }
  else {
    Log << "No indadmixturefile given\n";
  }

  //Set locusfortest if specified
  if( options->getLocusForTestIndicator() )
    _locusfortest = Loci->GetChrmAndLocus( options->getLocusForTest() );
  
  if(!options->getHapMixModelIndicator()){
    //draw initial values for individual admixture proportions
    for(unsigned i = 0; i < size; i++) _child[i]->drawInitialAdmixtureProps(alpha);
    if(TestInd)for(int i = 0; i < sizeTestInd; i++) TestInd[i]->drawInitialAdmixtureProps(alpha);
  }
  
  // allocate arrays for expected outcome vars and residuals in regression model
  if(NumOutcomes > 0){
    ExpectedY = alloc2D_d(NumOutcomes, NumInd);
    //SumResiduals = alloc2D_d(NumOutcomes, NumInd);
    //for(int j = 0; j < NumOutcomes; ++j)fill(SumResiduals[j], SumResiduals[j] + NumInd, 0.0);
  }
  
  // allocate array of sufficient statistics for update of population admixture parameters
  SumLogTheta = new double[ options->getPopulations()];
  
  //allocate array of sufficient statistics for update of locus-specific sumintensities
  if(options->getHapMixModelIndicator()) {
    SumAncestry = new int[Loci->GetNumberOfCompositeLoci()*2];
#ifdef PARALLEL
    if(workers_and_master.Get_rank()==0)GlobalSumAncestry = new int[Loci->GetNumberOfCompositeLoci()*2];
#endif
  }
}  

void IndividualCollection::LoadData(const AdmixOptions* const options, const InputData* const data_){

  if ( options->getNumberOfOutcomes()>0){
    delete[] OutcomeType;
    OutcomeType = new DataType[ options->getNumberOfOutcomes() ];
    data_->getOutcomeTypes(OutcomeType);

    if(strlen( options->getOutcomeVarFilename() ) != 0)LoadOutcomeVar(data_);
    LoadCovariates(data_, options);
  }
  if ( strlen( options->getReportedAncestryFilename() ) != 0 ){
    LoadRepAncestry(data_);
  }
}

void IndividualCollection::LoadCovariates(const InputData* const data_, const AdmixOptions* const options){
  if ( strlen( options->getCovariatesFilename() ) > 0 ){
    DataMatrix& CovData = (DataMatrix&)data_->getCovariatesMatrix();
    NumberOfInputCovariates = CovData.nCols();
    unsigned NumInds = CovData.nRows();//should already have been checked to be the same as in outcomevarfile

    if( !options->getTestForAdmixtureAssociation() && options->getPopulations() > 1 && !options->getHapMixModelIndicator()){
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
    getLabels(data_->getInputData()[0], CovariateLabels);
    
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
    if( !options->getTestForAdmixtureAssociation() && options->getPopulations() > 1 && !options->getHapMixModelIndicator() ){
      Covariates.setDimensions(NumInds, options->getPopulations());
      for(unsigned i = 0; i < NumInds; ++i)for(int j = 1; j < options->getPopulations(); ++j)
	Covariates.set(i, j, 1.0 / (double)options->getPopulations() );
    } else
      Covariates.setDimensions(NumInds, 1);//just an intercept
    for(unsigned i = 0; i < NumInds; ++i)
      Covariates.set(i,0, 1.0 );
  }
  if( options->getTestForAdmixtureAssociation()  || options->getHapMixModelIndicator())
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

void IndividualCollection::getLabels(const Vector_s& data, Vector_s& labels)
{
  for (size_t i = 0; i < data.size(); ++i) {
    labels.push_back(StringConvertor::dequote(data[i]));
  }
}

void IndividualCollection::SetExpectedY(int k, const double* const beta){
  //sets ExpectedY = X * Beta
  if(ExpectedY){
    matrix_product(Covariates.getData(), beta, ExpectedY[k], Covariates.nRows(), Covariates.nCols(), 1);
    if(OutcomeType[k] == Binary)
      //for binary outcome sets EY as logit^-1(X*beta)
      for(unsigned int i = 0; i < NumInd; i++ )
	ExpectedY[k][i] = 1 / ( 1 + exp( -ExpectedY[k][i] ) );
    }
}

void IndividualCollection::OpenExpectedYFile(const char* Filename, LogWriter & Log){
  if(ExpectedY){
    EYStream.open(Filename, ios::out);
    if( !EYStream.is_open() )
      {
	Log.setDisplayMode(On);
	Log<< "WARNING: Couldn't open expectedoutcomefile\n";
      }
    else{
      Log.setDisplayMode(Quiet);
      Log << "Writing expected values of outcome variable(s) to " << Filename << "\n";
      EYStream << "structure(.Data=c(" << endl;
    }
  }
}
void IndividualCollection::OutputExpectedY(int k){
  //output kth Expected Outcome to file
  if(EYStream.is_open()){
    for(unsigned i = worker_rank; i < size; i+= NumWorkers)
      EYStream << ExpectedY[k][i] << ",";
    EYStream << endl;
  }  
}

void IndividualCollection::HMMIsBad(bool b){
  if(TestInd)    for(int i = 0; i < sizeTestInd; ++i)TestInd[i]->HMMIsBad(b);
  for(unsigned i = worker_rank; i < size; i+= NumWorkers)
    _child[i]->HMMIsBad(b);
}

void IndividualCollection::resetStepSizeApproximators(int k) {
  for(unsigned i = worker_rank; i < size; i+= NumWorkers)
    _child[i]->resetStepSizeApproximator(k);
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
//       //broadcast current values of allele probs to workers 
//       const unsigned NumberOfStates = 2;//(*Loci)(locus)->GetNumberOfStates();
//       double* AlleleProbs;
//       if(rank_with_freqs == 0)AlleleProbs = (double*)(*Loci)(locus)->getAlleleProbs();
//       else AlleleProbs  = new double[NumberOfStates*Populations];
//       workers_and_freqs.Barrier();// wait till everyone is ready for next locus
//       workers_and_freqs.Bcast(AlleleProbs, NumberOfStates*Populations, MPI::DOUBLE, 0);
//get pointer to allele probs for this locus
      const double* AlleleProbs = A->GetAlleleFreqs(locus);//need to get alleleprobs from A as workers have no CompositeLocus objects
      for(unsigned int i = worker_rank; i < size; i+= NumWorkers ) {
	_child[i]->SetGenotypeProbs(j, jj, locus, AlleleProbs);
      }
      if(TestInd)
	for(int i = 0; i < sizeTestInd; ++i)
	  TestInd[i]->SetGenotypeProbs(j, jj, locus, AlleleProbs);
#else
      for(unsigned int i = worker_rank; i < size; i+= NumWorkers ) {
	_child[i]->SetGenotypeProbs(j, jj, locus, false);
      }
      if(TestInd)
	for(int i = 0; i < sizeTestInd; ++i)
	  TestInd[i]->SetGenotypeProbs(j, jj, locus, false);
#endif
      locus++;
    }
  }
}  

void IndividualCollection::annealGenotypeProbs(unsigned nchr, const double coolness, const double* Coolnesses){
  for(unsigned j = 0; j < nchr; ++j){
    
    if(TestInd) { // anneal test individual only
      for(int i = 0; i < sizeTestInd; ++i)
	TestInd[i]->AnnealGenotypeProbs(j, Coolnesses[i]);

    } else { // anneal all individuals
      for(unsigned int i = worker_rank; i < size; i+= NumWorkers) {
	_child[i]->AnnealGenotypeProbs(j, coolness);
      }
    }
  }
}

// ************** UPDATING **************
void IndividualCollection::UpdateIndivAdmixtureRandomWalk(int iteration, const AdmixOptions* const options,
					       const vector<Regression*> &R, const double* const poptheta,
					       const vector<vector<double> > &alpha, 
					       bool anneal=false){
  // Samples individual admixture proportions with random walk 
  bool _anneal = (anneal && !options->getTestOneIndivIndicator());
  vector<double> lambda; // regression precision
  vector<const double*> beta;
  int i0 = 0;
  double dispersion = 0.0; 
  if(options->getTestOneIndivIndicator()) { // anneal likelihood for test individual only 
    i0 = 1;
  }
  if(options->getNumberOfOutcomes() > 0) {
    dispersion = R[0]->getDispersion();
    for(int i = 0; i < options->getNumberOfOutcomes(); ++i){
      lambda.push_back( R[i]->getlambda());
      beta.push_back( R[i]->getbeta());
    }
  }
  
  // loop over individuals
  for(unsigned int i = worker_rank; i < size; i+=NumWorkers ){
    // ** update theta with random-walk update, requiring HMM likelihood
    if(Populations > 1 && !(iteration %2) && !options->getHapMixModelIndicator()) {
      _child[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType, lambda, NumCovariates, &Covariates, 
			     beta, poptheta, options, alpha,DerivativeInverseLinkFunction(i+i0), 
			     dispersion, true, _anneal);
      _child[i]->HMMIsBad(false);

    }
  }
}

/**
    (1) Samples individual admixture proportions on even-numbered iterations
    (2) Samples Locus Ancestry (after updating HMM)
    (3) accumulates sums of ancestry states in hapmixmodel
    (4) Samples Jump Indicators and accumulates sums of (numbers of arrivals) and (ancestry states where there is an arrival)
    (5) updates score, info and score squared for ancestry score tests
    coolness is not passed as argument to this function because annealing has already been implemented by 
    calling annealGenotypeProbs 
*/
void IndividualCollection::SampleLocusAncestry(int iteration, const AdmixOptions* const options,
					       const vector<Regression*> &R, const double* const poptheta,
					       const vector<vector<double> > &alpha, 
					       bool anneal=false){

  //first, some preliminaries
  //
  int Populations = options->getPopulations();
  vector<double> lambda;
  vector<const double*> beta;
  if(!options->getHapMixModelIndicator()){//required only for random walk update of individual admixture and ancestry scoretests
    for(int i = 0; i < options->getNumberOfOutcomes(); ++i){
      lambda.push_back( R[i]->getlambda());
      beta.push_back( R[i]->getbeta());
    }
  }

  //if( !options->getIndAdmixHierIndicator() ) alpha = admixtureprior;

  int i0 = 0;
  if(options->getTestOneIndivIndicator()) {// anneal likelihood for test individual only 
    i0 = 1;
  }

  fill(SumLogTheta, SumLogTheta+options->getPopulations(), 0.0);//reset to 0
  //reset arrays used in score test to 0. This must be done here as the B matrix is updated after sampling admixture
  if(iteration > options->getBurnIn())Individual::ResetScores(options);

  if(options->getHapMixModelIndicator()){
    fill(SumAncestry, SumAncestry + 2*NumCompLoci, 0);
#ifdef PARALLEL
    if(workers_and_master.Get_rank()==0)fill(GlobalSumAncestry, GlobalSumAncestry + 2*NumCompLoci, 0);
    if(worker_rank<(int)size)MPE_Log_event(15, iteration, "Sampleancestry");
#endif
  }
  bool _anneal = (anneal && !options->getTestOneIndivIndicator());
  double dispersion = 0.0; 
  if(!options->getHapMixModelIndicator() && R.size()>0) dispersion = R[0]->getDispersion();

  //now loop over individuals
  for(unsigned int i = worker_rank; i < size; i+=NumWorkers ){

    // ** set SumLocusAncestry and SumNumArrivals to 0
    if(!options->getHapMixModelIndicator())_child[i]->ResetSufficientStats();
    // ** update theta with random-walk proposal on even-numbered iterations
    if(Populations >1 && !(iteration %2) && !options->getHapMixModelIndicator())
      _child[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType, lambda, NumCovariates, &Covariates, 
			     beta, poptheta, options, alpha,DerivativeInverseLinkFunction(i+i0), 
			     dispersion, true, _anneal);
    // ** Run HMM forward recursions and sample locus ancestry
    if(Populations >1)_child[i]->SampleLocusAncestry(options);
    if(options->getHapMixModelIndicator())_child[i]->AccumulateAncestry(SumAncestry);
    // ** Sample JumpIndicators and update SumLocusAncestry and SumNumArrivals
    if(Populations >1 && !options->getHapMixModelIndicator())//no need in hapmixmodel
      _child[i]->SampleJumpIndicators((!options->isGlobalRho() || options->getHapMixModelIndicator()));
    // ** Update score, info and score^2 for ancestry score tests
    if(iteration > options->getBurnIn() && Populations >1 
       && (options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry()))
      _child[i]->UpdateScores(options, &Outcome, OutcomeType, &Covariates, DerivativeInverseLinkFunction(i+i0), 
			      dispersion, ExpectedY);

  }
#ifdef PARALLEL
  if(worker_rank<(int)size)MPE_Log_event(16, iteration, "Sampledancestry");
  if(options->getHapMixModelIndicator()){
    MPE_Log_event(1, iteration, "BarrierStart");
    workers_and_master.Barrier();
    MPE_Log_event(2, iteration, "Barrierend");
    MPE_Log_event(3, iteration, "RedAncStart");    
    workers_and_master.Reduce(SumAncestry, GlobalSumAncestry, 2*NumCompLoci, MPI::INT, MPI::SUM, 0); 
    MPE_Log_event(4, iteration, "RedAncEnd");
  }
#endif
}

void IndividualCollection::SampleLocusAncestry(int iteration, const AdmixOptions* const options,
					       const vector<Regression*> &R) {
  /*
    (1) Samples Locus Ancestry (after updating HMM)
    (2) accumulates sums of ancestry states in hapmixmodel
    (3) Samples Jump Indicators and accumulates sums of (numbers of arrivals) and (ancestry states where there is an arrival)
    (4) updates score, info and score squared for ancestry score tests
  */
  int i0 = 0;
  double dispersion = 0.0; 
  if(options->getTestOneIndivIndicator()) { // anneal likelihood for test individual only 
    i0 = 1;
  }
  if(R.size() > 0) {
    dispersion = R[0]->getDispersion();
  }

  fill(SumLogTheta, SumLogTheta+options->getPopulations(), 0.0);//reset to 0
  //reset arrays used in score test to 0. This must be done here as the B matrix is updated after sampling admixture
  if(iteration > options->getBurnIn())Individual::ResetScores(options);

  if(options->getHapMixModelIndicator()){
    fill(SumAncestry, SumAncestry + 2*NumCompLoci, 0);
#ifdef PARALLEL
    if(workers_and_master.Get_rank()==0)fill(GlobalSumAncestry, GlobalSumAncestry + 2*NumCompLoci, 0);
    if(worker_rank<(int)size)MPE_Log_event(15, iteration, "Sampleancestry");
#endif
  }

  //now loop over individuals
  for(unsigned int i = worker_rank; i < size; i+=NumWorkers ){
    // ** set SumLocusAncestry and SumNumArrivals to 0
    if(!options->getHapMixModelIndicator()) {
      _child[i]->ResetSufficientStats();
    }
    if(Populations > 1) {    // sample locus ancestry: requires forward recursion in HMM
      _child[i]->SampleLocusAncestry(options);
      if(options->getHapMixModelIndicator()) {
	_child[i]->AccumulateAncestry(SumAncestry);
      }

      if(!options->getHapMixModelIndicator() && (!(iteration %2) || !options->isGlobalRho())) { 
	// ** Sample JumpIndicators and update SumLocusAncestry and SumNumArrivals
	// jump indicators required only for conjugate update of theta (even-numbered iterations) or rho
	_child[i]->SampleJumpIndicators( !options->isGlobalRho() );
      }
    }

    // ** Update score, info and score^2 for ancestry score tests
    if(iteration > options->getBurnIn() && Populations >1 
       && (options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry()))
      _child[i]->UpdateScores(options, &Outcome, OutcomeType, &Covariates, DerivativeInverseLinkFunction(i+i0), 
			      dispersion, ExpectedY);
  }

#ifdef PARALLEL
  if(worker_rank<(int)size)MPE_Log_event(16, iteration, "Sampledancestry");
  if(options->getHapMixModelIndicator()){
    MPE_Log_event(1, iteration, "BarrierStart");
    workers_and_master.Barrier();
    MPE_Log_event(2, iteration, "Barrierend");
    MPE_Log_event(3, iteration, "RedAncStart");    
    workers_and_master.Reduce(SumAncestry, GlobalSumAncestry, 2*NumCompLoci, MPI::INT, MPI::SUM, 0); 
    MPE_Log_event(4, iteration, "RedAncEnd");
  }
#endif
}

/**
   Samples Haplotype pairs and upates allele/haplotype counts
*/
void IndividualCollection::SampleHapPairs(const AdmixOptions* const options, AlleleFreqs *A, const Genome* const Loci,
					  bool anneal=false){
  unsigned nchr = Loci->GetNumberOfChromosomes();
  unsigned locus = 0;
  for(unsigned j = 0; j < nchr; ++j){
    for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){
#ifdef PARALLEL
//       //broadcast current values of allele probs to workers 
//       const unsigned NumberOfStates = 2;//(*Loci)(locus)->GetNumberOfStates();
//       double* AlleleProbs;
//       if(rank_with_freqs == 0) AlleleProbs = (double*)(*Loci)(locus)->getAlleleProbs();
//       else AlleleProbs  = new double[NumberOfStates*Populations];
//       workers_and_freqs.Barrier();
//       workers_and_freqs.Bcast(AlleleProbs, NumberOfStates*Populations, MPI::DOUBLE, 0);
//get pointer to allele probs for this locus
      const double* AlleleProbs = A->GetAlleleFreqs(locus);
#endif
      
      // loop over individuals
      //condition for determining which allele freq sampler is in use
      bool condition = (anneal && options->getThermoIndicator() && !options->getTestOneIndivIndicator() );
      for(unsigned int i = worker_rank; i < size; i+=NumWorkers ){
	// ** Sample Haplotype Pair
#ifdef PARALLEL
	_child[i]->SampleHapPair(j, jj, locus, A, options->getHapMixModelIndicator(), condition, AlleleProbs);
#else
	_child[i]->SampleHapPair(j, jj, locus, A, options->getHapMixModelIndicator(), condition);//also updates allele counts
#endif
      }
      locus++;
// #ifdef PARALLEL
//       if(rank_with_freqs >0)delete[] AlleleProbs;
// #endif
    }
  }
}

///samples individual-level sumintensities and admixture
void IndividualCollection::SampleParameters(int iteration, const AdmixOptions* const options,
					    const vector<Regression*> &R, const double* const poptheta,
					    const vector<vector<double> > &alpha, double rhoalpha, double rhobeta,
					    bool anneal=false){
  //sufficient statistics have been stored in Individuals
  // coolness is not passed as argument to this function because annealing has already been implemented by 
  // calling annealGenotypeProbs 
  vector<double> lambda;
  vector<const double*> beta;
  double dispersion = 0.0; 
  if(worker_rank==(int)size || (worker_rank==0 && NumWorkers==1)){//master only
    for(unsigned i = 0; i < R.size(); ++i){
      lambda.push_back( R[i]->getlambda());
      beta.push_back( R[i]->getbeta());
    }
    if(R.size()>0) dispersion = R[0]->getDispersion();
  }

  // ** Test Individuals: all vars updated here as they contribute nothing to score tests or update of allele freqs
  // ---------------------------------------------------------------------------------------------
  int i0 = 0;
  if(options->getTestOneIndivIndicator()) {// anneal likelihood for test individual only 
    i0 = 1;
    for(int i = 0; i < sizeTestInd; ++i){
      // ** set SumLocusAncestry and SumNumArrivals to 0
      TestInd[i]->ResetSufficientStats();
      // ** update theta with random-walk proposal on even-numbered iterations
      if(Populations >1 && !(iteration %2) && !options->getHapMixModelIndicator())
	TestInd[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType, lambda, NumCovariates, &Covariates, 
				beta, poptheta, options, alpha, DerivativeInverseLinkFunction(0), 
				dispersion, true, anneal);
      // ** Run HMM forward recursions and Sample Locus Ancestry
      if(Populations >1)TestInd[i]->SampleLocusAncestry(options);
      // ** Sample JumpIndicators and update SumLocusAncestry and SumNumArrivals
      if(Populations >1)TestInd[i]->SampleJumpIndicators((!options->isGlobalRho() || options->getHapMixModelIndicator()));
      // ** Sample individual- or gamete-specific sumintensities
      if(Populations>1 && !options->getHapMixModelIndicator() && !options->isGlobalRho() ) 
	TestInd[i]->SampleRho( options, rhoalpha, rhobeta,   
			       (!anneal && iteration > options->getBurnIn()));
      // ** update admixture props with conjugate proposal on odd-numbered iterations
      if((iteration %2) && Populations >1 && !options->getHapMixModelIndicator() ) 
	TestInd[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType, lambda, NumCovariates, &Covariates, 
				beta, poptheta, options, alpha, DerivativeInverseLinkFunction(0), dispersion, false, anneal);
      // ** Sample missing values of outcome variable
      TestInd[i]->SampleMissingOutcomes(&Outcome, OutcomeType, ExpectedY, lambda);
      
    }//end loop over test individuals
  }
  // -----------------------------------------------------------------------------------------------------
  
  // ** Non-test individuals - conjugate updates only 
  for(unsigned int i = worker_rank; i < size; i+=NumWorkers ){
    // ** Sample individual- or gamete-specific sumintensities
    if(Populations>1 && !options->getHapMixModelIndicator() && !options->isGlobalRho() ) 
      _child[i]->SampleRho( options, rhoalpha, rhobeta,   
			    (!anneal && iteration > options->getBurnIn()));
    // ** update admixture props with conjugate proposal on odd-numbered iterations
    if((iteration %2) && Populations >1 && !options->getHapMixModelIndicator() ) 
      _child[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType, lambda, NumCovariates, &Covariates, 
			     beta, poptheta, options, alpha, DerivativeInverseLinkFunction(i+i0), dispersion, false, anneal);
    // ** Sample missing values of outcome variable
    _child[i]->SampleMissingOutcomes(&Outcome, OutcomeType, ExpectedY, lambda);
  }
}

void IndividualCollection::setChibNumerator(const AdmixOptions* const options,const vector<vector<double> > &alpha, 
				      double rhoalpha, double rhobeta, AlleleFreqs *A){
  _child[0]->setChibNumerator(options, alpha, rhoalpha, rhobeta, /*thetahat, rhohat*,*/ &MargLikelihood, A);
}

void IndividualCollection::updateChib(const AdmixOptions* const options,const vector<vector<double> > &alpha, 
				      double rhoalpha, double rhobeta, AlleleFreqs *A){
  _child[0]->updateChib(options, alpha, rhoalpha, rhobeta, /*thetahat, rhohat,*/ &MargLikelihood, A);
}

void IndividualCollection::FindPosteriorModes(const AdmixOptions* const options, 
					      const vector<Regression*> &R, 
					      const vector<vector<double> > &alpha, double rhoalpha, double rhobeta,
					      AlleleFreqs* A, 
					      const std::string* const PopulationLabels){
  //TODO: check this for hapmixmodel

  if(options->getDisplayLevel()>1)
    cout<< endl << "Finding posterior mode of individual parameters ..." << endl;
  //open output file and write header
  ofstream modefile(options->getIndAdmixModeFilename());
  modefile << "Individual \t";
  if(!options->isGlobalRho()){
    if(options->isRandomMatingModel())modefile << "rho0 \t rho1 \t";
    else modefile << "rho \t";
  }
  if(options->isRandomMatingModel()){
    for(int g = 0; g < 2; ++g)
      for(int k = 0; k < options->getPopulations(); ++k) modefile << "mu" <<g<<"."<<PopulationLabels[k]<<"\t";
  }
  else{
    for(int k = 0; k < options->getPopulations(); ++k)modefile << "mu"<<PopulationLabels[k]<<"\t";
  }
  modefile <<endl;

  fill(SumLogTheta, SumLogTheta+options->getPopulations(), 0.0);//reset to 0 as mode-finding function changes it
  //may be unecessary if SumLogTheta zeroed after call to this function(FindPosteriorModes)

  vector<double> lambda;
  vector<const double*> beta;
  for(unsigned i = 0; i < R.size(); ++i){
    lambda.push_back( R[i]->getlambda());
    beta.push_back( R[i]->getbeta());
  }
  if(options->getTestOneIndivIndicator()) {// find posterior mode for test individual only 
    TestInd[sizeTestInd-1]->FindPosteriorModes(options, alpha, rhoalpha, rhobeta, A,
					       modefile/*, thetahat, rhohat*/);
  }
  for(unsigned int i = worker_rank; i < size; i+= NumWorkers ){
    _child[i]->FindPosteriorModes(options, alpha, rhoalpha, rhobeta,A,
				  modefile/*, thetahat, rhohat*/);
    modefile << endl;
  }
  modefile.close();
}

// ************** ACCESSORS **************
int IndividualCollection::getSize()const {
  return size;
}

double IndividualCollection::GetSumrho()const
{
   double Sumrho = 0;
   for( unsigned int i = worker_rank; i < size; i+=NumWorkers )
      Sumrho += (*_child[i]).getSumrho();
   return Sumrho;
}

vector<double> IndividualCollection::getOutcome(int j)const{
  return Outcome.getCol(j);
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

const double* IndividualCollection::getCovariates()const{
  return Covariates.getData();
}

const std::string IndividualCollection::getCovariateLabels(int i)const{
  return CovariateLabels[i];
}
const Vector_s IndividualCollection::getCovariateLabels()const{
  return CovariateLabels;
}

double IndividualCollection::getExpectedY(int i)const{
  if(ExpectedY)
    return ExpectedY[0][i];
  else
    return 0.0;
}
double IndividualCollection::getExpectedY(int i, int k)const{
  if(ExpectedY)
    return ExpectedY[k][i];
  else
    return 0.0;
}
double IndividualCollection::getSumLogTheta(int i)const{
  return SumLogTheta[i];
}
const double* IndividualCollection::getSumLogTheta()const{
  return SumLogTheta;
}
const int* IndividualCollection::getSumAncestry()const{
#ifdef PARALLEL
  return GlobalSumAncestry;
#else
  return SumAncestry;
#endif
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
      if(ancestry[0] == pop)++counts[happair[0]];
      if(ancestry[1] == pop)++counts[happair[1]];
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

// returns a reference to the object MargLikelihood of class chib 
const chib* IndividualCollection::getChib()const{
  return &MargLikelihood;
}

///returns sample variance of jth outcome variable
double IndividualCollection::getSampleVarianceOfOutcome(int j)const{
  if(OutcomeType[j] == Continuous){
    double sum = 0.0, sumsq = 0.0, x = 0.0;
    for(unsigned i = 0; i < Outcome.nRows(); ++i){
      x = Outcome.get(i,j);
      sum += x;
      sumsq += x*x;
    }
    return (sumsq - sum*sum / (double)Outcome.nRows()) / (double)Outcome.nRows();
  }
  else return 1.0;
}
///returns sample variance of jth covariate
double IndividualCollection::getSampleVarianceOfCovariate(int j)const{
  if(j < NumberOfInputCovariates+1){
    double sum = 0.0, sumsq = 0.0, x = 0.0;
    for(unsigned i = 0; i < Covariates.nRows(); ++i){
      x = Covariates.get(i,j);
      sum += x;
      sumsq += x*x;
    }
    return (sumsq - sum*sum / (double)Covariates.nRows()) / (double)Covariates.nRows();
  }
  else return 1.0;

}
// ************** OUTPUT **************

double IndividualCollection::getDevianceAtPosteriorMean(const AdmixOptions* const options, vector<Regression *> &R, Genome* Loci,
							LogWriter &Log, const vector<double>& SumLogRho, unsigned numChromosomes
							, AlleleFreqs* A){
#ifdef PARALLEL
  const int rank = MPI::COMM_WORLD.Get_rank();
  //TODO: broadcast SumLogRho to workers
#else
  const int rank = -1;
#endif
  //SumRho = ergodic sum of global sumintensities
  int iterations = options->getTotalSamples()-options->getBurnIn();
  
  if(options->getDisplayLevel()>0)Log.setDisplayMode(On);
  else Log.setDisplayMode(Off);
  
  //update chromosomes using globalrho, for globalrho model
  if(options->getPopulations() >1 && (options->isGlobalRho() || options->getHapMixModelIndicator()) ){
    vector<double> RhoBar(Loci->GetNumberOfCompositeLoci());
    if(rank<1)//master only
    for(unsigned i = 0; i < Loci->GetNumberOfCompositeLoci(); ++i)RhoBar[i] = (exp(SumLogRho[i] / (double)iterations));
#ifdef PARALLEL
    workers_and_master.Barrier();
    workers_and_master.Bcast(&(*(RhoBar.begin())), RhoBar.size(), MPI::DOUBLE, 0);
#endif
    //set locus correlation
    if(rank<0 || rank >1){//workers only
	Loci->SetLocusCorrelation(RhoBar);
	    if(options->getHapMixModelIndicator())
	      for( unsigned int j = 0; j < numChromosomes; j++ )
	//set global state arrival probs in hapmixmodel
	//TODO: can skip this if xonly analysis with no females
	//KLUDGE: should use global theta as first arg here; Theta in Individual should be the same
	      Loci->getChromosome(j)->SetStateArrivalProbs(_child[worker_rank]->getAdmixtureProps(), options->isRandomMatingModel(), true);
						      }
  }
  
  //set haplotype pair probs to posterior means (in parallel version, sets AlleleProbs(Freqs) to posterior means
  if(rank==1 || rank==-1)
  for( unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ )
    (*Loci)(j)->SetHapPairProbsToPosteriorMeans(iterations);

#ifdef PARALLEL
  //broadcast allele freqs
  if(rank!=0)A->BroadcastAlleleFreqs(workers_and_freqs);
#endif
  
  //set genotype probs using happair probs calculated at posterior means of allele freqs 
  if(rank==-1 || rank>1)setGenotypeProbs(Loci, A);
  
  //accumulate deviance at posterior means for each individual
  // fix this to be test individual only if single individual
  double Lhat = 0.0; // Lhat = loglikelihood at estimates
  for(unsigned int i = worker_rank; i < size; i+= NumWorkers ){
    if(options->getHapMixModelIndicator())Lhat += _child[i]->getLogLikelihood(options, false, false);
    else Lhat += _child[i]->getLogLikelihoodAtPosteriorMeans(options);
  }
#ifdef PARALLEL
  if(rank!=1){
    double globalLhat = 0.0;
    workers_and_master.Barrier();
    workers_and_master.Reduce(&Lhat, &globalLhat, 1, MPI::DOUBLE, MPI::SUM, 0);
    Lhat = globalLhat;
  }
#endif

  if(rank <1){
    Log << "DevianceAtPosteriorMean(IndAdmixture)" << -2.0*Lhat << "\n";
    for(unsigned c = 0; c < R.size(); ++c){
      double RegressionLogL = R[c]->getLogLikelihoodAtPosteriorMeans(this, iterations);
      Lhat += RegressionLogL;
      Log << "DevianceAtPosteriorMean(Regression " << c << ")"
	  << -2.0*RegressionLogL << "\n";
    }
  }
  return(-2.0*Lhat);
}

void IndividualCollection::OutputIndAdmixture()
{
  indadmixoutput->visitIndividualCollection(*this);
  if(TestInd)
    indadmixoutput->visitIndividual(*(TestInd[sizeTestInd-1]), _locusfortest);
  for(unsigned int i = worker_rank; i < size; i+=NumWorkers){
    indadmixoutput->visitIndividual(*_child[i], _locusfortest);
  }
}

void IndividualCollection::OutputChibResults(LogWriter& Log) const {
  MargLikelihood.outputResults(Log);
//   Log.setDisplayMode(On);
//   Log << "\nCalculation of Chib algorithm using posterior mode of admixture and sum-intensities, prior mean of allele freqs"
//       << "\nDeviance\t" << -2.0*MargLikelihood.getLogLikelihood()
//       << "\nLogLikelihood\t" << MargLikelihood.getLogLikelihood()
//       << "\nLogPrior\t" << MargLikelihood.getLogPrior()
//       << "\nLogPosterior\t" << MargLikelihood.getLogPosterior()
//       << "\n\nLogMarginalLikelihoodFromChibAlgorithm\t" << MargLikelihood.getLogMarginalLikelihood()
//       << "\n";
} 

void IndividualCollection::getOnePopOneIndLogLikelihood(LogWriter &Log, const string* const PopulationLabels) {
  Log.setDisplayMode(On);
  Log << "Log-likelihood for unadmixed "  << (*PopulationLabels)[0] << ": "
      << _child[0]->getLogLikelihoodOnePop() << "\n";
}

double IndividualCollection::getEnergy(const AdmixOptions* const options, const vector<Regression*> &R, 
				       const bool & annealed) {
  // energy is minus the unnannealed log-likelihood summed over all individuals under study from both HMM and regression 
  // called every iteration after burnin, after update of genotype probs and before annealing
  // accumulates sums of deviance and squared deviance
  double LogLikHMM = 0.0;
  double LogLikRegression = 0.0;
  double Energy = 0.0;
  int global_rank = 0;
  // assume that HMM probs and stored loglikelihoods are bad, as this function is called after update of allele freqs  
  for(unsigned i = worker_rank; i < size; i+= NumWorkers) {
    LogLikHMM += _child[i]->getLogLikelihood(options, false, !annealed); // store result if not an annealed run
    // don't have to force an HMM update here - on even-numbered iterations with globalrho, stored loglikelihood is still valid
    
    if(annealed)  _child[i]->HMMIsBad(true); // HMM probs bad, stored loglikelihood bad
    else _child[i]->HMMIsBad(false); 
  }
#ifdef PARALLEL
  //send total to master
  double globalLogLikHMM = 0.0;
  workers_and_master.Barrier();
  workers_and_master.Reduce(&LogLikHMM, &globalLogLikHMM, 1, MPI::DOUBLE, MPI::SUM, 0);
  LogLikHMM = globalLogLikHMM;
  global_rank = MPI::COMM_WORLD.Get_rank();
#endif
  // get regression log-likelihood 
  if(global_rank==0)
    for(unsigned c = 0; c < R.size(); ++c) LogLikRegression += R[c]->getLogLikelihood(this);
  Energy = -(LogLikHMM + LogLikRegression);
  return Energy;
} 


void IndividualCollection::accumulateEnergyArrays(const AdmixOptions* const options) {
  double Energy = 0.0;
  for(int i = 0; i < sizeTestInd; ++i){ // loop over coolnesses - one copy of test individual at each coolness 
    Energy = -TestInd[i]->getLogLikelihood(options, true, false); // force HMM update, do not store result  
    SumEnergy[i] += Energy;
    SumEnergySq[i] += Energy*Energy;
    TestInd[i]->HMMIsBad(true); // HMM is bad, stored loglikelihood bad
  }
}

double* IndividualCollection::getSumEnergy(){
  return SumEnergy;
}
double* IndividualCollection::getSumEnergySq(){
  return SumEnergySq;
}
///returns Derivative of Inverse Link Function for individual i
double IndividualCollection::DerivativeInverseLinkFunction(int i)const{
  double DInvLink = 1.0;
  if(OutcomeType){//in case no regression model
    if(OutcomeType[0] == Binary)
      DInvLink = ExpectedY[0][i] * (1.0 - ExpectedY[0][i]);
    else if(OutcomeType[0] == Continuous)DInvLink = 1.0;
  }
  return DInvLink;    
}

void IndividualCollection::ResetChib(){
  MargLikelihood.Reset();
}

void IndividualCollection::OutputErgodicChib(std::ofstream *avgstream, bool fixedallelefreqs) {
  *avgstream << MargLikelihood.getLogPrior()<< " " << MargLikelihood.getLogPosterior() << " "
	     << _child[0]->getLogPosteriorTheta() << " " << _child[0]->getLogPosteriorRho()<< " ";
  if(!fixedallelefreqs)
    *avgstream  << _child[0]->getLogPosteriorAlleleFreqs() << " "; // not if fixedallelefreqs
  *avgstream << MargLikelihood.getLogMarginalLikelihood();
}
