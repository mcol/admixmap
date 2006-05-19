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
  SumLogTheta = 0;
  OutcomeType = 0;
  NumCovariates = 0;
  NumberOfInputCovariates = 0;
  ReportedAncestry = 0;
  SumDeviance = SumDevianceSq = 0.0;
  SumAncestry = 0;
}

IndividualCollection::IndividualCollection(const AdmixOptions* const options, const InputData* const Data, Genome* Loci) {
    worker_rank = 0;
    NumWorkers = 1;
#ifdef PARALLEL
    int global_rank = MPI::COMM_WORLD.Get_rank();
    GlobalSumAncestry = 0;
    //create communicator for workers and find size of and rank within this group
    workers = MPI::COMM_WORLD.Split( (global_rank>1), global_rank);
    NumWorkers = workers.Get_size();
    if(global_rank >1)
      worker_rank =workers.Get_rank();
    else worker_rank = size;//so that non-workers will not loop over Individuals

    //create communicator for messages between freqsampler and workers
    workers_and_freqs = MPI::COMM_WORLD.Split( (global_rank > 0), global_rank);
    rank_with_freqs = workers_and_freqs.Get_rank();

    workers_and_master = MPI::COMM_WORLD.Split((global_rank!=1), global_rank);
    Populations = options->getPopulations();
#endif
  OutcomeType = 0;
  NumOutcomes = 0;
  NumCovariates = 0;
  NumberOfInputCovariates = 0;
  indadmixoutput = 0;
  SumLogLikelihood = 0.0;
  SumDeviance = SumDevianceSq = 0.0;
  ExpectedY = 0;
  //SumResiduals = 0;
  SumLogTheta = 0;
  SumAncestry = 0;
  ReportedAncestry = 0;
  NumInd = Data->getNumberOfIndividuals();
  size = NumInd;
  NumCompLoci = Loci->GetNumberOfCompositeLoci();
  //sigma.resize(2);
  //sigma[0] = sigma[1] = 1.0;
  TestInd = 0;  // TestInd was declared as a pointer to an Individual object, defaults to 0 (null pointer)
  sizeTestInd = 0;
  SumEnergy = 0; SumEnergySq = 0; 

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

  _child = new Individual*[size];
  for (unsigned int i = worker_rank; i < size; i += NumWorkers) {
    _child[i] = new Individual(i+i0+1, options, Data, false);//NB: first arg sets Individual's number
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
  cout << flush;
  delete indadmixoutput;

  delete[] OutcomeType;
  free_matrix(ExpectedY, NumOutcomes);
  //free_matrix(SumResiduals, NumOutcomes);
  delete[] SumLogTheta;
  delete[] SumAncestry;
  delete[] ReportedAncestry;
#ifdef PARALLEL
  delete[] GlobalSumAncestry;
#endif
}

void IndividualCollection::FinishWritingEYAsRObject(unsigned NumIterations, const Vector_s Labels){
  //finish writing expected outcome as R object
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
				      const std::vector<std::vector<double> > &alpha, double rhoalpha, double rhobeta, 
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

//   // set priors on individual admixture to be passed to individual if no hierarchical model
//   admixtureprior.resize(2);
//   admixtureprior[0].resize(K);
//   admixtureprior[1].resize(K);
//   fill(admixtureprior[0].begin(), admixtureprior[0].end(), 1.0); // set default of flat prior on individual admixture
//   fill(admixtureprior[1].begin(), admixtureprior[1].end(), 1.0);
//   if( !options->getIndAdmixHierIndicator() && options->sizeInitAlpha() > 0 ) { 
//     admixtureprior[0] = options->getInitAlpha(0); // returns a vector of doubles of length K
//     admixtureprior[1] = options->getInitAlpha(options->sizeInitAlpha() - 1);
//   } // if hierarchical model, no need to set prior on individual admixture
  
  // allocate arrays for expected outcome vars and residuals in regression model
  if(options->getNumberOfOutcomes() > 0){
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

//   allocate and set initial values for estimates used in Chib algorithm
  if( options->getChibIndicator() )
    InitialiseMLEs(rhoalpha,rhobeta,options);
  //set to very large negative value (effectively -Inf) so the first value is guaranteed to be greater
  MaxLogLikelihood.assign(NumInd, -9999999 );
}


// // ** this function needs debugging
// required for chib algorithm but all necessary code could be moved to individual 
void IndividualCollection::InitialiseMLEs(double rhoalpha, double rhobeta, const AdmixOptions* const options){
  //set thetahat and rhohat, estimates of individual admixture and sumintensities
   size_t size_admix;
   int K = options->getPopulations();
   if( options->isRandomMatingModel() )
     size_admix = K*2; // double the size for 2 gametes in RMM
   else//assortative mating
     size_admix = K;

   thetahat = new double[size_admix];
   thetahatX = new double[size_admix];

   //initialise thetahat at initial values of individual admixture
   for(unsigned k = 0; k < size_admix; ++k)
     thetahat[k] = _child[0]->getAdmixtureProps()[k];

   if(!options->isGlobalRho())
     //initialise rhohat at initial value of Individual 1's sumintensities
      rhohat = _child[0]->getRho();
   else {
   //initialise rhohat at initial value of globalsumintensities ie prior mean
     vector<double> r(2, rhoalpha/rhobeta );
     rhohat = r;
     rhohatX = r;
   }
   //TODO: X objects
}

void IndividualCollection::LoadData(const AdmixOptions* const options, const InputData* const data_){
  LoadCovariates(data_, options);
  
  if ( strlen( options->getOutcomeVarFilename() ) != 0 ){
    LoadOutcomeVar(data_);
  }
  if ( strlen( options->getReportedAncestryFilename() ) != 0 ){
    LoadRepAncestry(data_);
  }

}

void IndividualCollection::LoadCovariates(const InputData* const data_, const AdmixOptions* const options){
  if ( strlen( options->getCovariatesFilename() ) > 0 ){
    DataMatrix& CovData = (DataMatrix&)data_->getCovariatesMatrix();
    NumberOfInputCovariates = CovData.nCols();

    if( !options->getTestForAdmixtureAssociation() && options->getPopulations() > 1 && !options->getHapMixModelIndicator()){
      Covariates.setDimensions(NumInd, CovData.nCols() + options->getPopulations());
      for(unsigned i = 0; i < NumInd; ++i)for(int j = 0; j < options->getPopulations()-1; ++j)
         Covariates.set(i, j + CovData.nCols() + 1,  1.0 / (double)options->getPopulations() );
      }
    else
      Covariates.setDimensions(NumInd, CovData.nCols()+1);//+1 for intercept
    for(unsigned i = 0; i < NumInd; ++i)for(unsigned j = 0; j < CovData.nCols(); ++j)
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
    //for(unsigned int i = 0; i < NumInd; i++ )
    //if(!Covariates.IsMissingValue(i,j+1)){
    // mean[j] += Covariates(i,j+1);
    //++count;
    //}
    //mean[j] /= (double)count;
    //}
    
    //for(unsigned int i = 0; i < NumInd; i++ )
    // for( int j = 0; j < NumberOfInputCovariates; j++ )
    //Covariates( i, j+1 ) -= mean[j];
  }
  else {//no covariatesfile
    if( !options->getTestForAdmixtureAssociation() && options->getPopulations() > 1 && !options->getHapMixModelIndicator() ){
      Covariates.setDimensions(NumInd, options->getPopulations());
        for(unsigned i = 0; i < NumInd; ++i)for(int j = 1; j < options->getPopulations(); ++j)
          Covariates.set(i, j, 1.0 / (double)options->getPopulations() );
    }
    else
      Covariates.setDimensions(NumInd, 1);//just an intercept
    for(unsigned i = 0; i < NumInd; ++i)
    Covariates.set(i,0, 1.0 );
    }

  if( options->getTestForAdmixtureAssociation()  || options->getHapMixModelIndicator())
       NumCovariates = NumberOfInputCovariates + 1;
  else
       NumCovariates = NumberOfInputCovariates + options->getPopulations();

}

void IndividualCollection::LoadOutcomeVar(const InputData* const data_){
  Outcome = data_->getOutcomeVarMatrix();
  NumOutcomes = Outcome.nCols();
 
  delete[] OutcomeType;
  OutcomeType = new DataType[ NumOutcomes ];
  data_->getOutcomeTypes(OutcomeType);
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

// void IndividualCollection::setAdmixtureProps(const double* const a, size_t thetasize)
// {
//   if(TestInd)
//     for(int i = 0; i < sizeTestInd; ++i)
//       TestInd[i]->setAdmixtureProps(a, thetasize);
//   for(unsigned int i = 0; i < size; i++){
//     _child[i]->setAdmixtureProps(a, thetasize);
//   }
// }

// void IndividualCollection::setAdmixturePropsX(const double* const a, size_t thetasize)
// {
//   if(TestInd)
//     for(int i = 0; i < sizeTestInd; ++i)
//       TestInd[i]->setAdmixtureProps(a, thetasize);
//   for(unsigned int i = 0; i < size; i++){
//     _child[i]->setAdmixturePropsX(a, thetasize);
//   }
// }

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
  EYStream.open(Filename, ios::out);
  if( !EYStream.is_open() )
    {
      Log<< "WARNING: Couldn't open expectedoutcomefile\n";
    }
  else{
    Log << "Writing expected values of outcome variable(s) to " << Filename << "\n";
    EYStream << "structure(.Data=c(" << endl;
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

// void IndividualCollection::UpdateSumResiduals(){
//   if(SumResiduals)
//     for(int k = 0; k < NumOutcomes; ++k)
//       for(unsigned i = 0; i < NumInd; ++i)
// 	SumResiduals[k][i] += Outcome.get(i,k) - ExpectedY[k][i];
// }
// void IndividualCollection::OutputResiduals(const char* ResidualFilename, const Vector_s Labels, int iterations){
//   std::ofstream ResidualStream(ResidualFilename, ios::out);
//   if( !ResidualStream )
//     {
//       cerr<< "WARNING: Couldn't open residualfile\n";
//     }
//   else{
//     for(int j = 0; j < NumOutcomes; ++j)
//       ResidualStream << Labels[j]<< "\t";
//     ResidualStream << endl;
//     for(unsigned i = 0; i < size; ++i){
//       for(int j = 0; j < NumOutcomes; ++j)
// 	ResidualStream << SumResiduals[j][i] / (double) iterations << "\t";
//       ResidualStream << endl;
//     }
//     ResidualStream.close();
//   }
// }


void IndividualCollection::HMMIsBad(bool b){
  if(TestInd)    for(int i = 0; i < sizeTestInd; ++i)TestInd[i]->HMMIsBad(b);
  for(unsigned i = worker_rank; i < size; i+= NumWorkers)
    _child[i]->HMMIsBad(b);
}

void IndividualCollection::resetStepSizeApproximators(int k) {
  for(unsigned i = worker_rank; i < size; i+= NumWorkers)
    _child[i]->resetStepSizeApproximator(k);
}

void IndividualCollection::setGenotypeProbs(const Genome* const Loci){
  unsigned nchr = Loci->GetNumberOfChromosomes();
  unsigned locus = 0;
  for(unsigned j = 0; j < nchr; ++j){
    for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){
#ifdef PARALLEL
      //broadcast current values of allele probs to workers 
      const unsigned NumberOfStates = (*Loci)(locus)->GetNumberOfStates();
      double* AlleleProbs;
      if(rank_with_freqs == 0)AlleleProbs = (double*)(*Loci)(locus)->getAlleleProbs();
      else AlleleProbs  = new double[NumberOfStates*Populations];
      workers_and_freqs.Barrier();// wait till everyone is ready for next locus
      workers_and_freqs.Bcast(AlleleProbs, NumberOfStates*Populations, MPI::DOUBLE, 0);

      for(unsigned int i = worker_rank; i < size; i+= NumWorkers ) {
	_child[i]->SetGenotypeProbs(j, jj, locus, AlleleProbs);
      }
      if(TestInd)
	for(int i = 0; i < sizeTestInd; ++i)
	  TestInd[i]->SetGenotypeProbs(j, jj, locus, AlleleProbs);
      if(rank_with_freqs >0)delete[] AlleleProbs;
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
void IndividualCollection::SampleLocusAncestry(int iteration, const AdmixOptions* const options,
					       const Regression* const R, const double* const poptheta,
					       const vector<vector<double> > &alpha, 
					       bool anneal=false){
  /*
    (1) Samples individual admixture proportions on even-numbered iterations
    (2) Samples Locus Ancestry (after updating HMM)
    (3) accumulates sums of ancestry states in hapmixmodel
    (4) Samples Jump Indicators and accumulates sums of (numbers of arrivals) and (ancestry states where there is an arrival)
    (5) updates score, info and score squared for ancestry score tests
    coolness is not passed as argument to this function because annealing has already been implemented by 
    calling annealGenotypeProbs 
*/

  //first, some preliminaries
  //
  int Populations = options->getPopulations();
  vector<double> lambda;
  vector<const double*> beta;
  for(int i = 0; i < options->getNumberOfOutcomes(); ++i){
    lambda.push_back( R[i].getlambda());
    beta.push_back( R[i].getbeta());
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
#endif
  }
  bool _anneal = (anneal && !options->getTestOneIndivIndicator());
  double dispersion = 0.0; 
  if(options->getNumberOfOutcomes()>0) dispersion = R[0].getDispersion();

  //now loop over individuals
  for(unsigned int i = worker_rank; i < size; i+=NumWorkers ){

    // ** set SumLocusAncestry and SumNumArrivals to 0
    _child[i]->ResetSufficientStats();
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
  if(options->getHapMixModelIndicator()){
    MPE_Log_event(11, iteration, "Barrier");
    workers_and_master.Barrier();
    MPE_Log_event(12, iteration, "BarrEnd");
    
    MPE_Log_event(3, iteration, "RedAncStart");
    workers_and_master.Reduce(SumAncestry, GlobalSumAncestry, 2*NumCompLoci, MPI::INT, MPI::SUM, 0); 
    MPE_Log_event(4, iteration, "RedAncEnd");
  }
#endif
}
void IndividualCollection::SampleHapPairs(const AdmixOptions* const options, AlleleFreqs *A, const Genome* const Loci,
					  bool anneal=false){
  /*
    Samples Haplotype pairs and upates allele/haplotype counts
  */
  unsigned nchr = Loci->GetNumberOfChromosomes();
  unsigned locus = 0;
  for(unsigned j = 0; j < nchr; ++j){
    for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){
#ifdef PARALLEL
      //broadcast current values of allele probs to workers 
      const unsigned NumberOfStates = (*Loci)(locus)->GetNumberOfStates();
      double* AlleleProbs;
      if(rank_with_freqs == 0) AlleleProbs = (double*)(*Loci)(locus)->getAlleleProbs();
      else AlleleProbs  = new double[NumberOfStates*Populations];
      workers_and_freqs.Barrier();
      workers_and_freqs.Bcast(AlleleProbs, NumberOfStates*Populations, MPI::DOUBLE, 0);
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
#ifdef PARALLEL
      if(rank_with_freqs >0)delete[] AlleleProbs;
#endif
    }
  }
}

void IndividualCollection::SampleParameters(int iteration, const AdmixOptions* const options,
					    const Regression* const R, const double* const poptheta,
					    const vector<vector<double> > &alpha, double rhoalpha, double rhobeta,
					    bool anneal=false){
  //samples individual-level sumintensities and admixture
  //sufficient statistics have been stored in Individuals

  // coolness is not passed as argument to this function because annealing has already been implemented by 
  // calling annealGenotypeProbs 
  int Populations = options->getPopulations();
  vector<double> lambda;
  vector<const double*> beta;
  for(int i = 0; i < options->getNumberOfOutcomes(); ++i){
    lambda.push_back( R[i].getlambda());
    beta.push_back( R[i].getbeta());
  }
  double dispersion = 0.0; 
  if(options->getNumberOfOutcomes()>0) dispersion = R[0].getDispersion();
  // ** Test Individuals
  // these are updated completely here as they contribute nothing to the score tests or update of allele freqs
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
  // ** Non-test individuals
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

void IndividualCollection::UpdateChib(int iteration, const AdmixOptions* const options,const vector<vector<double> > &alpha, 
				      double rhoalpha, double rhobeta, AlleleFreqs *A){
    _child[0]->Chib(iteration, //&SumLogLikelihood, &(MaxLogLikelihood[i]),
		    options, alpha, //globalrho, 
		    rhoalpha, rhobeta,
		    thetahat, thetahatX, rhohat, rhohatX, &MargLikelihood, A);

}

void IndividualCollection::FindPosteriorModes(const AdmixOptions* const options, 
					      const Regression* const R, 
					      const vector<vector<double> > &alpha, double rhoalpha, double rhobeta, 
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
  for(int i = 0; i < options->getNumberOfOutcomes(); ++i){
    lambda.push_back( R[i].getlambda());
    beta.push_back( R[i].getbeta());
  }
  if(options->getTestOneIndivIndicator()) {// find posterior mode for test individual only 
    TestInd[sizeTestInd-1]->FindPosteriorModes(options, alpha, rhoalpha, rhobeta, 
					       modefile, thetahat, thetahatX, rhohat, rhohatX);
  }
  for(unsigned int i = worker_rank; i < size; i+= NumWorkers ){
    _child[i]->FindPosteriorModes(options, alpha, rhoalpha, rhobeta,
				   modefile, thetahat, thetahatX, rhohat, rhohatX);
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
// const vector<int> IndividualCollection::getSumLocusAncestry(int K)const{
//   vector<int> sumlocusancestry(2*K, 0);
//   const int* indivsla;
//   for(unsigned i = 0; i < size; ++i){
//     indivsla = _child[i]->getSumLocusAncestry();
//     transform(indivsla, indivsla+2*K, sumlocusancestry.begin(), sumlocusancestry.begin(), std::plus<int>());
//   }

//   return sumlocusancestry;
// }
// const vector<int> IndividualCollection::getSumLocusAncestryX(int K)const{
//   vector<int> sumlocusancestry(2*K, 0);
//   const int* indivsla;
//   int length;
//   for(unsigned i = 0; i < size; ++i){
//     indivsla = _child[i]->getSumLocusAncestryX();
//     if(_child[i]->getSex()==female)length = 2*K; else length = K;//two gametes in females, otherwise one
//     transform(indivsla, indivsla+length, sumlocusancestry.begin(), sumlocusancestry.begin(), std::plus<int>());
//   }

//   return sumlocusancestry;
// }

// const vector<unsigned> IndividualCollection::getSumNumArrivals(){
//   //returns a vector of sums over gametes of numbers of arrivals between pairs of loci
//   vector<unsigned>sumNumArrivals(NumCompLoci, 0);
//   for(unsigned i = 0; i < size; ++i)
//     _child[i]->getSumNumArrivals(&sumNumArrivals);
//   return sumNumArrivals;
// }

unsigned IndividualCollection::GetSNPAlleleCounts(unsigned locus, int allele)const{
  /**
   * returns a count of the copies of allele a at a comp locus
   * only works for diallelic loci.
   * used in strat test to determine loci with given percentage of observed genotypes 
   */

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

int IndividualCollection::getNumberOfMissingGenotypes(unsigned locus)const{
  //count number of missing genotypes at locus
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

const chib* IndividualCollection::getChib()const{
  return &MargLikelihood;
}

double IndividualCollection::getSampleVarianceOfOutcome(int j)const{
  //returns sample variance of jth outcome variable
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
double IndividualCollection::getSampleVarianceOfCovariate(int j)const{
  //returns sample variance of jth covariate
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

double IndividualCollection::getDevianceAtPosteriorMean(const AdmixOptions* const options, Regression *R, Genome* Loci,
							LogWriter &Log, const vector<double>& SumLogRho, unsigned numChromosomes){
  // renamed from OutputDeviance

  //SumRho = ergodic sum of global sumintensities
  int iterations = options->getTotalSamples()-options->getBurnIn();
  
  if(options->getDisplayLevel()>0)Log.setDisplayMode(On);
  else Log.setDisplayMode(Off);
  
  //update chromosomes using globalrho, for globalrho model
  if(options->getPopulations() >1 && (options->isGlobalRho() || options->getHapMixModelIndicator()) ){
    vector<double> RhoBar;
    for(vector<double>::const_iterator i = SumLogRho.begin(); i < SumLogRho.end(); ++i)RhoBar.push_back(exp(*i / (double)iterations));
    //set locus correlation
    Loci->SetLocusCorrelation(RhoBar);
    if(options->getHapMixModelIndicator())
      for( unsigned int j = 0; j < numChromosomes; j++ )
	//set global state arrival probs in hapmixmodel
	//TODO: can skip this if xonly analysis with no females
	//KLUDGE: should use global theta as first arg here; Theta in Individual should be the same
	Loci->getChromosome(j)->SetStateArrivalProbs(_child[0]->getAdmixtureProps(), options->isRandomMatingModel(), true);
  }
  
  //set haplotype pair probs to posterior means
  for( unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ )
    (*Loci)(j)->SetHapPairProbsToPosteriorMeans(iterations);
  
  //set genotype probs using happair probs calculated at posterior means of allele freqs 
  setGenotypeProbs(Loci);
  
  //accumulate deviance at posterior means for each individual
  // fix this to be test individual only if single individual
  double Lhat = 0.0; // Lhat = loglikelihood at estimates
  for(unsigned int i = worker_rank; i < size; i+= NumWorkers ){
    if(options->getHapMixModelIndicator())Lhat += _child[i]->getLogLikelihood(options, false, false);
    else Lhat += _child[i]->getLogLikelihoodAtPosteriorMeans(options);
  }
  Log << "DevianceAtPosteriorMean(IndAdmixture)" << -2.0*Lhat << "\n";
  for(int c = 0; c < options->getNumberOfOutcomes(); ++c){
    double RegressionLogL = R[c].getLogLikelihoodAtPosteriorMeans(this, iterations);
    Lhat += RegressionLogL;
    Log << "DevianceAtPosteriorMean(Regression " << c << ")"
	<< -2.0*RegressionLogL << "\n";
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

void IndividualCollection::OutputChibEstimates(bool RandomMating, LogWriter &Log, int Populations)const{
  //Used only if chib = 1
  Log.setDisplayMode(Off);
  Log << "Parameter Values used for Chib Algorithm\t";

  for(int k = 0; k < Populations; ++k){
    Log << thetahat[k] << "\t";
  }
  if(RandomMating)
  for(int k = 0; k < Populations; ++k){
    Log << thetahat[Populations +k] << "\t";
  }
  Log << rhohat[0];
  if(RandomMating)Log <<"\t" << rhohat[1];
  Log << "\n";
}

void IndividualCollection::OutputChibResults(LogWriter& Log)const{
  Log.setDisplayMode(On);
  Log << "\nCalculation of Chib algorithm using posterior mode of admixture and sum-intensities, prior mean of allele freqs"
      << "\nDeviance\t" << -2.0*MargLikelihood.getLogLikelihood()
      << "\nLogLikelihood\t" << MargLikelihood.getLogLikelihood()
      << "\nLogPrior\t" << MargLikelihood.getLogPrior()
      << "\nLogPosterior\t" << MargLikelihood.getLogPosterior()
      << "\n\nLogMarginalLikelihoodFromChibAlgorithm\t" << MargLikelihood.getLogMarginalLikelihood()
      << "\n";
} 

void IndividualCollection::getOnePopOneIndLogLikelihood(LogWriter &Log, const string* const PopulationLabels) {
  Log.setDisplayMode(On);
  Log << "Log-likelihood for unadmixed "  << (*PopulationLabels)[0] << ": "
      << _child[0]->getLogLikelihoodOnePop() << "\n";
}

double IndividualCollection::getEnergy(const AdmixOptions* const options, const Regression* R, 
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
  // get regression log-likelihood 
  for(int c = 0; c < options->getNumberOfOutcomes(); ++c) LogLikRegression += R[c].getLogLikelihood(this);
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
//returns Derivative of Inverse Link Function for individual i
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

void IndividualCollection::OutputErgodicChib(std::ofstream *avgstream) {
  *avgstream << MargLikelihood.getLogPrior()<< " " << MargLikelihood.getLogPosterior() << " "
	     << _child[0]->getLogPosteriorTheta() << " " << _child[0]->getLogPosteriorRho()<< " " 
	     << _child[0]->getLogPosteriorAlleleFreqs() << " "
	     << MargLikelihood.getLogMarginalLikelihood();
}
