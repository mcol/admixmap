/** 
 *   ADMIXMAP
 *   IndividualCollection.cc 
 *   Class to hold an array of Individuals
 *   Copyright (c) 2002-2006 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#include "IndividualCollection.h"
#include "StringSplitter.h"
#include "Regression.h"

using namespace std;

// **** CONSTRUCTORS  ****
IndividualCollection::IndividualCollection() {
  SumLogTheta = 0;
  OutcomeType = 0;
  CovariateLabels = 0;
  NumCovariates = 0;
  NumberOfInputCovariates = 0;
  ReportedAncestry = 0;
  SumDeviance = SumDevianceSq = 0.0;
}

IndividualCollection::IndividualCollection(const AdmixOptions* const options, const InputData* const Data, const Genome& Loci, 
					   const Chromosome* const* chrm) {
  OutcomeType = 0;
  NumOutcomes = 0;
  NumCovariates = 0;
  NumberOfInputCovariates = 0;
  indadmixoutput = 0;
  SumLogLikelihood = 0.0;
  SumDeviance = SumDevianceSq = 0.0;
  CovariateLabels = 0;
  ExpectedY = 0;
  SumResiduals = 0;
  SumLogTheta = 0;
  ReportedAncestry = 0;
  NumInd = Data->getNumberOfIndividuals();
  size = NumInd;
  NumCompLoci = Loci.GetNumberOfCompositeLoci();
  //sigma.resize(2);
  //sigma[0] = sigma[1] = 1.0;
  TestInd = 0;  // TestInd was declared as a pointer to an Individual object, defaults to 0 (null pointer)
  sizeTestInd = 0;
  SumEnergy = 0; SumEnergySq = 0; 

  Individual::SetStaticMembers(&Loci, options);
  
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
      TestInd[i] = new Individual(1, options, Data, Loci, chrm, true); 
    }     
    ++i0;
    --size;
  }

  _child = new Individual*[size];
  for (unsigned int i = 0; i < size; ++i) {
    _child[i] = new Individual(i+i0+1, options, Data, Loci, chrm, false);//NB: first arg sets Individual's number
    }
}

// ************** DESTRUCTOR **************
IndividualCollection::~IndividualCollection() {
  cout << "Deleting individual objects\n" << flush;
  Individual::DeleteStaticMembers();
  for(unsigned int i = 0; i < size; i++){
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
  free_matrix(SumResiduals, NumOutcomes);
  delete[] SumLogTheta;
  delete[] CovariateLabels;
  delete[] ReportedAncestry;
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
  
  //draw initial values for individual admixture proportions
  for(unsigned i = 0; i < size; i++) _child[i]->drawInitialAdmixtureProps(alpha);
  if(TestInd)for(int i = 0; i < sizeTestInd; i++) TestInd[i]->drawInitialAdmixtureProps(alpha);

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
    SumResiduals = alloc2D_d(NumOutcomes, NumInd);
    for(int j = 0; j < NumOutcomes; ++j)fill(SumResiduals[j], SumResiduals[j] + NumInd, 0.0);
  }
  
  // allocate array of sufficient statistics for update of population admixture parameters
  SumLogTheta = new double[ options->getPopulations()];

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

    if( !options->getTestForAdmixtureAssociation() && options->getPopulations() > 1 ){
      Covariates.setDimensions(NumInd, CovData.nCols() + options->getPopulations());
      for(unsigned i = 0; i < NumInd; ++i)for(int j = 0; j < options->getPopulations()-1; ++j)
         Covariates.set(i, j + CovData.nCols() + 1,  1.0 / (double)options->getPopulations() );
      }
    else
      Covariates.setDimensions(NumInd, CovData.nCols()+1);
    for(unsigned i = 0; i < NumInd; ++i)for(unsigned j = 0; j < CovData.nCols(); ++j)
      {
      Covariates.set(i,0, 1.0); //for intercept
      Covariates.set(i,j+1, CovData.get(i+1, j) );
      Covariates.isMissing(i,j+1, CovData.isMissing(i+1,j));
      }
    if ( Covariates.hasMissing() ) Covariates.SetMissingValuesToColumnMeans();
    CovariateLabels = new string[ NumberOfInputCovariates ];
    getLabels(data_->getInputData()[0], CovariateLabels);

    //vector<double> mean(NumberOfInputCovariates);
    
    //centre covariates about their means
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
  else {
    if( !options->getTestForAdmixtureAssociation() && options->getPopulations() > 1 ){
      Covariates.setDimensions(NumInd, options->getPopulations());
        for(unsigned i = 0; i < NumInd; ++i)for(int j = 1; j < options->getPopulations(); ++j)
          Covariates.set(i, j, 1.0 / (double)options->getPopulations() );
    }
    else
      Covariates.setDimensions(NumInd, 1);
    for(unsigned i = 0; i < NumInd; ++i)
    Covariates.set(i,0, 1.0 );
    }

  if( options->getTestForAdmixtureAssociation() )
       NumCovariates = NumberOfInputCovariates + 1;
  else
       NumCovariates = NumberOfInputCovariates + options->getPopulations();

}

void IndividualCollection::LoadOutcomeVar(const InputData* const data_){
  DataMatrix& OutcomeVarData = (DataMatrix&)data_->getOutcomeVarMatrix();
  NumOutcomes = OutcomeVarData.nCols();
  Outcome = data_->getOutcomeVarMatrix();;
 
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
void IndividualCollection::getLabels(const Vector_s& data, string *labels)
{
  for (size_t i = 0, index = 0; i < data.size(); ++i) {
    labels[index++] = data[i];
  }
}

// shouldn't need this function any more
void IndividualCollection::setAdmixtureProps(const double* const a, size_t thetasize)
{
  if(TestInd)
    for(int i = 0; i < sizeTestInd; ++i)
      TestInd[i]->setAdmixtureProps(a, thetasize);
  for(unsigned int i = 0; i < size; i++){
    _child[i]->setAdmixtureProps(a, thetasize);
  }
}

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

void IndividualCollection::UpdateSumResiduals(){
  if(SumResiduals)
    for(int k = 0; k < NumOutcomes; ++k)
      for(unsigned i = 0; i < NumInd; ++i)
	SumResiduals[k][i] += Outcome.get(i,k) - ExpectedY[k][i];
}

void IndividualCollection::HMMIsBad(bool b){
  if(TestInd)    for(int i = 0; i < sizeTestInd; ++i)TestInd[i]->HMMIsBad(b);
  for(unsigned i = 0; i < size; ++i)
    _child[i]->HMMIsBad(b);
}

void IndividualCollection::resetStepSizeApproximators(int k) {
  for(unsigned i = 0; i < size; ++i)
    _child[i]->resetStepSizeApproximator(k);
}

// ************** UPDATING **************
// pass reference to alpha: either dirichlet parameters from Latent, or admixtureprior set at initialization of this object  
void IndividualCollection::Update(int iteration, const AdmixOptions* const options, Chromosome **chrm, AlleleFreqs *A,
				  const Regression* const R, const double* const poptheta,
				  const std::string* const PopulationLabels,
				  const vector<vector<double> > &alpha, //double globalrho,
				  double rhoalpha, double rhobeta, //LogWriter &Log, 
				  bool anneal=false){
  // coolness is not passed as argument to this function because annealing has already been implemented by 
  // calling annealGenotypeProbs 
  // but we need a similar function to anneal outcome data on each individual for regressions 

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
    for(int i = 0; i < sizeTestInd; ++i){
      TestInd[i]->SampleParameters(SumLogTheta, A, iteration , &Outcome, OutcomeType, ExpectedY,
				lambda, NumCovariates, &Covariates, beta, poptheta, options,
				   chrm, alpha, rhoalpha, rhobeta, //sigma,  
				DerivativeInverseLinkFunction(0),
				R[0].getDispersion(), anneal, true, true, true );
    }
  }
  //next 2 lines go here to prevent test individual contributing to SumLogTheta or sum of scores
  fill(SumLogTheta, SumLogTheta+options->getPopulations(), 0.0);//reset to 0
  if(iteration > options->getBurnIn())Individual::ResetScores(options);

  //posterior modes of individual admixture
  //search at start of burn-in, once annealing is finished
  if(!anneal && iteration == 0 && (options->getChibIndicator() || strlen(options->getIndAdmixModeFilename()))) {
    FindPosteriorModes(options, chrm, A, R, poptheta, alpha, rhoalpha, rhobeta, PopulationLabels);
  }

  for(unsigned int i = 0; i < size; i++ ){
    int prev = i-1;
    if(i==0) prev = size-1;
    cout << flush;
    _child[i]->SampleParameters(SumLogTheta, A, iteration , &Outcome, OutcomeType, ExpectedY,
				lambda, NumCovariates, &Covariates, beta, poptheta, options,
				chrm, alpha, rhoalpha, rhobeta, //sigma,  
				DerivativeInverseLinkFunction(i+i0),
				R[0].getDispersion(), (anneal && !options->getTestOneIndivIndicator()), true, true, true );

    if(size > 1)_child[prev]->HMMIsBad(true); // HMM is bad - see above


    
    if( options->getChibIndicator() && (i == 0) && !anneal ) // if chib option and first individual and not an annealing run
      _child[i]->Chib(iteration, //&SumLogLikelihood, &(MaxLogLikelihood[i]),
		      options, chrm, alpha, //globalrho, 
		      rhoalpha, rhobeta,
		      thetahat, thetahatX, rhohat, rhohatX, &MargLikelihood, A);
  }
}

void IndividualCollection::setGenotypeProbs(Chromosome** C, unsigned nchr){
  if(TestInd)
    for(int i = 0; i < sizeTestInd; ++i)
      for(unsigned j = 0; j < nchr; ++j)
	TestInd[i]->SetGenotypeProbs(j, C[j],false);
  for(unsigned int i = 0; i < size; i++ ) {
    for(unsigned j = 0; j < nchr; ++j)
      _child[i]->SetGenotypeProbs(j, C[j],false);
  }
}  

void IndividualCollection::annealGenotypeProbs(Chromosome** C, unsigned nchr, const double coolness, const double* Coolnesses){
  if(TestInd) { // anneal test individual only
    for(int i = 0; i < sizeTestInd; ++i)
      for(unsigned j = 0; j < nchr; ++j) TestInd[i]->AnnealGenotypeProbs(j, C[j], Coolnesses[i]);

    } else { // anneal all individuals
    for(unsigned int i = 0; i < size; ++i) {
      for(unsigned j = 0; j < nchr; ++j) _child[i]->AnnealGenotypeProbs(j, C[j], coolness);
    }
  }
}

void IndividualCollection::FindPosteriorModes(const AdmixOptions* const options, Chromosome **chrm, AlleleFreqs *A,
					      const Regression* const R, const double* const poptheta,
					      const vector<vector<double> > &alpha, double rhoalpha, double rhobeta, 
					      const std::string* const PopulationLabels){
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
  int i0 = 0;
  if(options->getTestOneIndivIndicator()) {// find posterior mode for test individual only 
    i0 = 1;
    TestInd[sizeTestInd-1]->FindPosteriorModes(SumLogTheta, A, &Outcome, OutcomeType, ExpectedY,
				lambda, NumCovariates, &Covariates, beta, poptheta, options,
					       chrm, alpha, rhoalpha, rhobeta, //sigma,  
				DerivativeInverseLinkFunction(0),
				R[0].getDispersion(), modefile, 
				thetahat, thetahatX, rhohat, rhohatX);
  }
  for(unsigned int i = 0; i < size; i++ ){
    _child[i]->FindPosteriorModes(SumLogTheta, A, &Outcome, OutcomeType, ExpectedY,
				  lambda, NumCovariates, &Covariates, beta, poptheta, options,
				  chrm, alpha, rhoalpha, rhobeta, //sigma,  
				  DerivativeInverseLinkFunction(i+i0),
				  R[0].getDispersion(), modefile, 
				  thetahat, thetahatX, rhohat, rhohatX);
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
   for( unsigned int i = 0; i < size; i++ )
      Sumrho += (*_child[i]).getSumrho();
   return Sumrho;
}

vector<double> IndividualCollection::getOutcome(int j)const{
  return Outcome.getCol(j);
}

double IndividualCollection::getOutcome(int j, int ind)const{
    return Outcome.get(ind, j);
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
const std::string *IndividualCollection::getCovariateLabels()const{
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
const vector<int> IndividualCollection::getSumLocusAncestry(int K)const{
  vector<int> sumlocusancestry(2*K, 0);
  const int* indivsla;
  for(unsigned i = 0; i < size; ++i){
    indivsla = _child[i]->getSumLocusAncestry();
    transform(indivsla, indivsla+2*K, sumlocusancestry.begin(), sumlocusancestry.begin(), std::plus<int>());
  }

  return sumlocusancestry;
}
const vector<int> IndividualCollection::getSumLocusAncestryX(int K)const{
  vector<int> sumlocusancestry(2*K, 0);
  const int* indivsla;
  int length;
  for(unsigned i = 0; i < size; ++i){
    indivsla = _child[i]->getSumLocusAncestryX();
    if(_child[i]->getSex()==female)length = 2*K; else length = K;//two gametes in females, otherwise one
    transform(indivsla, indivsla+length, sumlocusancestry.begin(), sumlocusancestry.begin(), std::plus<int>());
  }

  return sumlocusancestry;
}

const chib* IndividualCollection::getChib()const{
  return &MargLikelihood;
}
// ************** OUTPUT **************

double IndividualCollection::getDevianceAtPosteriorMean(const AdmixOptions* const options, Chromosome** C, Regression *R, 
					  LogWriter &Log, double SumLogRho, unsigned numChromosomes){
  // renamed from OutputDeviance

  //SumRho = ergodic sum of global sumintensities
  int iterations = options->getTotalSamples()-options->getBurnIn();
  
  if(options->getDisplayLevel()>0)Log.setDisplayMode(On);
  else Log.setDisplayMode(Off);
  
  //update chromosomes using globalrho, for globalrho model
  if(options->isGlobalRho())
    for( unsigned int j = 0; j < numChromosomes; j++ )
      C[j]->SetLociCorr(exp(SumLogRho / (double)iterations));
  
  //set haplotype pair probs to posterior means
  for( unsigned int j = 0; j < numChromosomes; j++ )
    for(unsigned jj = 0; jj < C[j]->GetNumberOfCompositeLoci(); ++jj)
      (*C[j])(jj)->SetHapPairProbsToPosteriorMeans(iterations);
  
  //set genotype probs using happair probs calculated at posterior means of allele freqs 
  setGenotypeProbs(C, numChromosomes);
  
  //accumulate deviance at posterior means for each individual
  // fix this to be test individual only if single individual
  double Lhat = 0.0; // Lhat = loglikelihood at estimates
  for(unsigned int i = 0; i < size; i++ ){
    Lhat += _child[i]->getLogLikelihoodAtPosteriorMeans(options, C);
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
  for(unsigned int i = 0; i < size; i++){
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

void IndividualCollection::OutputResiduals(const char* ResidualFilename, const Vector_s Labels, int iterations){
  std::ofstream ResidualStream(ResidualFilename, ios::out);
  if( !ResidualStream )
    {
      cerr<< "WARNING: Couldn't open residualfile\n";
    }
  else{
    for(int j = 0; j < NumOutcomes; ++j)
      ResidualStream << Labels[j]<< "\t";
    ResidualStream << endl;
    for(unsigned i = 0; i < size; ++i){
      for(int j = 0; j < NumOutcomes; ++j)
	ResidualStream << SumResiduals[j][i] / (double) iterations << "\t";
      ResidualStream << endl;
    }
    ResidualStream.close();
  }
}

void IndividualCollection::getOnePopOneIndLogLikelihood(LogWriter &Log, const string* const PopulationLabels) {
  Log.setDisplayMode(On);
  Log << "Log-likelihood for unadmixed "  << (*PopulationLabels)[0] << ": "
      << _child[0]->getLogLikelihoodOnePop() << "\n";
}

double IndividualCollection::getEnergy(const AdmixOptions* const options, Chromosome **C, const Regression* R, 
				       const bool & annealed) {
  // energy is minus the unnannealed log-likelihood summed over all individuals under study from both HMM and regression 
  // called every iteration after burnin, after update of genotype probs and before annealing
  // accumulates sums of deviance and squared deviance
  double LogLikHMM = 0.0;
  double LogLikRegression = 0.0;
  double Energy = 0.0;
  // assume that HMM probs and stored loglikelihoods are bad, as this function is called after update of allele freqs  
  for(unsigned i = 0; i < size; ++i) {
    LogLikHMM += _child[i]->getLogLikelihood(options, C, false, !annealed); // store result if not an annealed run
    // don't have to force an HMM update here - on even-numbered iterations with globalrho, stored loglikelihood is still valid
    
    if(annealed)  _child[i]->HMMIsBad(true); // HMM probs bad, stored loglikelihood bad
    else _child[i]->HMMIsBad(false); 
  }
  // get regression log-likelihood 
  for(int c = 0; c < options->getNumberOfOutcomes(); ++c) LogLikRegression += R[c].getLogLikelihood(this);
  Energy = -(LogLikHMM + LogLikRegression);
  return Energy;
} 


void IndividualCollection::accumulateEnergyArrays(const AdmixOptions* const options, Chromosome **C) {
  double Energy = 0.0;
  for(int i = 0; i < sizeTestInd; ++i){ // loop over coolnesses - one copy of test individual at each coolness 
    Energy = -TestInd[i]->getLogLikelihood(options, C, true, false); // force HMM update, do not store result  
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
