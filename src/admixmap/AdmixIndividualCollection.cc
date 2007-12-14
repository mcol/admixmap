/** 
 *   ADMIXMAP
 *   AdmixIndividualCollection.cc
 *   Class to hold an array of AdmixedIndividuals
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "AdmixIndividualCollection.h"
#include "InputAdmixData.h"
#include "bclib/Regression.h"

using namespace std;

// **** CONSTRUCTORS  ****
// AdmixIndividualCollection::AdmixIndividualCollection() {
//     SetNullValues();
// }
void AdmixIndividualCollection::SetNullValues(){
  IndividualCollection::SetNullValues();
  AdmixedChild = 0;
  TestInd = 0;
  sizeTestInd = 0;
  indadmixoutput = 0;
  SumLogTheta = 0;
}

AdmixIndividualCollection::AdmixIndividualCollection(const AdmixOptions* const options, 
						     const InputAdmixData* const Data, Genome* Loci):
  IndividualCollection(Data->getNumberOfIndividuals(), options->getPopulations(), Loci->GetNumberOfCompositeLoci())
 {
  SetNullValues();
  //  Populations = options->getPopulations();
  //NumInd = Data->getNumberOfIndividuals();
  size = NumInd;
  //NumCompLoci = Loci->GetNumberOfCompositeLoci();

  AdmixedIndividual::SetStaticMembers(Loci, options);
  
  unsigned i0 = 0; // used to offset numbering of other individuals (not the one under test)
  // Fill separate individuals.
  if(options->getTestOneIndivIndicator()) {
    sizeTestInd = options->getNumAnnealedRuns()+1;
    //create array of copies of TestInd
    TestInd = new AdmixedIndividual*[sizeTestInd];
    SumEnergy = new double[sizeTestInd];
    fill(SumEnergy, SumEnergy+sizeTestInd, 0.0);
    SumEnergySq = new double[sizeTestInd];
    fill(SumEnergySq, SumEnergySq+sizeTestInd, 0.0);

    for(int i = 0; i < sizeTestInd; ++i){
      TestInd[i] = new AdmixedIndividual(1, options, Data, true); 
    }     
    ++i0;
    --size;
    if(!TestInd[0]->isHaploidIndividual())++NumDiploidIndividuals;
  }

  AdmixedChild = new AdmixedIndividual*[size];
  for (unsigned int i = 0; i < size; i++ ) {
    AdmixedChild[i] = new AdmixedIndividual(i+i0+1, options, Data, false);
    if(!AdmixedChild[i]->isHaploidIndividual())
      ++NumDiploidIndividuals;
  }
  
  _child = (Individual**)AdmixedChild;

}

// ************** DESTRUCTOR **************
AdmixIndividualCollection::~AdmixIndividualCollection() {
  //cout << "\n Deleting individual objects\n" << flush;
  //  AdmixedIndividual::DeleteStaticMembers();

  if(TestInd){
    for(int i = 0; i < sizeTestInd; ++i)delete TestInd[i];
  }

  delete[] TestInd;
  delete indadmixoutput;
}


// ************** INITIALISATION AND LOADING OF DATA **************

void AdmixIndividualCollection::Initialise(const AdmixOptions& options, const Genome& Loci,
				      const Vector_s& PopulationLabels, bclib::LogWriter &Log){
  Log.setDisplayMode(bclib::Quiet);
  //Open indadmixture file  
  if ( strlen( options.getIndAdmixtureFilename() ) ){
    Log << "Writing individual-level parameters to " << options.getIndAdmixtureFilename() <<"\n";
    indadmixoutput = new IndAdmixOutputter(options, Loci, PopulationLabels);
  }
  else {
    Log << "No indadmixturefile given\n";
  }

  //Set locusfortest if specified
  if( options.getLocusForTestIndicator() )
    _locusfortest = Loci.GetChrmAndLocus( options.getLocusForTest() );
  
  // allocate array of sufficient statistics for update of population admixture parameters
  SumLogTheta = new double[ options.getPopulations()];
  
  //allocate array of sufficient statistics for update of locus-specific sumintensities

}  

void AdmixIndividualCollection::DrawInitialAdmixture(const std::vector<std::vector<double> > &alpha){
    //draw initial values for individual admixture proportions
  for(unsigned int i = 0; i < size; i++) 
    AdmixedChild[i]->drawInitialAdmixtureProps(alpha);
    if(TestInd)
      for(int i = 0; i < sizeTestInd; i++) 
	TestInd[i]->drawInitialAdmixtureProps(alpha);

}

void AdmixIndividualCollection::resetStepSizeApproximators(int k) {
  for(unsigned i = 0; i < size; i++)
    AdmixedChild[i]->resetStepSizeApproximator(k);
}

void AdmixIndividualCollection::LoadData(const AdmixOptions* const options, const InputAdmixData* const data_){

  IndividualCollection::LoadData(options, data_, 
				 (!options->getTestForAdmixtureAssociation() && options->getPopulations() > 1));
  if ( strlen( options->getReportedAncestryFilename() ) != 0 ){
    LoadRepAncestry(data_);
  }
}

void AdmixIndividualCollection::LoadRepAncestry(const InputAdmixData* const data_){
  ReportedAncestry = new bclib::DataMatrix[NumInd];
  bclib::DataMatrix& temporary = (bclib::DataMatrix&)data_->getReportedAncestryMatrix();
  for( unsigned i = 0; i < temporary.nRows() / 2; i++ )
    ReportedAncestry[i] = temporary.SubMatrix( 2*i, 2*i + 1, 0, temporary.nCols() - 1 );
 
}

// ************** UPDATING **************
//TODO: ?? have next two functions use those in base class and have ones here only operate on test individuals
void AdmixIndividualCollection::setGenotypeProbs(const Genome* const Loci){
  unsigned nchr = Loci->GetNumberOfChromosomes();
  unsigned locus = 0; // absolute locus number
  for(unsigned j = 0; j < nchr; ++j){
    for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){
      for(unsigned int i = 0; i < size; i++ ) {
	AdmixedChild[i]->SetGenotypeProbs(j, jj, locus, false);
      }
      if(TestInd)
	for(int i = 0; i < sizeTestInd; ++i)
	  TestInd[i]->SetGenotypeProbs(j, jj, locus, false);
      locus++;
    }
  }
}  

void AdmixIndividualCollection::annealGenotypeProbs(unsigned nchr, const double coolness, const double* Coolnesses){
  for(unsigned j = 0; j < nchr; ++j){
    if(TestInd) { // anneal test individual only
      for(int i = 0; i < sizeTestInd; ++i)
	TestInd[i]->AnnealGenotypeProbs(j, Coolnesses[i]);

    } else { // anneal all individuals
      for(unsigned int i = 0; i < size; i++) {
	AdmixedChild[i]->AnnealGenotypeProbs(j, coolness);
      }
    }
  }
}

/**
   (1)Samples admixture with random walk, on even-numbered iterations
   (2) Samples Locus Ancestry (after updating HMM)
   (3) Samples Jump Indicators and accumulates sums of (numbers of arrivals) and (ancestry states where there is an arrival)
   (4) updates score, info and score squared for ancestry score tests
   coolness is not passed as argument to this function because annealing has already been implemented by 
   calling annealGenotypeProbs 
*/
// void AdmixIndividualCollection::SampleAdmixtureWithRandomWalk(int iteration, const AdmixOptions* const options,
// 							 const vector<bclib::Regression*> &R, const double* const poptheta,
// 							 const vector<vector<double> > &alpha, CopyNumberAssocTest& ancestryAssocTest, bool anneal){
//   vector<double> lambda;
//   vector<const double*> beta;
  
//   for(int i = 0; i < options->getNumberOfOutcomes(); ++i){
//     lambda.push_back( R[i]->getlambda());
//     beta.push_back( R[i]->getbeta());
//   }
  
//   double dispersion = 0.0; 
//   if(R.size()>0) dispersion = R[0]->getDispersion();

//   //if( !options->getIndAdmixHierIndicator() ) alpha = admixtureprior;
//   int i0 = 0;
//   if(options->getTestOneIndivIndicator()) {// anneal likelihood for test individual only 
//     i0 = 1;
//   }
//   fill(SumLogTheta, SumLogTheta+options->getPopulations(), 0.0);//reset to 0
//   if(iteration >= options->getBurnIn())
//     ancestryAssocTest.Reset();

//   bool _anneal = (anneal && !options->getTestOneIndivIndicator());
//     for(unsigned int i = 0; i < size; i++ ){
// 	double DinvLink = 1.0;
// 	if(R.size())DinvLink = R[0]->DerivativeInverseLinkFunction(i+i0);
// // ** update theta with random-walk proposal on even-numbered iterations
// 	AdmixedChild[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType, lambda, NumCovariates,
// 			       &Covariates, beta, poptheta, options, alpha, DinvLink, 
// 			       dispersion, ancestryAssocTest, true, _anneal);
//     }
// }

void AdmixIndividualCollection::HMMUpdates(int iteration, const AdmixOptions& options,
					   const vector<bclib::Regression*> &R, const double* const poptheta,
                                                    const vector<vector<double> > &alpha, 
                                                    AffectedsOnlyTest& affectedsOnlyTest, CopyNumberAssocTest& ancestryAssocTest, bool anneal){
  vector<double> lambda;
  vector<const double*> beta;
  double dispersion = 0.0; 
  const bool even_numbered_iteration = !(iteration %2);
  const bool _anneal = (anneal && !options.getTestOneIndivIndicator());
  const bool updateScoreTests = iteration > options.getBurnIn()  
    && (options.getTestForAffectedsOnly() || options.getTestForLinkageWithAncestry());

  if(even_numbered_iteration){//lambda, beta and dispersion are only required for random-walk update of admixture
    for(int i = 0; i < options.getNumberOfOutcomes(); ++i){
      lambda.push_back( R[i]->getlambda());
      beta.push_back( R[i]->getbeta());
    }
    if(R.size()>0) dispersion = R[0]->getDispersion();
  }

  //if( !options.getIndAdmixHierIndicator() ) alpha = admixtureprior;
  int i0 = 0;
  if(options.getTestOneIndivIndicator()) {//skip test individuals when obtaining derivative inverse-link (test indivs are not included in regression)
    i0 = 1;
  }
  //reset sufficient stats for update of pop admixture params to 0
  //  if((iteration %2))
  fill(SumLogTheta, SumLogTheta+options.getPopulations(), 0.0);

  //reset arrays used in score test to 0. This must be done here as the B matrix is updated after sampling admixture
  if(iteration >= options.getBurnIn()){
    ancestryAssocTest.Reset();
    affectedsOnlyTest.Reset();
  }

  //now loop over individuals
  for(unsigned int i = 0; i < size; i++ ){
    // ** set SumLocusAncestry and SumNumArrivals to 0
    AdmixedChild[i]->ResetSufficientStats();

    if(Populations >1){//no updates required for single-population model
      // ** update theta with random-walk proposal on even-numbered iterations
      if(even_numbered_iteration){
        double DinvLink = 1.0;
        if(R.size())DinvLink = R[0]->DerivativeInverseLinkFunction(i+i0);
        AdmixedChild[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType, lambda, NumCovariates,
                                     &Covariates, beta, poptheta, options, alpha, DinvLink, 
                                     dispersion, ancestryAssocTest, true, _anneal);
      }
      
      // ** Run HMM forward recursions, if required, and sample locus ancestry
      _child[i]->SampleHiddenStates(options);
      // ** Sample JumpIndicators and update SumLocusAncestry and SumNumArrivals
      AdmixedChild[i]->SampleJumpIndicators(!options.isGlobalRho());
    
      // ** Update score, info and varscore for ancestry score tests
      if(updateScoreTests)
        AdmixedChild[i]->UpdateScores(options, &Outcome, &Covariates, R, affectedsOnlyTest, ancestryAssocTest);
    }

  }
}

///samples individual-level sumintensities and admixture
void AdmixIndividualCollection::SampleParameters(int iteration, const AdmixOptions& options,
						 const vector<bclib::Regression*> &R, const double* const poptheta,
					    const vector<vector<double> > &alpha, double rhoalpha, double rhobeta,
					    CopyNumberAssocTest& ancestryAssocTest, bool anneal=false){
  //sufficient statistics have been stored in Individuals
  // coolness is not passed as argument to this function because annealing has already been implemented by 
  // calling annealGenotypeProbs 
  vector<double> lambda;
  vector<const double*> beta;
  double dispersion = 0.0; 

  for(unsigned i = 0; i < R.size(); ++i){
    lambda.push_back( R[i]->getlambda());
    beta.push_back( R[i]->getbeta());
  }
  if(R.size()>0) dispersion = R[0]->getDispersion();


  // ** Test Individuals: all vars updated here as they contribute nothing to score tests or update of allele freqs
  // ---------------------------------------------------------------------------------------------
//TODO: move this to separate function
  int i0 = 0;
  if(options.getTestOneIndivIndicator()) {// anneal likelihood for test individual only 
    i0 = 1;
    for(int i = 0; i < sizeTestInd; ++i){
      // ** set SumLocusAncestry and SumNumArrivals to 0
      TestInd[i]->ResetSufficientStats();
      // ** update theta with random-walk proposal on even-numbered iterations
      if(Populations >1 && !(iteration %2))
	TestInd[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType, lambda, NumCovariates, &Covariates, 
				beta, poptheta, options, alpha, 0.0, 
				dispersion, ancestryAssocTest, true, anneal);
      // ** Run HMM forward recursions and Sample Locus Ancestry
      if(Populations >1)TestInd[i]->SampleHiddenStates(options);
      // ** Sample JumpIndicators and update SumLocusAncestry and SumNumArrivals
      if(Populations >1)TestInd[i]->SampleJumpIndicators(!options.isGlobalRho());
      // ** Sample individual- or gamete-specific sumintensities
      if(Populations>1 && !options.isGlobalRho() ) 
	TestInd[i]->SampleRho( options, rhoalpha, rhobeta,   
			       (!anneal && iteration > options.getBurnIn()));
      // ** update admixture props with conjugate proposal on odd-numbered iterations
      if((iteration %2) && Populations >1 ) 
	TestInd[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType, lambda, NumCovariates, &Covariates, 
				beta, poptheta, options, alpha, 0.0, dispersion, ancestryAssocTest, false, anneal);
      // ** Sample missing values of outcome variable
      //TestInd[i]->SampleMissingOutcomes(&Outcome, R);
      
    }//end loop over test individuals
  }
  // -----------------------------------------------------------------------------------------------------
  
  // ** Non-test individuals - conjugate updates only 
  for(unsigned int i = 0; i < size; i++ ){
    // ** Sample individual- or gamete-specific sumintensities
    if(Populations>1 && !options.isGlobalRho() ) 
       AdmixedChild[i]->SampleRho( options, rhoalpha, rhobeta,
				   (!anneal && iteration > options.getBurnIn()));
    // ** update admixture props with conjugate proposal on odd-numbered iterations
     if((iteration %2) && Populations >1 ){
           double DinvLink = 1.0;
       if(R.size())DinvLink = R[0]->DerivativeInverseLinkFunction(i+i0);
        AdmixedChild[i]->SampleTheta(iteration, SumLogTheta, &Outcome, OutcomeType, lambda, NumCovariates, &Covariates,
			       beta, poptheta, options, alpha, DinvLink, dispersion, ancestryAssocTest, false, anneal);
            }
    // ** Sample missing values of outcome variable
    _child[i]->SampleMissingOutcomes(&Outcome, R);
  }
}

void AdmixIndividualCollection::setChibNumerator(const AdmixOptions& options,const vector<vector<double> > &alpha,
				      double rhoalpha, double rhobeta, AlleleFreqs *A){
   AdmixedChild[0]->setChibNumerator(options, alpha, rhoalpha, rhobeta, /*thetahat, rhohat*,*/ &MargLikelihood, A);
}

void AdmixIndividualCollection::updateChib(const AdmixOptions& options,const vector<vector<double> > &alpha,
				      double rhoalpha, double rhobeta, AlleleFreqs *A){
   AdmixedChild[0]->updateChib(options, alpha, rhoalpha, rhobeta, /*thetahat, rhohat,*/ &MargLikelihood, A);
}

void AdmixIndividualCollection::FindPosteriorModes(const AdmixOptions& options,
						   const vector<bclib::Regression*> &R, 
					      const vector<vector<double> > &alpha, double rhoalpha, double rhobeta,
					      AlleleFreqs* A, 
					      const Vector_s& PopulationLabels){
  if(options.getDisplayLevel()>1)
    cout<< endl << "Finding posterior mode of individual parameters ..." << endl;
  //open output file and write header
  ofstream modefile(options.getIndAdmixModeFilename());
  modefile << "Individual \t";
  if(!options.isGlobalRho()){
    if(options.isRandomMatingModel())modefile << "rho0 \t rho1 \t";
    else modefile << "rho \t";
  }
  if(options.isRandomMatingModel()){
    for(int g = 0; g < 2; ++g)
      for(int k = 0; k < options.getPopulations(); ++k) modefile << "mu" <<g<<"."<<PopulationLabels[k]<<"\t";
  }
  else{
    for(int k = 0; k < options.getPopulations(); ++k)modefile << "mu"<<PopulationLabels[k]<<"\t";
  }
  modefile <<endl;

  fill(SumLogTheta, SumLogTheta+options.getPopulations(), 0.0);//reset to 0 as mode-finding function changes it
  //may be unecessary if SumLogTheta zeroed after call to this function(FindPosteriorModes)

  vector<double> lambda;
  vector<const double*> beta;
  for(unsigned i = 0; i < R.size(); ++i){
    lambda.push_back( R[i]->getlambda());
    beta.push_back( R[i]->getbeta());
  }
  if(options.getTestOneIndivIndicator()) {// find posterior mode for test individual only 
    TestInd[sizeTestInd-1]->FindPosteriorModes(options, alpha, rhoalpha, rhobeta, A,
					       modefile/*, thetahat, rhohat*/);
  }
  for(unsigned int i = 0; i < size; i++ ){
     AdmixedChild[i]->FindPosteriorModes(options, alpha, rhoalpha, rhobeta,A,
				  modefile/*, thetahat, rhohat*/);
    modefile << endl;
  }
  modefile.close();
}

// ************** ACCESSORS **************

double AdmixIndividualCollection::GetSumrho()const
{
   double Sumrho = 0;
   for( unsigned int i = 0; i < size; i++ )
      Sumrho += ( AdmixedChild[i])->getSumrho();
   return Sumrho;
}

double AdmixIndividualCollection::getSumLogTheta(int i)const{
  return SumLogTheta[i];
}
const double* AdmixIndividualCollection::getSumLogTheta()const{
  return SumLogTheta;
}

// returns a reference to the object MargLikelihood of class chib
const chib* AdmixIndividualCollection::getChib()const{
  return &MargLikelihood;
}
// ************** OUTPUT **************



void AdmixIndividualCollection::OutputIndAdmixture()
{
  indadmixoutput->visitIndividualCollection(*this);
  if(TestInd)
    indadmixoutput->visitIndividual(*(TestInd[sizeTestInd-1]), _locusfortest);
  for(unsigned int i = 0; i < size; i++){
    indadmixoutput->visitIndividual(* AdmixedChild[i], _locusfortest);
  }
}

void AdmixIndividualCollection::OutputChibResults(bclib::LogWriter& Log) const {
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

void AdmixIndividualCollection::getOnePopOneIndLogLikelihood(bclib::LogWriter &Log, const Vector_s& PopulationLabels) {
  Log.setDisplayMode(bclib::On);
  Log << "Log-likelihood for unadmixed "  << PopulationLabels[0] << ": "
      <<  AdmixedChild[0]->getLogLikelihoodOnePop() << "\n";
}

void AdmixIndividualCollection::accumulateEnergyArrays(const Options& options) {
  double Energy = 0.0;
  for(int i = 0; i < sizeTestInd; ++i){ // loop over coolnesses - one copy of test individual at each coolness 
    Energy = -TestInd[i]->getLogLikelihood(options, true, false); // force HMM update, do not store result  
    SumEnergy[i] += Energy;
    SumEnergySq[i] += Energy*Energy;
    TestInd[i]->HMMIsBad(true); // HMM is bad, stored loglikelihood bad
  }
}

void AdmixIndividualCollection::HMMIsBad(bool b){
  if(TestInd)    
    for(int i = 0; i < sizeTestInd; ++i)
      TestInd[i]->HMMIsBad(b);
  for(unsigned i = 0; i < size; i++)
    _child[i]->HMMIsBad(b);
}

void AdmixIndividualCollection::ResetChib(){
  MargLikelihood.Reset();
}

void AdmixIndividualCollection::OutputErgodicChib(std::ofstream *avgstream, bool fixedallelefreqs) {
  *avgstream << MargLikelihood.getLogPrior()<< " " << MargLikelihood.getLogPosterior() << " "
	     <<  AdmixedChild[0]->getLogPosteriorTheta() << " " <<   AdmixedChild[0]->getLogPosteriorRho()<< " ";
  if(!fixedallelefreqs)
    *avgstream  <<  AdmixedChild[0]->getLogPosteriorAlleleFreqs() << " "; // not if fixedallelefreqs
  *avgstream << MargLikelihood.getLogMarginalLikelihood();
}
double* AdmixIndividualCollection::getSumEnergy()const{
    return SumEnergy;
}
double* AdmixIndividualCollection::getSumEnergySq()const{
    return SumEnergySq;
}
double AdmixIndividualCollection::getDevianceAtPosteriorMean(const Options& options, vector<bclib::Regression *> &R, Genome* Loci,
							     bclib::LogWriter &Log, const vector<double>& SumLogRho, unsigned,  AlleleFreqs* A){
  //SumRho = ergodic sum of global sumintensities
  int iterations = options.getTotalSamples()-options.getBurnIn();
  
  //update chromosomes using globalrho, for globalrho model
  if(options.getPopulations() >1 && options.isGlobalRho() ){
    vector<double> RhoBar(Loci->GetNumberOfCompositeLoci());
    for(unsigned i = 0; i < Loci->GetNumberOfCompositeLoci(); ++i)RhoBar[i] = (exp(SumLogRho[i] / (double)iterations));
    //set locus correlation
    Loci->SetLocusCorrelation(RhoBar);
  }
  
  //set haplotype pair probs to posterior means 
  for( unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ )
    (*Loci)(j)->SetHapPairProbsToPosteriorMeans(iterations);
  
  //set genotype probs using happair probs calculated at posterior means of allele freqs 
  setGenotypeProbs(Loci);

  //accumulate deviance at posterior means for each individual
  // fix this to be test individual only if single individual
  double Lhat = 0.0; // Lhat = loglikelihood at estimates
  for(unsigned int i = 0; i < size; i++ ){
    Lhat += _child[i]->getLogLikelihoodAtPosteriorMeans(options);
  }

  Log << bclib::Quiet << "DevianceAtPosteriorMean(IndAdmixture)" << -2.0*Lhat << "\n";
  for(unsigned c = 0; c < R.size(); ++c){
    double RegressionLogL = R[c]->getLogLikelihoodAtPosteriorMeans(iterations, getOutcome(c));
    Lhat += RegressionLogL;
    Log << "DevianceAtPosteriorMean(Regression " << c+1 << ")"
	<< -2.0*RegressionLogL << "\n";
  }

  return(-2.0*Lhat);
}

#include "AdmixFilenames.h"
///write posterior means of individual admixture params to file
void AdmixIndividualCollection::WritePosteriorMeans(const AdmixOptions& options, const vector<string>& PopLabels)const{
  ofstream meanfile((options.getResultsDir() + "/" + IND_ADMIXTURE_POSTERIOR_MEANS).c_str());
  //set 3 decimal places
  meanfile << std::setfill(' ');
  meanfile.setf(std::ios::fixed); 
  meanfile.precision(3);
  meanfile.width(3);

  //write header
  for(int k = 0; k < options.getPopulations(); ++k){
    meanfile << PopLabels[k];
    if(options.isRandomMatingModel())
      meanfile << "1";
    meanfile << "\t";
  }
  if(options.isRandomMatingModel())
    for(int k = 0; k < options.getPopulations(); ++k){
      meanfile << PopLabels[k] << "2\t";
    }
  if(!options.isGlobalRho()){
    meanfile << "rho";
    if(options.isRandomMatingModel()){
      meanfile << "1\trho2";
    }
  }
  meanfile << endl;

  const unsigned samples = options.getTotalSamples() - options.getBurnIn();
  for(unsigned int i = 0; i < size; i++ ){
    AdmixedChild[i]->WritePosteriorMeans(meanfile, samples, options.isGlobalRho());
    meanfile << endl;
  }
  meanfile.close();
}
