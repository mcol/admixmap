/* 
 *   HAPMIXMAP
 *   HapMixModel.h
 *   
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "HapMixModel.h"
#include "HapMixIndividual.h"
#include "HapMixFilenames.h"
#include "bclib/LogWriter.h"

#define SCORETEST_UPDATE_EVERY 2
int numdiploidIndivs = 0;

using namespace bclib;

HapMixModel::HapMixModel(){
  L = 0;
  HMIC = 0;
  IC = 0; 
}

HapMixModel::~HapMixModel(){
  delete L;
}

void HapMixModel::Initialise(HapMixOptions& options, InputHapMixData& data,  LogWriter& Log){

  InitialiseGenome(Loci, options, data, Log);
  //check allelefreq files, initialises allele freqs and finishes setting up Composite Loci
  A.Initialise(&options, &data, &Loci, Log); 
  pA = &A;//set pointer to AlleleFreqs object

  L = new PopHapMix(options, Loci);
  L->Initialise(data.getUnitOfDistanceAsString(), Log);

  //create HapMixIndividualCollection object
  //NB call after A Initialise, and before R Initialise
  HMIC = new HapMixIndividualCollection(&options, &data, &Loci, L->getGlobalMixtureProps());
  //set IC pointer in base class to the same address as HMIC
  IC = HMIC;

  //load Outcome and Covariate data into IC
  IC->LoadData(&options, &data, false);

  // set unannealed probs
  //IC->setGenotypeProbs(&Loci, &A);


  //get number of diploid individuals by summing over individuals
  numdiploidIndivs = IC->getNumDiploidIndividuals();
  Loci.SetHMMDimensions(options.getNumberOfBlockStates(), (bool)(numdiploidIndivs !=0));

  //if there are any diploid individuals, allocate space for diploid genotype probs
  if(numdiploidIndivs){
    A.AllocateDiploidGenotypeProbs();
    A.SetDiploidGenotypeProbs();
  }
  A.setSampler(options.getThermoIndicator(), 
	       ( !strlen(options.getPriorAlleleFreqFilename()) 
		 && !(options.getInitialAlleleFreqFilename().size()) ));

  HapMixIndividual::SetGenotypeProbs(&Loci, A.getHaploidGenotypeProbs(), A.getDiploidGenotypeProbs());

  const int numindivs = data.getNumberOfIndividuals();
  if(numindivs > 1){
    Log.setDisplayMode(Quiet);
    //Log << numindivs << " individuals\n";
    if(numdiploidIndivs > 0){
      Log << numdiploidIndivs << " diploid "; 
      if(numdiploidIndivs < numindivs)Log<< "and ";
    }
    if(numdiploidIndivs < numindivs)Log << numindivs- numdiploidIndivs<< " haploid ";
    Log << "individuals\n\n";
  }
  
  A.PrintPrior(Log);
  
  if( options.getNumberOfBlockStates()>1 ) {
    L->SetHMMStateArrivalProbs((bool)(numdiploidIndivs !=0));
  }
  
  //initialise regression objects
  if (options.getNumberOfOutcomes()>0 ){
    InitialiseRegressionObjects(options, data, Log);
  }

  if(options.getMHTest())
    MHTest.Initialise(options.getPopulations(), &Loci, options.getResultsDir());
}

void HapMixModel::Iterate(const double* Coolnesses, unsigned coolness,
			  Options & options, InputData & data, LogWriter& Log, 
			  double & SumEnergy, double & SumEnergySq, 
			  bool AnnealedRun){

  double AISz = 0.0;
  if(!AnnealedRun) cout << endl;

  const unsigned NumberOfRuns = ((HapMixOptions&)options).GetNumStarts();
  const int samples = options.getSamples();
  const int burnin  = options.getBurnIn();

  for(unsigned run = 0; run < NumberOfRuns; ++run){
    if(run > 0){
      Log << Quiet << "\n** Restarting **\n";
      ReadInitialValuesFromFile(run, (HapMixOptions&)options, Log);
    }

    for( int iteration = 0; iteration <= samples; iteration++ ) {
      
      //Sample Parameters
      UpdateParameters(iteration, options, Log, data, Coolnesses, coolness, AnnealedRun, SumEnergy, SumEnergySq, AISz);
      
      //Output Parameters
      if(!AnnealedRun){    
	// output parameters every 'getSampleEvery()' iterations
	if(!(iteration % options.getSampleEvery()) )
	  OutputParameters(iteration, options, Log);
	
	if( iteration > burnin &&
	    // output every 'getSampleEvery() * 10' iterations after burnin
	    !( (iteration - burnin) % (options.getSampleEvery() * 10))
	    ){ 
	  //Output ergodic averages
	  OutputErgodicAverages( iteration - burnin, SumEnergy, SumEnergySq);
	  //Output Test results
	  OutputTests((HapMixOptions&)options, data, Log);	
	}   
      }
      
    }// end loop over iterations
    //use Annealed Importance Sampling to calculate marginal likelihood
    if(coolness>0) AISsumlogz += log(AISz /= (double)(samples-burnin));

  }

}
void HapMixModel::UpdateParameters(int iteration, const Options& _options, LogWriter&, 
                                   const InputData & data, const double* Coolnesses, unsigned coolness_index, bool anneal, 
				   double & SumEnergy, double & SumEnergySq, double& AISz){
  const double coolness = Coolnesses[coolness_index];
  //cast Options pointer to HapMixOptions for access to HAPMIXMAP options
  const HapMixOptions& options = (const HapMixOptions&) _options;
  
  A.ResetAlleleCounts(options.getPopulations());
  
  ///////////////////////////////////////////////////////////////
  // Update individual-level parameters, sampling hidden states
  ///////////////////////////////////////////////////////////////
    
    HMIC->SampleHiddenStates(options, iteration);
    //accumulate energy
    if( iteration > options.getBurnIn()) 
      AccumulateEnergy(Coolnesses, coolness_index, _options, SumEnergy, SumEnergySq, AISz, anneal, iteration );
  
    //Write Iteration Number to screen
    if( !anneal &&  !(iteration % options.getSampleEvery()) ) {
      WriteIterationNumber(iteration, (int)log10((double)options.getTotalSamples()+1 ), options.getDisplayLevel());
    }

  IC->AccumulateAlleleCounts(options, &A, &Loci, anneal); 
    
  if( options.getHWTestIndicator() || options.getMHTest() || options.getTestForResidualAllelicAssoc() ) {
    // loops over individuals to sample hap pairs, not skipping missing genotypes. 
    // Does not update counts since done already
    IC->SampleHapPairs(options, &A, &Loci, false, anneal, false); 
  }

    //accumulate conditional genotype probs for masked individuals at masked loci
  if(options.OutputCGProbs() && iteration > options.getBurnIn())
    HMIC->AccumulateConditionalGenotypeProbs(options, (const InputHapMixData &)data, Loci.GetNumberOfCompositeLoci());
 
  ////////////////////////////////////////////////////////////////
  //score tests
  ///////////////////////////////////////////////////////////////
  if(  !anneal && iteration > options.getBurnIn() ){
    //update score tests every SCORETEST_UPDATE_EVERY iterations after burn-in
    if( options.getTestForAllelicAssociation()&& !(iteration%SCORETEST_UPDATE_EVERY) ){
      AllelicAssocTest.Update(HMIC, R[0], Loci);
    }
    if(options.getTestForResidualAllelicAssoc()/*&& !(iteration%SCORETEST_UPDATE_EVERY)*/){
      ResidualAllelicAssocScoreTest.Reset();
      ResidualAllelicAssocScoreTest.Update(A.GetAlleleFreqs(), true);
    }
  }
  if(options.getMHTest() && !anneal && (iteration > options.getBurnIn()))
    MHTest.Update(IC, Loci);//update Mantel-Haenszel test

  //////////////////////////////////////////////////////////////////
  // update allele frequencies conditional on locus ancestry states
  ///////////////////////////////////////////////////////////////////
  if( !options.getFixedAlleleFreqs()){
    A.Update(IC, (iteration > options.getBurnIn() && !anneal), coolness);

    A.SetDiploidGenotypeProbs();
  }//end if random allele freqs
  
  //////////////////////////////////////////////////////////////////////////
  // update of allele freqs sets HMM probs and stored loglikelihoods as bad.
  // next update of stored loglikelihoods will be from getEnergy 
  //////////////////////////////////////////////////////////////////////////
  IC->HMMIsBad(true); 

  
  //////////////////////////////////////////////////////////////////////////////
  //Sample arrival rates with Hamiltonian Sampler, using sampled hidden states. 
  //NB: requires accumulating of StateArrivalCounts in IC
  //////////////////////////////////////////////////////////////////////////////
  

  L->SampleArrivalRate(HMIC->getConcordanceCounts(), (!anneal && iteration > options.getBurnIn() 
						      && options.getPopulations() > 1) );
  
  //sample mixture proportions with conjugate update
  if(!options.getFixedMixtureProps())
    L->SampleMixtureProportions(HMIC->getSumArrivalCounts());
  
  //Set global StateArrivalProbs in HMM objects. Do not force setting of mixture props (if fixed)
  L->SetHMMStateArrivalProbs( (numdiploidIndivs>0));

  ///////////////////////////////////////////////////////////////////////////
  // update regression parameters (if regression model)
  //////////////////////////////////////////////////////////////////////////

    bool condition = (!anneal && iteration > options.getBurnIn());
    for(unsigned r = 0; r < R.size(); ++r){
    R[r]->Update(condition, IC->getOutcome(r), coolness );
    //output expected values of outcome variables to file every 'every' iterations after burnin
    if(condition && !(iteration % options.getSampleEvery()) ) {
      R[r]->OutputExpectedY();
    }
  }
}

///write ergodic averages to file
void HapMixModel::OutputErgodicAverages(int samples, double & SumEnergy, double & SumEnergySq){
  if ( avgstream.is_open() ){
    L->OutputErgodicAvg(samples, avgstream);//average arrival rates
    A.OutputErgodicAvg(samples, &avgstream);//freq Dirichlet param prior rate
    
    for(unsigned r = 0; r < R.size(); ++r)//regression params
      R[r]->OutputErgodicAvg(samples, avgstream);
    
    OutputErgodicAvgDeviance(samples, SumEnergy, SumEnergySq);
    avgstream << endl;
  }
}

///Write score test output
void HapMixModel::OutputTests(HapMixOptions& options, InputData & data, LogWriter& Log  ){
  if(options.getTestForResidualAllelicAssoc())
    ResidualAllelicAssocScoreTest.Output(data.getLocusLabels());

  if( options.getTestForAllelicAssociation() )    {
    AllelicAssocTest.Output(Loci);
  }
  
  if(options.getMHTest() ){
    MHTest.Output(data.getLocusLabels());
  }
}

void HapMixModel::OutputParameters(int iteration, const Options& options, LogWriter& Log){

  Log.setDisplayMode(Quiet);

  bclib::Delimitedstdout ScreenWriter(' ');

  //output sample mean and variance of arrival rates and the prior parameters
  if(options.getPopulations() > 1){ 
    //write to Robject
    if( iteration > options.getBurnIn() ){
      //write to file
      L->OutputParams();

      //write to R object
      L->OutputParams(paramstream);
    }
    if( options.getDisplayLevel() > 2 )
      L->OutputParams(ScreenWriter);
  }

  if( !options.getFixedAlleleFreqs() ){
    //output Allele freq prior params to file if after burnin and to screen if displaylevel >2
    if(iteration > options.getBurnIn()){
      A.OutputPriorParams();
      A.OutputPriorParams(paramstream);
    }

    if(options.getDisplayLevel()>2)
      A.OutputPriorParams(ScreenWriter);
  }

  // ** regression parameters
  for(unsigned r = 0; r < R.size(); ++r){
    //output regression parameters to file if after burnin and to screen if displaylevel >2
    if( iteration > options.getBurnIn() ){
      R[r]->OutputParams(paramstream);
      R[r]->OutputParams(options.getNumberOfOutcomes());
    }
    if(options.getDisplayLevel()>2){
      if(R.size()>1)
	cout << "\nRegression " << r +1 << " ";
      R[r]->OutputParams(ScreenWriter);
     }
  }
  if( iteration > options.getBurnIn() )
    paramstream << bclib::newline;

  //if( options->getDisplayLevel()>2 ) cout << endl;
  // ** new line in log file but not on screen 
  if( iteration == 0 ) {
    Log << Off << "\n"  << Quiet;
  }
  // cout << endl;
}

void HapMixModel::PrintAcceptanceRates(const Options& options, LogWriter& Log,
                                       const Vector_s& /* PopulationLabels */) {

  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
  else Log.setDisplayMode(On);

  L->printAcceptanceRates(Log);

  A.OutputAlleleFreqSamplerAcceptanceRates((options.getResultsDir() + "/AlleleFreqSamplerAcceptanceRates.txt").c_str());
  if(!options.getFixedAlleleFreqs() ){
    Log << "Average expected Acceptance rate in allele frequency prior parameter sampler:\n" << A.getAcceptanceRate()
	<< "\nwith average final step size of " << A.getStepSize() << "\n";
  }
}

/// Final tasks to perform just before exit
void HapMixModel::Finalize(const Options& _options, LogWriter& Log, const InputData& data ){
  //cast Options object to HapMixOptions for access to HAPMIXMAP options
  const HapMixOptions& options = (const HapMixOptions&) _options;

  //Output results of Mantel-Haentszel Test
  if(options.getMHTest()){
    MHTest.WriteFinalTable(options.getResultsDir(), data.getLocusLabels(), Log);
  }
  //Write final score test tables  
  if(options.getTestForResidualAllelicAssoc())
    ResidualAllelicAssocScoreTest.WriteFinalTable((options.getResultsDir() + "/" + RESIDUAL_LD_TEST_FINAL).c_str(), data.getLocusLabels(), Log);
  
  if( options.getTestForAllelicAssociation() )    {
    AllelicAssocTest.WriteFinalTable(options.getResultsDir(), Loci, (InputHapMixData&)data, Log);
  }
  
  //output posterior means of lambda (expected number of arrivals)
  L->OutputArrivalRatePosteriorMeans(options.getArrivalRateOutputFilename(), options.getTotalSamples()-options.getBurnIn(), 
				     data.getUnitOfDistanceAsString());
  //output final values of arrival rates and their prior params
  L->OutputArrivalRates(options.getFinalLambdaFilename().c_str());
  //output final values of mixture proportions
  L->OutputMixtureProps(options.getFinalMixturePropsFilename().c_str());
  
  //output final values of allele freq prior
  A.OutputFinalValues(options.getFinalFreqPriorFilename().c_str(), Log);
  
  //output final values of allelefreqs
  A.OutputAlleleFreqs(options.getFinalAlleleFreqFilename().c_str(), Log);
  
  //output posterior means of allele freq precision
  if(options.OutputAlleleFreqPrior())
    A.OutputPosteriorMeans(options.getAlleleFreqPriorOutputFilename(), Log);
  
  
  //output posterior means of predictive genotype probs at masked loci for masked individuals
  if(options.OutputCGProbs()){
    std::string s = options.getResultsDir();
    s.append("/PPGenotypeProbs.txt");

    const Vector_s& LocusLabels = data.getLocusLabels();
    Vector_s MaskedLocusLabels;
    for(unsigned locus = 0; locus < Loci.GetNumberOfCompositeLoci(); ++locus)
      if(!((InputHapMixData&)data).isTypedLocus(locus))
       MaskedLocusLabels.push_back(LocusLabels[locus]);

    HMIC->OutputCGProbs(s.c_str(), MaskedLocusLabels);
  }
  //if(options.outputParams())
    {
      WriteParamsAsRObjectDimensions(options, data);
    }

}

void HapMixModel::WriteParamsAsRObjectDimensions(const HapMixOptions& options, const InputData& data){
  //write dimensions and dimnames of R object holding sampled parameters

  //arrival rate and mixture props summaries
  vector<vector<string> > dimnames(1);
  //if(strlen(options.getParameterFilename()))
  {
    dimnames[0].push_back("Arrivals.perMb.shapeParam");
    if(!L->fixedRateParameter())
      dimnames[0].push_back("Arrivals.perMb.rateParam");
    dimnames[0].push_back("Arrivals.perMb.Mean");
  }
  
  //residual allelic diversity summaries
  //if(strlen(options.getFreqPrecisionOutputFilename()))
  {
    dimnames[0].push_back("RAD.Mean");
    dimnames[0].push_back("RAD.Var");
    if(options.isFreqPrecisionHierModel())
      dimnames[0].push_back("RAD.Prior.Rate");
  }
  
  //regression parameters
  //if(strlen(options.getRegressionOutputFilename()))
  {
    for(unsigned r = 0; r < R.size(); ++r){
      dimnames[0].push_back("intercept");
      const vector<string>& CovariateLabels = data.getCovariateLabels();
      for(vector<string>::const_iterator i = CovariateLabels.begin(); i < CovariateLabels.end(); ++i)
	dimnames[0].push_back(*i);
      if(R[r]->getRegressionType()==Linear)dimnames[0].push_back("precision");
    }
  }
    
  vector<int> dims;
  dims.push_back(dimnames[0].size());
  dims.push_back((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery());

  paramstream.close(dims, dimnames);
}

void HapMixModel::InitialiseTests(Options& options, const InputData& data, LogWriter& Log){
  if( options.getTestForAllelicAssociation() ){
    AllelicAssocTest.Initialise(options.getResultsDir(), Loci.GetNumberOfCompositeLoci());
  }
  if(options.getTestForResidualAllelicAssoc())
    ResidualAllelicAssocScoreTest.Initialise((options.getResultsDir() + "/" + RESIDUAL_LD_TEST_PVALUES).c_str(), IC, &Loci);
  
  //  InitializeErgodicAvgFile(&options, Log, data.GetHiddenStateLabels(), data.getCovariateLabels());
}
//this function is here because three different objects have to write to avgstream
void HapMixModel::InitializeErgodicAvgFile(const Options&  _options, LogWriter &Log,  
                                           const Vector_s& , const Vector_s& CovariateLabels){
  Log.setDisplayMode(Quiet);
  const HapMixOptions& options = (const HapMixOptions& )_options;

  //Open ErgodicAverageFile  
  if( strlen( options.getErgodicAverageFilename() ) ) {
    avgstream.open( options.getErgodicAverageFilename(), ios::out );
    if( ! avgstream ){
      Log.setDisplayMode(On);
      Log << "ERROR: Couldn't open Ergodic Average file\n";
      //exit( 1 );
    } else {
      Log << "Writing ergodic averages of parameters to "
	  << options.getErgodicAverageFilename() << "\n\n";
    }
    
    // Header line of ergodicaveragefile
    if(options.getPopulations()>1){
      avgstream << "ArrivalRate.mean\t";
      
    }//end if hierarchical model

    //rate parameter of prior on frequency Dirichlet prior params
    if(!options.getFixedAlleleFreqs() && options.isFreqPrecisionHierModel()){
      avgstream << "FreqPrecisionPriorRate\t"; 
    }

    // Regression parameters
    if( options.getNumberOfOutcomes() > 0 ){
      for(int r = 0; r < options.getNumberOfOutcomes(); ++r){
	if(IC->getOutcomeType(r)!=CoxData)
          avgstream << "intercept\t";

	//write covariate labels to header
        copy(CovariateLabels.begin(), CovariateLabels.end(), ostream_iterator<string>(avgstream, "\t")); 

	if( IC->getOutcomeType(r)==Continuous )//linear regression
	  avgstream << "precision\t";
      }
    }

    avgstream << "MeanDeviance\tVarDeviance\t";
    avgstream << "\n";
  } else {
    Log << "Not writing ergodic averages to file\n";
  }
}

double HapMixModel::getDevianceAtPosteriorMean(const Options& options, LogWriter& Log){
  //return HMIC->getDevianceAtPosteriorMean(options, R, &Loci, Log, L->getGlobalMixtureProps(), L->getSumLogRho());
  return 0.0;
}

void HapMixModel::ReadInitialValuesFromFile(unsigned startnum, const HapMixOptions& options, LogWriter& Log){
  L->ReadInitialArrivalRatesFromFile(options.getInitialArrivalRateFilename(startnum).c_str(), Log);
  if(!options.getFixedMixtureProps())
    L->ReadInitialMixturePropsFromFile(options.getInitialMixturePropsFilename(startnum).c_str(), Log);

  if(!options.getFixedAlleleFreqs()){
    A.LoadInitialAlleleFreqs(options.getInitialAlleleFreqFilename(startnum).c_str(), Log);
    A.ReadInitialPriorParamsFromFile(options.getInitialFreqPriorFilename(startnum).c_str(), Log);
  }
}
