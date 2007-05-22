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
#include "EventLogger.hh"
#include "HapMixIndividual.h"

#define SCORETEST_UPDATE_EVERY 2
int numdiploidIndivs = 0;

// HapMixModel::HapMixModel(){

// }

HapMixModel::~HapMixModel(){
  delete L;
}

void HapMixModel::Initialise(HapMixOptions& options, InputData& data,  LogWriter& Log){
  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  InitialiseGenome(Loci, options, data, Log);
  //check allelefreq files, initialises allele freqs and finishes setting up Composite Loci
  A.Initialise(&options, &data, &Loci, Log); 
  pA = &A;//set pointer to AlleleFreqs object

  L = new PopHapMix(&options, &Loci);
  if(isMaster || isWorker)L->Initialise(data.getUnitOfDistanceAsString(), Log);

  //create HapMixIndividualCollection object
  //NB call after A Initialise, and before R Initialise
  HMIC = new HapMixIndividualCollection(&options, &data, &Loci, L->getGlobalMixtureProps());
  //set IC pointer in base class to the same address as HMIC
  IC = HMIC;

  if(options.CheckData())
    data.CheckForMonomorphicLoci(Log);

  //load Outcome and Covariate data into IC
  if(isMaster || isWorker) IC->LoadData(&options, &data, false);

  // set unannealed probs
  //if(isWorker)IC->setGenotypeProbs(&Loci, &A);


  if(isMaster || isWorker){
    //get number of diploid individuals by summing over individuals
    numdiploidIndivs = IC->getNumDiploidIndividuals();
  }
  if(isWorker){
    Loci.SetHMMDimensions(options.getNumberOfBlockStates(), (bool)(numdiploidIndivs !=0));
  }

#ifdef PARALLEL
  //broadcast number of diploid individuals
  MPI::COMM_WORLD.Bcast(&numdiploidIndivs, 1, MPI::INT, 0);
#endif

  if(isFreqSampler)
    A.setSampler(options.getThermoIndicator(), (!numdiploidIndivs), 
                 ( !strlen(options.getPriorAlleleFreqFilename()) && !strlen(options.getInitialAlleleFreqFilename()) ));

  //if there are any diploid individuals, allocate space for diploid genotype probs
  if(numdiploidIndivs && (isFreqSampler || isWorker)){
    A.AllocateDiploidGenotypeProbs();
    A.SetDiploidGenotypeProbs();
  }
  HapMixIndividual::SetGenotypeProbs(&Loci, A.getHaploidGenotypeProbs(), A.getDiploidGenotypeProbs());

  if(isMaster){
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
  }
  
  if(isFreqSampler)A.PrintPrior(Log);
  
  if( options.getNumberOfBlockStates()>1 && isWorker) {
    L->SetHMMStateArrivalProbs((bool)(numdiploidIndivs !=0));
  }
  
  //initialise regression objects
  if (options.getNumberOfOutcomes()>0 && (isMaster || isWorker)){
    InitialiseRegressionObjects(options, data, Log);
  }

  if(options.getMHTest())
    MHTest.Initialise(options.getPopulations(), &Loci, options.getMHTestFilename(), Log);
}

void HapMixModel::Iterate(const int & samples, const int & burnin, const double* Coolnesses, unsigned coolness,
			  Options & options, InputData & data, LogWriter& Log, 
			  double & SumEnergy, double & SumEnergySq, 
			  bool AnnealedRun){
  const bool isMaster = Comms::isMaster();

  double AISz = 0.0;
  if(isMaster && !AnnealedRun) cout << endl;

  for( int iteration = 0; iteration <= samples; iteration++ ) {
     //Write Iteration Number to screen
    if( isMaster && !AnnealedRun &&  !(iteration % options.getSampleEvery()) ) {
      WriteIterationNumber(iteration, (int)log10((double) samples+1 ), options.getDisplayLevel());
    }

    //Sample Parameters
    UpdateParameters(iteration, &options, Log, data.GetPopLabels(), Coolnesses, coolness, AnnealedRun, SumEnergy, SumEnergySq, AISz);
    SubIterate(iteration, burnin, options, data, Log, SumEnergy, SumEnergySq,
	       AnnealedRun);
	
  }// end loop over iterations
  //use Annealed Importance Sampling to calculate marginal likelihood
  if(coolness>0) AISsumlogz += log(AISz /= (double)(samples-burnin));

}
void HapMixModel::UpdateParameters(int iteration, const Options * _options, LogWriter&, 
                                   const Vector_s&, const double* Coolnesses, unsigned coolness_index, bool anneal, 
				   double & SumEnergy, double & SumEnergySq, double& AISz){
  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();
  const double coolness = Coolnesses[coolness_index];
  //cast Options pointer to HapMixOptions for access to HAPMIXMAP options
  const HapMixOptions* options = (const HapMixOptions*) _options;
  
  
  A.ResetAlleleCounts(options->getPopulations());

  ///////////////////////////////////////////////////////////////
  // Update individual-level parameters, sampling hidden states
  ///////////////////////////////////////////////////////////////
    
  if(isMaster || isWorker){
   HMIC->SampleHiddenStates(options, iteration);
  }
  //accumulate energy
  if( (isMaster || isWorker) && iteration > options->getBurnIn()) 
    GetEnergy(Coolnesses, coolness_index, *_options, SumEnergy, SumEnergySq, AISz, anneal, iteration );
  
  if(isWorker || isFreqSampler){
    IC->AccumulateAlleleCounts(options, &A, &Loci, anneal); 
    
    if( options->getHWTestIndicator() || options->getMHTest() || options->getTestForResidualAllelicAssoc() ) {
    EventLogger::LogEvent(13, iteration, "sampleHapPairs");
    // loops over individuals to sample hap pairs, not skipping missing genotypes. 
    // Does not update counts since done already
    IC->SampleHapPairs(options, &A, &Loci, false, anneal, false); 
    EventLogger::LogEvent(14, iteration, "sampledHapPairs");
    }
  }

  if(isWorker){
    //accumulate conditional genotype probs for masked individuals at masked loci
    if(options->OutputCGProbs() && iteration > options->getBurnIn())
      HMIC->AccumulateConditionalGenotypeProbs(options, Loci);
  }
  
  ////////////////////////////////////////////////////////////////
  //score tests
  ///////////////////////////////////////////////////////////////
  if( (isMaster || isWorker) && !anneal && iteration > options->getBurnIn() ){
    //update score tests every SCORETEST_UPDATE_EVERY iterations after burn-in
    if( options->getTestForAllelicAssociation()&& !(iteration%SCORETEST_UPDATE_EVERY) ){
      AllelicAssocTest.Update(HMIC, R[0], Loci);
    }
    if(options->getTestForResidualAllelicAssoc()/*&& !(iteration%SCORETEST_UPDATE_EVERY)*/){
      EventLogger::LogEvent(21, iteration, "ResLDTestStart");
      ResidualAllelicAssocScoreTest.Reset();
      ResidualAllelicAssocScoreTest.Update(A.GetAlleleFreqs(), true);
      EventLogger::LogEvent(22, iteration, "ResLDTestEnd");
    }
  }
  if(options->getMHTest() && !anneal && (iteration > options->getBurnIn()))
    MHTest.Update(IC, Loci);//update Mantel-Haenszel test

  //////////////////////////////////////////////////////////////////
  // update allele frequencies conditional on locus ancestry states
  ///////////////////////////////////////////////////////////////////
  if( !options->getFixedAlleleFreqs()){
    if(isFreqSampler){
      EventLogger::LogEvent(7, iteration, "SampleFreqs");
      A.Update(IC, (iteration > options->getBurnIn() && !anneal), coolness, (numdiploidIndivs==0));
      EventLogger::LogEvent(8, iteration, "SampledFreqs");
    }
#ifdef PARALLEL
    if(isWorker || isFreqSampler ) { 
      A.BroadcastAlleleFreqs();
    }
#endif
    A.SetDiploidGenotypeProbs();
  }//end if random allele freqs
  
  //////////////////////////////////////////////////////////////////////////
  // update of allele freqs sets HMM probs and stored loglikelihoods as bad.
  // next update of stored loglikelihoods will be from getEnergy 
  //////////////////////////////////////////////////////////////////////////
  if(isWorker ) {   
    IC->HMMIsBad(true); 
  } 
  
  //////////////////////////////////////////////////////////////////////////////
  //Sample arrival rates with Hamiltonian Sampler, using sampled hidden states. 
  //NB: requires accumulating of StateArrivalCounts in IC
  //////////////////////////////////////////////////////////////////////////////
  
  if(isMaster || isWorker){
    L->SampleArrivalRate(HMIC->getConcordanceCounts(), (!anneal && iteration > options->getBurnIn() 
							&& options->getPopulations() > 1) );
    
    //sample mixture proportions with conjugate update
    if(!options->getFixedMixtureProps())
      L->SampleMixtureProportions(HMIC->getSumArrivalCounts());
    
    //Set global StateArrivalProbs in HMM objects. Do not force setting of mixture props (if fixed)
    if( isWorker) {
      L->SetHMMStateArrivalProbs( (numdiploidIndivs>0));
    }
    
  }
  ///////////////////////////////////////////////////////////////////////////
  // update regression parameters (if regression model)
  //////////////////////////////////////////////////////////////////////////
  if(isMaster){
    bool condition = (!anneal && iteration > options->getBurnIn() && isMaster);
    for(unsigned r = 0; r < R.size(); ++r){
      R[r]->Update(condition, IC->getOutcome(r), coolness );
      //output expected values of outcome variables to file every 'every' iterations after burnin
      if(condition && !(iteration % options->getSampleEvery()) ) {
	R[r]->OutputExpectedY();
      }
    }
  }
#ifdef PARALLEL
  if(isMaster || isWorker){    //broadcast regression parameters to workers
    for(unsigned r = 0; r < R.size(); ++r){
      Comms::BroadcastRegressionParameters(R[r]->getbeta(), R[r]->getNumCovariates());
      //if(R[r]->getRegressionType()==Linear) Comms::BroadcastRegressionPrecision(R[r]->getlambda());
    }
  }
#endif
}

void HapMixModel::SubIterate(int iteration, const int & burnin, Options & _options, InputData & data, 
			     LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
			     bool AnnealedRun){
  const bool isMaster = Comms::isMaster();
  //const bool isFreqSampler = Comms::isFreqSampler();
  //  const bool isWorker = Comms::isWorker();
  //cast Options object to HapMixOptions for access to HAPMIXMAP options
  HapMixOptions& options = (HapMixOptions&) _options;

  //Output parameters to file and to screen
  Log.setDisplayMode(Quiet);
  if(!AnnealedRun){    
    // output every 'getSampleEvery()' iterations
    if(!(iteration % options.getSampleEvery()) && (isMaster/* || isFreqSampler*/))
      OutputParameters(iteration, &options, Log);
    
    //Updates and Output after BurnIn     
    if( !AnnealedRun && iteration > burnin && isMaster){
		
      // output every 'getSampleEvery() * 10' iterations (still after BurnIn)
      if (!( (iteration - burnin) % (options.getSampleEvery() * 10))){    
        //Ergodic averages
        Log.setDisplayMode(On);
        if ( strlen( options.getErgodicAverageFilename() ) ){
          int samples = iteration - burnin;

          L->OutputErgodicAvg(samples, avgstream);//average arrival rates
          A.OutputErgodicAvg(samples, &avgstream);//freq Dirichlet param prior rate

          for(unsigned r = 0; r < R.size(); ++r)//regression params
            R[r]->OutputErgodicAvg(samples,&avgstream);

          OutputErgodicAvgDeviance(samples, SumEnergy, SumEnergySq);
          avgstream << endl;
        }
        //Score Test output
	ResidualAllelicAssocScoreTest.Output(false, data.getLocusLabels());
	if( options.getTestForAllelicAssociation() )    {
	  string filename; //ignored if !final
	  AllelicAssocTest.Output(data.GetPopLabels(), Loci, false, filename.c_str());
	}

        if(options.getMHTest() && Comms::isMaster()){
          MHTest.Output(options.getMHTestFilename(), data.getLocusLabels(), false);
        }
      }//end "if every'*10" block
    }//end "if after BurnIn" block
  } // end "if not AnnealedRun" block
}

void HapMixModel::OutputParameters(int iteration, const Options *options, LogWriter& Log){

  // fix so that params can be output to console  
  Log.setDisplayMode(Quiet);

  //output sample mean and variance of arrival rates and the prior parameters
  if(Comms::isMaster() && options->getPopulations() > 1) L->OutputParams(iteration, Log);

  if( Comms::isFreqSampler() &&  !options->getFixedAlleleFreqs() )
    //output Allele freq prior params to file if after burnin and to screen if displaylevel >2
    A.OutputPriorParams((iteration > options->getBurnIn()), (bool)(options->getDisplayLevel()>2));

  // ** regression parameters
  if(Comms::isMaster())
    for(unsigned r = 0; r < R.size(); ++r)
      //output regression parameters to file if after burnin and to screen if displaylevel >2
      R[r]->Output(options->getNumberOfOutcomes(), (bool)(options->getDisplayLevel()>2), (bool)(iteration > options->getBurnIn()) );
  
  //if( options->getDisplayLevel()>2 ) cout << endl;
  // ** new line in log file but not on screen 
  if( iteration == 0 ) {
    Log << Off << "\n"  << Quiet;
  }
  // cout << endl;
}

void HapMixModel::PrintAcceptanceRates(const Options& options, LogWriter& Log){

  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
  else Log.setDisplayMode(On);

  if(Comms::isMaster()){
    L->printAcceptanceRates(Log);
  }
  A.OutputAlleleFreqSamplerAcceptanceRates((options.getResultsDir() + "/AlleleFreqSamplerAcceptanceRates.txt").c_str());
  if(!options.getFixedAlleleFreqs() && Comms::isFreqSampler()){
    Log << "Average expected Acceptance rate in allele frequency prior parameter sampler:\n" << A.getAcceptanceRate()
	<< "\nwith average final step size of " << A.getStepSize() << "\n";
  }
}

/// Final tasks to perform just before exit
void HapMixModel::Finalize(const Options& _options, LogWriter& Log, const InputData& data ){
  //cast Options object to HapMixOptions for access to HAPMIXMAP options
  const HapMixOptions& options = (const HapMixOptions&) _options;

  //Output results of Mantel-Haentszel Test
  if(options.getMHTest() && Comms::isMaster()){
    std::string s = options.getResultsDir();
    s.append("/MHTestFinal.txt");
    MHTest.Output(s.c_str(), data.getLocusLabels(), true);
  }
  //Write final score test tables  
  if(Comms::isMaster()){
    ResidualAllelicAssocScoreTest.ROutput();
    ResidualAllelicAssocScoreTest.Output(true, data.getLocusLabels());

    if( options.getTestForAllelicAssociation() )    {
      string filename; //ignored if !final
      filename = (options.getResultsDir());
      filename.append("/AllelicAssocTestsFinal.txt");
      AllelicAssocTest.Output(data.GetPopLabels(), Loci, true, filename.c_str());
      AllelicAssocTest.ROutput();
    }

    //output posterior means of lambda (expected number of arrivals)
    L->OutputArrivalRatePosteriorMeans(options.getArrivalRateOutputFilename(), options.getTotalSamples()-options.getBurnIn(), 
				       data.getUnitOfDistanceAsString());
    //output final values of arrival rates and their prior params
    L->OutputArrivalRates(options.getFinalLambdaFilename());
    //output final values of mixture proportions
    L->OutputMixtureProps(options.getFinalMixturePropsFilename());
  }
  if(Comms::isFreqSampler()){
    //output final values of allele freq prior
    A.OutputFinalValues(options.getFinalFreqPriorFilename(), Log);

    //output final values of allelefreqs
    A.OutputAlleleFreqs(options.getAlleleFreqOutputFilename(), Log);

    //output posterior means of allele freq precision
    if(options.OutputAlleleFreqPrior())
      A.OutputPosteriorMeans(options.getAlleleFreqPriorOutputFilename(), Log);
  }
  
  //output posterior means of predictive genotype probs at masked loci for masked individuals
  //TODO: parallelize - requires reducing sums in GenotypeProbOutputter class
  if(options.OutputCGProbs()){
    std::string s = options.getResultsDir();
    s.append("/PPGenotypeProbs.txt");
    const vector<unsigned>& maskedLoci = options.getMaskedLoci();
    const Vector_s& LocusLabels = data.getLocusLabels();
    Vector_s MaskedLocusLabels;
    for(vector<unsigned>::const_iterator l = maskedLoci.begin(); l!=maskedLoci.end(); ++l)
      // Masked loci indices are 1-based, need to offset by one
      MaskedLocusLabels.push_back(LocusLabels[(*l) - 1]);

    HMIC->OutputCGProbs(s.c_str(), MaskedLocusLabels);
  }
}
void HapMixModel::InitialiseTests(Options& options, const InputData& data, LogWriter& Log){
  const bool isMaster = Comms::isMaster();
  //const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();
  if(isMaster || isWorker){
    if( options.getTestForAllelicAssociation() ){
      AllelicAssocTest.Initialise(options.getAllelicAssociationScoreFilename(), 1, Loci.GetNumberOfCompositeLoci(), Log, false);
    }
    ResidualAllelicAssocScoreTest.Initialise(&options, IC, &Loci, Log);
  }
  if(isMaster)
    InitializeErgodicAvgFile(&options, Log, data.GetPopLabels(), data.getCovariateLabels());
}
//this function is here because three different objects have to write to avgstream
void HapMixModel::InitializeErgodicAvgFile(const Options* const _options, LogWriter &Log,  
                                           const Vector_s& , const Vector_s& CovariateLabels){
  Log.setDisplayMode(Quiet);
  const HapMixOptions* options = (const HapMixOptions* )_options;

  //Open ErgodicAverageFile  
  if( strlen( options->getErgodicAverageFilename() ) ) {
    avgstream.open( options->getErgodicAverageFilename(), ios::out );
    if( ! avgstream ){
      Log.setDisplayMode(On);
      Log << "ERROR: Couldn't open Ergodic Average file\n";
      //exit( 1 );
    } else {
      Log << "Writing ergodic averages of parameters to "
	  << options->getErgodicAverageFilename() << "\n\n";
    }
    
    // Header line of ergodicaveragefile
    if(options->getPopulations()>1){
      avgstream << "ArrivalRate.mean\t";
      
    }//end if hierarchical model

    //rate parameter of prior on frequency Dirichlet prior params
    if(!options->getFixedAlleleFreqs() && options->isFreqPrecisionHierModel()){
      avgstream << "FreqPrecisionPriorRate\t"; 
    }

    // Regression parameters
    if( options->getNumberOfOutcomes() > 0 ){
      for(int r = 0; r < options->getNumberOfOutcomes(); ++r){
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

double HapMixModel::getDevianceAtPosteriorMean(const Options* const options, LogWriter& Log){
  return HMIC->getDevianceAtPosteriorMean(options, R, &Loci, Log, L->getGlobalMixtureProps(), L->getSumLogRho(), Loci.GetNumberOfChromosomes(), &A);
}
