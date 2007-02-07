#include "AdmixMapModel.h"

// AdmixMapModel::AdmixMapModel(){

// }

AdmixMapModel::~AdmixMapModel(){
    delete L;
}

void AdmixMapModel::Initialise(AdmixOptions& options, InputData& data,  LogWriter& Log){
  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  InitialiseLoci(options, data, Log);

  A.Initialise(&options, &data, &Loci, Log, options.getChibIndicator()); //checks allelefreq files, initialises allele freqs and finishes setting up Composite Loci
  pA = &A;//set pointer to AlleleFreqs object
  
  AdmixedIndividuals = new AdmixIndividualCollection(&options, &data, &Loci);//NB call after A Initialise;//and before L and R Initialise
  IC = (IndividualCollection*) AdmixedIndividuals;
  if(isMaster || isWorker)IC->LoadData(&options, &data, (!options.getTestForAdmixtureAssociation() && options.getPopulations() > 1));    
  if(isWorker)AdmixedIndividuals->setGenotypeProbs(&Loci, &A); // sets unannealed probs
  if(isMaster){
    const int numdiploid = IC->getNumDiploidIndividuals();
    const int numindivs = data.getNumberOfIndividuals();
    if(numindivs > 1){
      Log.setDisplayMode(Quiet);
      //Log << numindivs << " individuals\n";
      if(numdiploid > 0){
	Log << numdiploid << " diploid "; 
	if(numdiploid < numindivs)Log<< "and ";
      }
      if(numdiploid < numindivs)Log << numindivs- numdiploid<< " haploid ";
      Log << "individuals\n\n";
    }
  }
  
  L = new PopAdmix(&options, &Loci);    
  if(isMaster || isWorker)L->Initialise(IC->getSize(), data.GetPopLabels(), Log);
  if(isFreqSampler)A.PrintPrior(data.GetPopLabels(), Log);
  
  if( options.getPopulations()>1 && (options.isGlobalRho()) && isWorker) {
    Loci.SetLocusCorrelation(L->getrho());
  }
  
  if(isMaster || isWorker){
      AdmixedIndividuals->Initialise(&options, &Loci, data.GetPopLabels(), Log);
      AdmixedIndividuals->DrawInitialAdmixture(L->getalpha());
  }
  
  //initialise regression objects
  if (options.getNumberOfOutcomes()>0 && (isMaster || isWorker)){
    InitialiseRegressionObjects(options, data, Log);
  }
}

void AdmixMapModel::UpdateParameters(int iteration, const Options *_options, LogWriter& Log, 
		      const Vector_s& PopulationLabels, double coolness, bool anneal){
  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  //cast Options pointer to AdmixOptions pointer for access to ADMIXMAP options
  //TODO: change pointer to reference
  const AdmixOptions* options = (const AdmixOptions*) _options;

  A.ResetAlleleCounts(options->getPopulations());

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if(isMaster || isWorker){
    // ** update global sumintensities conditional on genotype probs and individual admixture proportions
     if((options->getPopulations() > 1) && options->getIndAdmixHierIndicator() && 
        (Loci.GetLengthOfGenome() + Loci.GetLengthOfXchrm() > 0.0))
       L->UpdateGlobalSumIntensities(AdmixedIndividuals, (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1));
    // leaves individuals with HMM probs bad, stored likelihood ok
    // this function also sets locus correlations in Chromosomes
  }

  //find posterior modes of individual admixture at end of burn-in
  //set Chib numerator
  if(!anneal && iteration == options->getBurnIn() && (options->getChibIndicator() || strlen(options->getIndAdmixModeFilename()))) {
    AdmixedIndividuals->FindPosteriorModes(options, R, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), &A, PopulationLabels);
    if( options->getChibIndicator() ) {
      AdmixedIndividuals->setChibNumerator(options, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), &A);
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
     HMMUpdates does 4 things:
     (1)Sample individual admixture on even-numbered iterations
     (2)sample locus ancestry states
     (3)sample jump indicators, accumulating sufficient stats for updates of individual admixture and sumintensities (where required).
     (4)updates score, varscore and info in ancestry assoc tests, if required
     HMM forward updates are invoked as necessary. 
     This is when computing loglikelihood for random walk update of admixture (once each for current and proposal) then again if the proposal is rejected.
     When not using random-walk sampler (odd iterations), it is when sampling ancestry.
     Backward updates are invoked just before the score test update but only if required.
  */
  if(isMaster || isWorker){
    AdmixedIndividuals->HMMUpdates(iteration, options, R, L->getpoptheta(), L->getalpha(), Scoretests.getAffectedsOnlyTest(), Scoretests.getAncestryAssocTest(), anneal);
  }

  if(isWorker || isFreqSampler) {
#ifdef PARALLEL
    MPE_Log_event(13, iteration, "sampleHapPairs");
#endif
// loops over individuals to sample hap pairs then increment allele counts, skipping missing genotypes
    AdmixedIndividuals->SampleHapPairs(options, &A, &Loci, true, anneal, true);
#ifdef PARALLEL
    MPE_Log_event(14, iteration, "sampledHapPairs");
#endif
  }
  
#ifdef PARALLEL
  if(isWorker || isFreqSampler){
    A.SumAlleleCountsOverProcesses(options->getPopulations());
  }
#endif
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if( (isMaster || isWorker) && !anneal && iteration > options->getBurnIn() ){
    //score tests
    if( options->getScoreTestIndicator() ){
#ifdef PARALLEL
      MPE_Log_event(21, iteration, "ScoreTestUpdatestart");
#endif
      //TODO: broadcast dispersion in linear regression
      Scoretests.Update(R);//score tests evaluated for first outcome var only
#ifdef PARALLEL
      MPE_Log_event(22, iteration, "ScoreTestUpdateEnd");
#endif
    }
    if(options->getTestForResidualAllelicAssoc()){
#ifdef PARALLEL
      MPE_Log_event(21, iteration, "ResLDTestStart");
#endif
      Scoretests.UpdateScoresForResidualAllelicAssociation(A.GetAlleleFreqs());
#ifdef PARALLEL
      MPE_Log_event(22, iteration, "ResLDTestEnd");
#endif
    }
  }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // sample individual admixture and sum-intensities 
  AdmixedIndividuals->SampleParameters(iteration, options, R, L->getpoptheta(), L->getalpha(),
		       L->getrhoalpha(), L->getrhobeta(), Scoretests.getAncestryAssocTest(),anneal);
  // stored HMM likelihoods will now be bad if the sum-intensities are set at individual level
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if( options->getChibIndicator() && !anneal && iteration >= options->getBurnIn() )
    AdmixedIndividuals->updateChib(options, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), &A);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // update allele frequencies conditional on locus ancestry states
  // TODO: this requires fixing to anneal allele freqs for historicallelefreq model
  if( !options->getFixedAlleleFreqs()){
    if(isFreqSampler){
#ifdef PARALLEL
      MPE_Log_event(7, iteration, "SampleFreqs");
#endif
      A.Update(IC, (iteration > options->getBurnIn() && !anneal), coolness);
#ifdef PARALLEL
    MPE_Log_event(8, iteration, "SampledFreqs");
#endif
    }
#ifdef PARALLEL
    if(isWorker || isFreqSampler ) { 
      A.BroadcastAlleleFreqs();
    }
#endif
  }//end if random allele freqs
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // even for fixed allele freqs, must reset annealed genotype probs as unnannealed

  if(isWorker && (!options->getFixedAlleleFreqs() || anneal)) {   
#ifdef PARALLEL
    MPE_Log_event(11, iteration, "setGenotypeProbs"); 
#endif
    AdmixedIndividuals->setGenotypeProbs(&Loci, &A); // sets unannealed probs ready for getEnergy
    AdmixedIndividuals->HMMIsBad(true); // update of allele freqs sets HMM probs and stored loglikelihoods as bad
#ifdef PARALLEL
    MPE_Log_event(12, iteration, "GenotypeProbsSet"); 
#endif
  } // update of allele freqs sets HMM probs and stored loglikelihoods as bad
  
  // next update of stored loglikelihoods will be from getEnergy if not annealing run, from updateRhowithRW if globalrho, 
  // or from update of individual-level parameters otherwise
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(isMaster || isWorker){
    //update population admixture Dirichlet parameters conditional on individual admixture
    L->UpdatePopAdmixParams(iteration, AdmixedIndividuals, Log);
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // ** update regression parameters (if regression model) conditional on individual admixture
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
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
void AdmixMapModel::ResetStepSizeApproximators(int resetk){
  Model::ResetStepSizeApproximators(resetk);
 L->resetStepSizeApproximator(resetk);
}

void AdmixMapModel::SubIterate(int iteration, const int & burnin, Options& _options, InputData & data, 
			       LogWriter& Log, double& SumEnergy, double& SumEnergySq, 
			       bool AnnealedRun){
  const bool isMaster = Comms::isMaster();
  //const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();
  //cast Options object to AdmixOptions for access to ADMIXMAP options
  AdmixOptions& options = (AdmixOptions&) _options;

    //Output parameters to file and to screen
    Log.setDisplayMode(Quiet);
     if(!AnnealedRun){    
      // output every 'getSampleEvery()' iterations
      if(!(iteration % options.getSampleEvery()) && isMaster)
	OutputParameters(iteration, &options, Log);
      
      // ** set merged haplotypes for allelic association score test 
      if( (isMaster || isWorker) && iteration == options.getBurnIn() ){
#ifndef PARALLEL
	if(options.getTestForAllelicAssociation())
	  Scoretests.SetAllelicAssociationTest(L->getalpha0());
#endif
	if( isMaster && options.getStratificationTest() )
	  StratTest.Initialize( &options, Loci, IC, Log);
      }
	    
      //Updates and Output after BurnIn     
      if( !AnnealedRun && iteration > burnin && isMaster){
	//dispersion test
	if( options.getTestForDispersion() )DispTest.TestForDivergentAlleleFrequencies(&A, IC);
	//stratification test
	if( options.getStratificationTest() )StratTest.calculate(IC, A.GetAlleleFreqs(), Loci.GetChrmAndLocus(), 
								 options.getPopulations());
	//tests for mis-specified allelefreqs
	if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
	  AlleleFreqTest.Update(IC, &A, &Loci);
	//test for Hardy-Weinberg eq
	if( options.getHWTestIndicator() )
	  HWtest.Update(IC, &Loci);
		
	// output every 'getSampleEvery() * 10' iterations (still after BurnIn)
	if (!( (iteration - burnin) % (options.getSampleEvery() * 10))){    
	  //Ergodic averages
	  Log.setDisplayMode(On);
	  if ( strlen( options.getErgodicAverageFilename() ) ){
	    int samples = iteration - burnin;
	    if( options.getIndAdmixHierIndicator() ){
	      L->OutputErgodicAvg(samples,&avgstream);//pop admixture params, pop (mean) sumintensities
	      A.OutputErgodicAvg(samples, &avgstream);//dispersion parameter in dispersion model, freq Dirichlet param prior rate in hapmixmodel
	    }
	    for(unsigned r = 0; r < R.size(); ++r)//regression params
	      R[r]->OutputErgodicAvg(samples,&avgstream);

	    OutputErgodicAvgDeviance(samples, SumEnergy, SumEnergySq);
	    if(options.getChibIndicator()) AdmixedIndividuals->OutputErgodicChib(&avgstream, options.getFixedAlleleFreqs());
	    avgstream << endl;
	  }
	  //Score Test output
	  if( options.getScoreTestIndicator() )  Scoretests.Output(iteration-burnin, data.GetPopLabels(), data.getLocusLabels(), false);
	}//end "if every'*10" block
      }//end "if after BurnIn" block
    } // end "if not AnnealedRun" block
}


void AdmixMapModel::OutputParameters(int iteration, const AdmixOptions *options, LogWriter& Log){
  // fix so that params can be output to console  
  Log.setDisplayMode(Quiet);
  if(options->getIndAdmixHierIndicator()  ){
    //output population-level parameters only when there is a hierarchical model on indadmixture
    // ** pop admixture, sumintensities
    if(options->getPopulations() > 1) L->OutputParams(iteration, Log);
    // ** dispersion parameter (if dispersion model)
    A.OutputEta(iteration, options, Log);
  }
 
  // ** regression parameters
  for(unsigned r = 0; r < R.size(); ++r)
    R[r]->Output(options->getNumberOfOutcomes(), (bool)(options->getDisplayLevel()>2), (bool)(iteration > options->getBurnIn()) );
  
  // ** new line in log file but not on screen 
  if( iteration == 0 && (options->getIndAdmixHierIndicator() || options->getNumberOfOutcomes())) {
    Log.setDisplayMode(Off);
    Log << "\n";
    Log.setDisplayMode(Quiet);
  }
  
  //if( options->getDisplayLevel()>2 ) cout << endl;
  if( iteration > options->getBurnIn() ){
    // output individual and locus parameters every 'getSampleEvery()' iterations after burnin
    if( strlen( options->getIndAdmixtureFilename() ) ) AdmixedIndividuals->OutputIndAdmixture();
    if(options->getOutputAlleleFreq()){
	A.OutputAlleleFreqs();
    }
  }
  // cout << endl;
}

void AdmixMapModel::PrintAcceptanceRates(const Options& _options, LogWriter& Log){
  const AdmixOptions& options = (const AdmixOptions&)_options;
 if( options.getIndAdmixHierIndicator() ){
   
   if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
   else Log.setDisplayMode(On);
   
   if(options.getPopulations() > 1 && Comms::isMaster()){
     L->printAcceptanceRates(Log);
   }
   if(options.getCorrelatedAlleleFreqs()){
     Log<< "Expected acceptance rates in sampler for allele frequency proportion parameters: \n";
     for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci(); i++){
       if(Loci(i)->GetNumberOfStates()>2)
         Log << A.getAlphaSamplerAcceptanceRate(i) << " ";
     }
     Log<< "Expected acceptance rate in sampler for allele frequency dispersion parameter: \n";
     Log << A.getEtaSamplerAcceptanceRate(0)
         << "\nwith final step sizes of \n";
     for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci(); i++){
       if(Loci(i)->GetNumberOfStates()>2)
         Log << A.getAlphaSamplerStepsize(i) << " " ;
     }
     Log <<  A.getEtaSamplerStepsize(0) << "\n" ;
   }
   
   if( strlen( options.getHistoricalAlleleFreqFilename() )){
     Log << "Expected acceptance rates in allele frequency dispersion parameter samplers:\n ";
     for(int k = 0; k < options.getPopulations(); ++k){Log << A.getEtaSamplerAcceptanceRate(k)<< " " ;}
     Log << "\nwith final step sizes of ";
     for(int k = 0; k < options.getPopulations(); ++k){Log <<  A.getEtaSamplerStepsize(k) << " ";}
     Log << "\n";
   }
   A.OutputAlleleFreqSamplerAcceptanceRates((options.getResultsDir() + "/AlleleFreqSamplerAcceptanceRates.txt").c_str());
 }
}

void AdmixMapModel::Finalize(const Options& _options, LogWriter& Log, const InputData& data){
  const AdmixOptions& options = (const AdmixOptions&)_options;
  if( options.getChibIndicator()) {
    //AdmixedIndividuals->OutputChibEstimates(options.isRandomMatingAdmixMapModel(), Log, options.getPopulations());
    //MLEs of admixture & sumintensities used in Chib algorithm to estimate marginal likelihood
    if(AdmixedIndividuals->getSize()==1) AdmixedIndividuals->OutputChibResults(Log);
  }
  //FST
  if( strlen( options.getHistoricalAlleleFreqFilename()) && Comms::isFreqSampler()  ){
    A.OutputFST();
  }
  if(Comms::isFreqSampler())
    A.CloseOutputFile((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery(), data.GetPopLabels());

  //stratification test
  if( options.getStratificationTest() ) StratTest.Output(Log);
  //dispersion test
  if( options.getTestForDispersion() )  DispTest.Output(options.getTotalSamples() - options.getBurnIn(), Loci, data.GetPopLabels());
  //tests for mis-specified allele frequencies
  if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
    AlleleFreqTest.Output(options.getTotalSamples() - options.getBurnIn(), &Loci, data.GetPopLabels()); 
  //test for H-W eq
  if( options.getHWTestIndicator() )
    HWtest.Output(data.getLocusLabels()); 

  if( options.getScoreTestIndicator() && Comms::isMaster() ) {
    //finish writing score test output as R objects
    Scoretests.ROutput();
    //write final tables
    Scoretests.Output(options.getTotalSamples() - options.getBurnIn(), data.GetPopLabels(), data.getLocusLabels(), true);
  }
  //output to likelihood ratio file
  if(options.getTestForAffectedsOnly())
    //Individual::OutputLikRatios(options.getLikRatioFilename(), options.getTotalSamples()-options.getBurnIn(), data.GetPopLabels());
    Scoretests.OutputLikelihoodRatios(options.getLikRatioFilename(), options.getTotalSamples()-options.getBurnIn(), 
				      data.GetPopLabels());	
}
void AdmixMapModel::InitialiseTests(Options& _options, const InputData& data, LogWriter& Log){
  const bool isMaster = Comms::isMaster();
  //const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  //cast Options object to AdmixOptions for access to ADMIXMAP options
  AdmixOptions& options = (AdmixOptions&) _options;
  if( options.getScoreTestIndicator() && (isMaster || isWorker)){
    Scoretests.Initialise(&options, IC, &Loci, data.GetPopLabels(), Log);
  }

  //if(isMaster){
  if( options.getTestForDispersion() ){
    DispTest.Initialise(&options, Log, Loci.GetNumberOfCompositeLoci());    
  }
  if(options.getStratificationTest())
      StratTest.OpenOutputFile(options.getStratTestFilename(), Log);
  if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
    AlleleFreqTest.Initialise(&options, &Loci, Log );  
  if( options.getHWTestIndicator() )
    HWtest.Initialise(&options, Loci.GetTotalNumberOfLoci(), Log);
  if(isMaster)
    InitializeErgodicAvgFile(&options, Log, data.GetPopLabels(), data.getCovariateLabels());
  //}

}
//this function is here because three different objects have to write to avgstream
void AdmixMapModel::InitializeErgodicAvgFile(const AdmixOptions* const options, LogWriter &Log,  
				     const Vector_s& PopulationLabels, const Vector_s& CovariateLabels){
  Log.setDisplayMode(Quiet);
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
    if( options->getIndAdmixHierIndicator() ){
      if(options->getPopulations()>1){
	  for( int i = 0; i < options->getPopulations(); i++ ){
	    avgstream << PopulationLabels[i] << "\t";
	  }
	if( options->isGlobalRho() ) avgstream << "sumIntensities\t";
	else avgstream << "sumIntensities.mean\t";
      }
      
      // dispersion parameters
      if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
	for( int k = 0; k < options->getPopulations(); k++ ){
	  avgstream << "eta" << k << "\t";
	}
      }
    }//end if hierarchical model

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
    if(options->getChibIndicator()){// chib calculation
      avgstream << "LogPrior\tLogPosterior\tLogPosteriorAdmixture\tLogPosteriorSumIntensities\t"
		 << "LogPosteriorAlleleFreqs\tLogMarginalLikelihood";
    }
    avgstream << "\n";
  } else {
    Log << "Not writing ergodic averages to file\n";
  }
}


double AdmixMapModel::getDevianceAtPosteriorMean(const Options* const options, LogWriter& Log){
  return AdmixedIndividuals->getDevianceAtPosteriorMean(options, R, &Loci, Log, L->getSumLogRho(), Loci.GetNumberOfChromosomes(), &A);
}
double* AdmixMapModel::getSumEnergy()const{
    return AdmixedIndividuals->getSumEnergy();
}
double* AdmixMapModel::getSumEnergySq()const{
    return AdmixedIndividuals->getSumEnergySq();
}
