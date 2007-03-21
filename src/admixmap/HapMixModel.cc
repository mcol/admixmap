#include "HapMixModel.h"
//#include "MixturePropsWrapper.hh"
#include "EventLogger.hh"
// HapMixModel::HapMixModel(){

// }

HapMixModel::~HapMixModel(){
  delete L;
}

void HapMixModel::Initialise(HapMixOptions& options, InputData& data,  LogWriter& Log){
  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  InitialiseLoci(options, data, Log);
  A.Initialise(&options, &data, &Loci, Log); //checks allelefreq files, initialises allele freqs and finishes setting up Composite Loci
  pA = &A;//set pointer to AlleleFreqs object

  L = new PopHapMix(&options, &Loci);
  if(isMaster || isWorker)L->Initialise(data.getUnitOfDistanceAsString(), Log);

  //create HapMixIndividualCollection object
  //NB call after A Initialise, and before R Initialise
  HMIC = new HapMixIndividualCollection(&options, &data, &Loci, &A, L->getGlobalTheta());
  //set IC pointer in base class to the same address as HMIC
  IC = HMIC;

  //load Outcome and Covariate data into IC
  if(isMaster || isWorker) IC->LoadData(&options, &data, false);

  // set unannealed probs
  //if(isWorker)IC->setGenotypeProbs(&Loci, &A);

  int numdiploid = 0;
  if(isMaster || isWorker){
    //get number of diploid individuals by summing over individuals
    numdiploid = IC->getNumDiploidIndividuals();
  }

#ifdef PARALLEL
  //broadcast number of diploid individuals
  MPI::COMM_WORLD.Bcast(&numdiploid, 1, MPI::INT, 0);
#endif

  if(isFreqSampler)
    A.setSampler(options.getThermoIndicator(), (!numdiploid), 
                 ( !strlen(options.getPriorAlleleFreqFilename()) && !strlen(options.getInitialAlleleFreqFilename()) ));

  //if there are any diploid individuals, allocate space for diploid genotype probs
  if(numdiploid && (isFreqSampler || isWorker)){
    A.AllocateDiploidGenotypeProbs();
    A.SetDiploidGenotypeProbs();
  }

  if(isMaster){
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
  
  if(isFreqSampler)A.PrintPrior(Log);
  
  if( options.getNumberOfBlockStates()>1 && isWorker) {
    L->SetHMMStateArrivalProbs(true);
  }
  
  //initialise regression objects
  if (options.getNumberOfOutcomes()>0 && (isMaster || isWorker)){
    InitialiseRegressionObjects(options, data, Log);
  }

  if(options.getMHTest())
    MHTest.Initialise(options.getPopulations(), &Loci, options.getMHTestFilename(), Log);
}

void HapMixModel::UpdateParameters(int iteration, const Options * _options, LogWriter&, 
                                   const Vector_s&, const double* , double coolness, bool anneal){
  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();
  //cast Options pointer to HapMixOptions for access to HAPMIXMAP options
  const HapMixOptions* options = (const HapMixOptions*) _options;


  A.ResetAlleleCounts(options->getPopulations());

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Update individual-level parameters, sampling hidden states
  if(isMaster || isWorker){
    HMIC->SampleHiddenStates(options, iteration);
  }

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

  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if( (isMaster || isWorker) && !anneal && iteration > options->getBurnIn() ){
    //score tests
    if( options->getScoreTestIndicator() ){
      EventLogger::LogEvent(21, iteration, "ScoreTestUpdatestart");
      //TODO: broadcast dispersion in linear regression
      Scoretests.Update(R);//score tests evaluated for first outcome var only
      EventLogger::LogEvent(22, iteration, "ScoreTestUpdateEnd");
    }
    if(options->getTestForResidualAllelicAssoc()){
      EventLogger::LogEvent(21, iteration, "ResLDTestStart");
      Scoretests.UpdateScoresForResidualAllelicAssociation(A.GetAlleleFreqs());
      EventLogger::LogEvent(22, iteration, "ResLDTestEnd");
    }
  }
  if(options->getMHTest() && !anneal && (iteration > options->getBurnIn()))
    MHTest.Update(IC, Loci);//update Mantel-Haenszel test
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  // update allele frequencies conditional on locus ancestry states
  // TODO: this requires fixing to anneal allele freqs for historicallelefreq model
  if( !options->getFixedAlleleFreqs()){
    if(isFreqSampler){
      EventLogger::LogEvent(7, iteration, "SampleFreqs");
      A.Update(IC, (iteration > options->getBurnIn() && !anneal), coolness);
      EventLogger::LogEvent(8, iteration, "SampledFreqs");
    }
#ifdef PARALLEL
    if(isWorker || isFreqSampler ) { 
      A.BroadcastAlleleFreqs();
    }
#endif
    A.SetDiploidGenotypeProbs();
  }//end if random allele freqs
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(isWorker ) {   
    IC->HMMIsBad(true); 
  } // update of allele freqs sets HMM probs and stored loglikelihoods as bad
  // next update of stored loglikelihoods will be from getEnergy 
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(isMaster || isWorker){
    //L->UpdateGlobalTheta(iteration, IC);

    //Sample arrival rates with Hamiltonian Sampler, using sampled hidden states. 
    //NB: requires accumulating of StateArrivalCounts in IC
    L->SampleArrivalRate(HMIC->getConcordanceCounts(), (!anneal && iteration > options->getBurnIn() 
						       && options->getPopulations() > 1) );

    //sample mixture proportions with conjugate update
    if(!options->getFixedMixtureProps())
      L->SampleMixtureProportions(HMIC->getSumArrivalCounts());

    //Set global StateArrivalProbs in HMM objects. Do not force setting of mixture props (if fixed)
    if( isWorker) {
      L->SetHMMStateArrivalProbs(false);
    }

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
        if( options.getScoreTestIndicator() )  Scoretests.Output(iteration-burnin, data.GetPopLabels(), data.getLocusLabels(), false);
        if(options.getMHTest() && Comms::isMaster()){
          MHTest.Output(options.getMHTestFilename(), options.getTotalSamples() - options.getBurnIn(), data.getLocusLabels(), false);
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
    MHTest.Output(s.c_str(), options.getTotalSamples() - options.getBurnIn(), data.getLocusLabels(), true);
  }
  //Write final score test tables  
  if(Comms::isMaster()){
    if( options.getScoreTestIndicator() ) {
      //finish writing score test output as R objects
      Scoretests.ROutput();
      //write final tables
      Scoretests.Output(options.getTotalSamples() - options.getBurnIn(), data.GetPopLabels(), data.getLocusLabels(), true);
    }
    //output posterior means of lambda (expected number of arrivals)
    L->OutputLambdaPosteriorMeans(options.getArrivalRateOutputFilename(), options.getTotalSamples()-options.getBurnIn());
    //output final values of lambda and its prior params
    L->OutputLambda(options.getFinalLambdaFilename());
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
    HMIC->OutputCGProbs(s.c_str());
  }
}
void HapMixModel::InitialiseTests(Options& options, const InputData& data, LogWriter& Log){
  const bool isMaster = Comms::isMaster();
  //const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  if( options.getScoreTestIndicator() && (isMaster || isWorker)){
    Scoretests.Initialise(&options, IC, &Loci, data.GetPopLabels(), Log);
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
  return IC->getDevianceAtPosteriorMean(options, R, &Loci, Log, L->getSumLogRho(), Loci.GetNumberOfChromosomes(), &A);
}
