#include "HapMixModel.h"

// HapMixModel::HapMixModel(){

// }

HapMixModel::~HapMixModel(){
    delete L;
}

void HapMixModel::Initialise(AdmixOptions& options, InputData& data,  LogWriter& Log){
  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  InitialiseLoci(options, data, Log);
  A.Initialise(&options, &data, &Loci, Log); //checks allelefreq files, initialises allele freqs and finishes setting up Composite Loci
  pA = &A;//set pointer to AlleleFreqs object

  IC = new HapMixIndividualCollection(&options, &data, &Loci);//NB call after A Initialise;
  if(isMaster || isWorker) IC->LoadData(&options, &data, false);    //and before L and R Initialise
  if(isWorker)IC->setGenotypeProbs(&Loci, &A); // sets unannealed probs
  if(isMaster || isWorker){
    const int numdiploid = IC->getNumDiploidIndividuals();
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
  }

  L = new PopHapMix(&options, &Loci);
  if(isMaster || isWorker)L->Initialise(data.getUnitOfDistanceAsString(), Log);
  if(isFreqSampler)A.PrintPrior(Log);
  
  if( options.getPopulations()>1 && isWorker) {
    Loci.SetLocusCorrelation(L->getlambda());
    for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ ){
      //set global state arrival probs in hapmixmodel
      //TODO: can skip this if xonly analysis with no females
      //NB: assumes always diploid in hapmixmodel
      Loci.getChromosome(j)->SetHMMTheta(L->getGlobalTheta(), options.isRandomMatingModel(), true);
      Loci.getChromosome(j)->SetStateArrivalProbs(options.isRandomMatingModel(), true);
    }
  }
  
  //initialise regression objects
  if (options.getNumberOfOutcomes()>0 && (isMaster || isWorker)){
    InitialiseRegressionObjects(options, data, Log);
  }

  if(options.getMHTest())
    MHTest.Initialise(options.getPopulations(), &Loci, options.getMHTestFilename(), Log);
}

void HapMixModel::UpdateParameters(int iteration, const AdmixOptions *options, LogWriter&, 
		      const Vector_s&, double coolness, bool anneal){
  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  A.ResetAlleleCounts(options->getPopulations());

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//   if(iteration%2){
//     if(isMaster || isWorker){
//       L->UpdateSumIntensitiesByRandomWalk(IC,  
// 			      (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1) );
//     }
//   }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Update individual-level parameters, sampling hidden states
  if(isMaster || isWorker){
    IC->SampleLocusAncestry(options);
   }

  if(isWorker || isFreqSampler){
    IC->AccumulateAlleleCounts(options, &A, &Loci, anneal); 

    if( options->getHWTestIndicator() || options->getMHTest() || options->getTestForResidualAllelicAssoc() ) {
#ifdef PARALLEL
      MPE_Log_event(13, iteration, "sampleHapPairs");
#endif
      // loops over individuals to sample hap pairs, not skipping missing genotypes. Does not update counts since done already
      IC->SampleHapPairs(options, &A, &Loci, false, anneal, false); 
#ifdef PARALLEL
      MPE_Log_event(14, iteration, "sampledHapPairs");
#endif
    }
  }

  if(isWorker){
    //accumulate conditional genotype probs for masked individuals at masked loci
    if(options->OutputCGProbs() && iteration > options->getBurnIn())
      IC->AccumulateConditionalGenotypeProbs(options, Loci);
  }

  
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
  if(options->getMHTest() && !anneal && (iteration > options->getBurnIn()))
    MHTest.Update(IC, Loci);//update Mantel-Haenszel test
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
    IC->setGenotypeProbs(&Loci, &A); // sets unannealed probs ready for getEnergy
    IC->HMMIsBad(true); // update of allele freqs sets HMM probs and stored loglikelihoods as bad
#ifdef PARALLEL
    MPE_Log_event(12, iteration, "GenotypeProbsSet"); 
#endif
  } // update of allele freqs sets HMM probs and stored loglikelihoods as bad
  
  // next update of stored loglikelihoods will be from getEnergy if not annealing run, from updateRhowithRW if globalrho, 
  // or from update of individual-level parameters otherwise
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(isMaster || isWorker){
      //L->UpdateGlobalTheta(iteration, IC);
      //Hamiltonian Sampler, using sampled ancestry states. NB: requires accumulating of SumAncestry in IC
      L->SampleHapMixLambda(IC->getSumAncestry(),  
			    (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1) );
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

void HapMixModel::SubIterate(int iteration, const int & burnin, AdmixOptions & options, InputData & data, 
			     LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
			     bool AnnealedRun){
  const bool isMaster = Comms::isMaster();
  //const bool isFreqSampler = Comms::isFreqSampler();
//  const bool isWorker = Comms::isWorker();


    //Output parameters to file and to screen
    Log.setDisplayMode(Quiet);
     if(!AnnealedRun){    
      // output every 'getSampleEvery()' iterations
      if(!(iteration % options.getSampleEvery()) && isMaster)
	OutputParameters(iteration, &options, Log);
      
      //Updates and Output after BurnIn     
      if( !AnnealedRun && iteration > burnin && isMaster){
		
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

void HapMixModel::OutputParameters(int iteration, const AdmixOptions *options, LogWriter& Log){
  // fix so that params can be output to console  
  Log.setDisplayMode(Quiet);
  if(options->getIndAdmixHierIndicator()  ){
    //output population-level parameters only when there is a hierarchical model on indadmixture
    // ** pop admixture, sumintensities
    if(options->getPopulations() > 1) L->OutputParams(iteration, Log);
  }
  if((options->getDisplayLevel() > 2) && !options->getFixedAlleleFreqs())A.OutputPriorParams(cout, false);
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
      A.OutputPriorParams();

    if(options->getOutputAlleleFreq() && !options->getHapMixModelIndicator()){
	A.OutputAlleleFreqs();
    }
  }
  // cout << endl;
}

void HapMixModel::PrintAcceptanceRates(const AdmixOptions& options, LogWriter& Log){
  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
  else Log.setDisplayMode(On);

  if(options.getPopulations() > 1 && Comms::isMaster()){
    L->printAcceptanceRates(Log);
  }
  A.OutputAlleleFreqSamplerAcceptanceRates((options.getResultsDir() + "/AlleleFreqSamplerAcceptanceRates.txt").c_str());
  if(!options.getFixedAlleleFreqs() && Comms::isFreqSampler()){
    Log << "Average expected Acceptance rate in allele frequency prior parameter sampler:\n" << A.getAcceptanceRate()
	<< "\nwith average final step size of " << A.getStepSize() << "\n";
  }
}

/// Final tasks to perform just before exit
void HapMixModel::Finalize(const AdmixOptions& options, LogWriter& Log, const InputData& data ){
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
    std::string s = options.getResultsDir();
    s.append("/lambdaPosteriorMeans.txt");
    L->OutputLambdaPosteriorMeans(s.c_str(), options.getTotalSamples()-options.getBurnIn());
    //output final values of lambda
    const char* ss = options.getHapMixLambdaOutputFilename();
    if(strlen(ss))
      L->OutputLambda(ss);
  }
  if(Comms::isFreqSampler()){
    //output final values of allelefreqs
    const char* ss = options.getAlleleFreqOutputFilename();
    if(strlen(ss))
      A.OutputAlleleFreqs(ss, Log);

    //output posterior means of allele freq dispersion
    if(options.OutputAlleleFreqPrior())
      A.OutputPosteriorMeans(options.getAlleleFreqPriorOutputFilename(), Log);
  }
  
  //output posterior means of predictive genotype probs at masked loci for masked individuals
  //TODO: parallelize - requires reducing sums in GenotypeProbOutputter class
  if(options.OutputCGProbs()){
    std::string s = options.getResultsDir();
    s.append("/PPGenotypeProbs.txt");
    IC->OutputCGProbs(s.c_str());
  }
}
void HapMixModel::InitialiseTests(AdmixOptions& options, const InputData& data, LogWriter& Log){
  const bool isMaster = Comms::isMaster();
  //const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

    if( options.getScoreTestIndicator() && (isMaster || isWorker)){
  Scoretests.Initialise(&options, IC, &Loci, data.GetPopLabels(), Log);
  }
  if(isMaster)
    InitializeErgodicAvgFile(&options, Log, data.GetPopLabels(), data.getCovariateLabels());
  //}
}
//this function is here because three different objects have to write to avgstream
void HapMixModel::InitializeErgodicAvgFile(const AdmixOptions* const options, LogWriter &Log,  
				     const Vector_s& , const Vector_s& CovariateLabels){
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
	avgstream << "sumIntensities.mean\t";
      }
      
    }//end if hierarchical model

    //rate parameter of prior on frequency Dirichlet prior params
    if(!options->getFixedAlleleFreqs()){
      avgstream << "FreqPriorRate\t"; 
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
    if(options->getChibIndicator()){// chib calculation
      avgstream << "LogPrior\tLogPosterior\tLogPosteriorAdmixture\tLogPosteriorSumIntensities\t"
		 << "LogPosteriorAlleleFreqs\tLogMarginalLikelihood";
    }
    avgstream << "\n";
  } else {
    Log << "Not writing ergodic averages to file\n";
  }
}

double HapMixModel::getDevianceAtPosteriorMean(const AdmixOptions* const options, LogWriter& Log){
  return IC->getDevianceAtPosteriorMean(options, R, &Loci, Log, L->getSumLogRho(), Loci.GetNumberOfChromosomes(), &A);
}
