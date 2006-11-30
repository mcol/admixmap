#include "admixmap.h"

// HapMixModel::HapMixModel(){

// }

HapMixModel::~HapMixModel(){
    delete L;
}

void HapMixModel::Initialise(Genome& Loci, AdmixOptions& options, InputData& data,  LogWriter& Log){
  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  Model::Initialise(Loci, options, data, Log);

  L = new PopHapMix(&options, &Loci);
  if(isMaster || isWorker)L->Initialise(IC->getSize(), data.GetPopLabels(), Log);
  if(isFreqSampler)A.PrintPrior(data.GetPopLabels(), Log);
  
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
  
  if(isMaster || isWorker)
    IC->Initialise(&options, &Loci, data.GetPopLabels(), Log);
  
  //initialise regression objects
  if (options.getNumberOfOutcomes()>0 && (isMaster || isWorker)){
    InitialiseRegressionObjects(options, data, Log);
  }
}

void HapMixModel::UpdateParameters(int iteration, const AdmixOptions *options, const Genome *Loci, LogWriter&, 
		      const Vector_s&, double coolness, bool anneal){
  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  A.ResetAlleleCounts(options->getPopulations());

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//   if(options->getHapMixModelIndicator() && (iteration%2)){
//     if(isMaster || isWorker){
//       L->UpdateSumIntensitiesByRandomWalk(IC,  
// 			      (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1) );
//     }
//   }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Update individual-level parameters, sampling locus ancestry states
  // then update jump indicators (+/- num arrivals if required for conjugate update of admixture or rho
  if(isMaster || isWorker){
  IC->SampleLocusAncestry(iteration, options, R, Scoretests.getAffectedsOnlyTest(), anneal);
   }

  if(isWorker || isFreqSampler) {
#ifdef PARALLEL
    MPE_Log_event(13, iteration, "sampleHapPairs");
#endif
    IC->SampleHapPairs(options, &A, Loci, anneal); // loops over individuals to sample hap pairs then increment allele counts
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
      Scoretests.UpdateScoresForResidualAllelicAssociation(A.GetAlleleFreqs(), true);
#ifdef PARALLEL
      MPE_Log_event(22, iteration, "ResLDTestEnd");
#endif
    }
  }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // sample individual admixture and sum-intensities 
  //  IC->SampleParameters(iteration, options, R, L->getpoptheta(), L->getalpha(),  
  //	       L->getrhoalpha(), L->getrhobeta(), anneal);
  // stored HMM likelihoods will now be bad if the sum-intensities are set at individual level
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // update allele frequencies conditional on locus ancestry states
  // TODO: this requires fixing to anneal allele freqs for historicallelefreq model
  if( !options->getFixedAlleleFreqs()){
    if(isFreqSampler){
#ifdef PARALLEL
      MPE_Log_event(7, iteration, "SampleFreqs");
#endif
      A.Update(IC, (iteration > options->getBurnIn() && !anneal), coolness, true);
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
    IC->setGenotypeProbs(Loci, &A); // sets unannealed probs ready for getEnergy
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
			     const Genome & , LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
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
	    if(options.getChibIndicator()) IC->OutputErgodicChib(&avgstream, options.getFixedAlleleFreqs());
	    avgstream << endl;
	  }
	  //Score Test output
	  if( options.getScoreTestIndicator() )  Scoretests.Output(iteration-burnin, data.GetPopLabels(), data.getLocusLabels(), false);
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
    // output individual and locus parameters every 'getSampleEvery()' iterations after burnin
    if( strlen( options->getIndAdmixtureFilename() ) ) IC->OutputIndAdmixture();
      A.OutputPriorParams();

    if(options->getOutputAlleleFreq() && !options->getHapMixModelIndicator()){
	A.OutputAlleleFreqs();
    }
  }
  // cout << endl;
}

void HapMixModel::PrintAcceptanceRates(const AdmixOptions& options, const Genome& , LogWriter& Log){
  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
  else Log.setDisplayMode(On);

  if(options.getPopulations() > 1 && Comms::isMaster()){
    L->printAcceptanceRates(Log);
  }
  A.OutputAlleleFreqSamplerAcceptanceRates((options.getResultsDir() + "/AlleleFreqSamplerAcceptanceRates.txt").c_str());
  if(!options.getFixedAlleleFreqs() && Comms::isFreqSampler()){
    Log << "Average expected Acceptance rate in allele frequency prior parameter sampler:\n" << A.getHapMixPriorSamplerAcceptanceRate()
	<< "\nwith average final step size of " << A.getHapMixPriorSamplerStepSize() << "\n";
  }
}

void HapMixModel::Finalize(const AdmixOptions& options, LogWriter& , const InputData& data, const Genome& ){
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
//output posterior means of lambda (expected number of arrivals)
  std::string s = options.getResultsDir();
  s.append("/lambdaPosteriorMeans.txt");
  L->OutputLambdaPosteriorMeans(s.c_str(), options.getTotalSamples()-options.getBurnIn());
//output final values of lambda
  const char* ss = options.getHapMixLambdaOutputFilename();
  if(strlen(ss))
      L->OutputLambda(ss);

//output final values of allelefreqs
  ss = options.getAlleleFreqOutputFilename();
  if(strlen(ss))
      A.OutputAlleleFreqs(ss);

}
void HapMixModel::InitialiseTests(AdmixOptions& options, const InputData& data, const Genome& Loci, LogWriter& Log){
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
    Log << "No ergodicaveragefile given\n";
  }
}

double HapMixModel::getDevianceAtPosteriorMean(const AdmixOptions* const options, Genome* Loci, LogWriter& Log){
  return IC->getDevianceAtPosteriorMean(options, R, Loci, Log, L->getSumLogRho(), Loci->GetNumberOfChromosomes(), &A);
}
