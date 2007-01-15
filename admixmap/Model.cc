#include "Model.h"
#include "Annealer.h"
#include "defs.h"

Model::Model(){

}

Model::~Model(){
  for(unsigned r = 0; r < R.size(); ++r)delete R[r];
  R.clear();
  delete IC;
  if(avgstream.is_open())avgstream.close();
}

void Model::InitialiseLoci(const Options& options, InputData& data, LogWriter& Log){
  Loci.Initialise(&data, options.getPopulations(), Log);//reads locusfile and creates CompositeLocus objects
  if(Comms::isFreqSampler()){
    //print table of loci for R script to read
    string locustable = options.getResultsDir();
    locustable.append("/LocusTable.txt");
    Loci.PrintLocusTable(locustable.c_str(), data.getLocusMatrix().getCol(1));
    locustable.clear();
  }
}

void Model::InitialiseRegressionObjects(AdmixOptions& options, InputData& data,  LogWriter& Log){
  if( Comms::isMaster() ){
    Regression::OpenOutputFile(options.getNumberOfOutcomes(), options.getRegressionOutputFilename(), Log);  
    Regression::OpenExpectedYFile(options.getEYFilename(), Log);
  }
  for(int r = 0; r < options.getNumberOfOutcomes(); ++r){
    //determine regression type and allocate regression objects
    if( data.getOutcomeType(r)== Binary ) R.push_back( new LogisticRegression() );
    else if( data.getOutcomeType(r)== Continuous ) R.push_back( new LinearRegression());
    else if( data.getOutcomeType(r)== CoxData ) R.push_back(new CoxRegression());
    
    if(Comms::isMaster()) {
      if(R[r]->getRegressionType()==Cox)
	R[r]->Initialise(r, options.getRegressionPriorPrecision(), IC->getCovariatesMatrix(),data.getCoxOutcomeVarMatrix(), Log);
      else
	R[r]->Initialise(r, options.getRegressionPriorPrecision(), IC->getCovariatesMatrix(), IC->getOutcomeMatrix(), Log);
    }
    else R[r]->Initialise(r, IC->GetNumCovariates());
    R[r]->InitializeOutputFile(data.getCovariateLabels(), options.getNumberOfOutcomes());
  }
}

void Model::Run(AdmixOptions& options, InputData& data, LogWriter& Log){
  const bool isMaster = Comms::isMaster();
#ifdef PARALLEL
  double t1 = 0, t2 = 0, t3 = 0, t0 = MPI::Wtime();
#endif
#ifdef __ADMIXMAP__
   //  ******** single individual, one population, fixed allele frequencies  ***************************
    if( getNumIndividuals() == 1 && options.getPopulations() == 1 && strlen(options.getAlleleFreqFilename()) )
      // nothing to do except calculate likelihood
      getOnePopOneIndLogLikelihood(Log, data.GetPopLabels());
    else {
#endif
          // ******************* Initialize test objects and ergodicaveragefile *******************************
      InitialiseTests(options, data, Log);

      //open file to output loglikelihood
      string s = options.getResultsDir()+"/loglikelihoodfile.txt";
      ofstream loglikelihoodfile(s.c_str());
    
      // ******************* Set annealing schedule ************************************************
      double SumEnergy = 0.0, SumEnergySq = 0.0;
      int samples = options.getTotalSamples();
      int burnin = options.getBurnIn();
      int NumAnnealedRuns = options.getNumAnnealedRuns();
      if( options.getTestOneIndivIndicator() )NumAnnealedRuns = 0;
      double coolness = 1.0; // default

      s = options.getResultsDir()+"/annealmon.txt";
      Annealer A(options.getThermoIndicator(), NumAnnealedRuns, samples, burnin, s.c_str());
      double AISsumlogz = 0.0; //for computing marginal likelihood by Annealed Importance Sampling
    
       // set number of samples : 1 for annealing runs, "samples" option otherwise. Overriden for final, unannealed run with "thermo" option

      //if( options.getTestOneIndivIndicator() )NumAnnealedRuns = 0;
      A.PrintRunLengths(Log, options.getTestOneIndivIndicator());
    
      A.SetAnnealingSchedule();
    
      //Write initial values
      //     if(options.getIndAdmixHierIndicator()  ){
      //       if(options.getDisplayLevel()>2)Log.setDisplayMode(On);
      //       else Log.setDisplayMode(Quiet);
      //       //Log << "InitialParameterValues:\n"
      //       //OutputParameters(-1, &L, &A, R, &options, Log);
      //       //Log << "\n";
      //     }
#ifdef PARALLEL
      t1 = MPI::Wtime()-t0;
#endif
#ifdef __ADMIXMAP__
      if(!options.getTestOneIndivIndicator()) { 
#endif
	for(int run=0; run <= NumAnnealedRuns; ++run) { //loop over coolnesses from 0 to 1
	  // should call a posterior mode-finding algorithm before last run at coolness of 1
	  //resets for start of each run
	  SumEnergy = 0.0;//cumulative sum of modified loglikelihood
	  SumEnergySq = 0.0;//cumulative sum of square of modified loglikelihood
	  bool AnnealedRun = A.SetRunLengths(run, &samples, &burnin, &coolness);


	  if(NumAnnealedRuns > 0) {
	    cout <<"\rSampling at coolness of " << coolness << "       "<< flush;
	    // reset approximation series in step size tuners
	    int resetk = NumAnnealedRuns; //   
	    if(samples < NumAnnealedRuns) {// samples=1 if annealing without thermo integration
	      resetk = 1 + run;
	    }
	    ResetStepSizeApproximators(resetk);
	    
	  }
	  // accumulate scalars SumEnergy and SumEnergySq at this coolness
	  // array Coolnesses is not used unless TestOneIndivIndicator is true
	  Iterate(samples, burnin, A.GetCoolnesses(), run, options, data, Log, SumEnergy, SumEnergySq, AISsumlogz,
			 AnnealedRun, loglikelihoodfile);

#ifdef PARALLEL
	  t2 = MPI::Wtime()-t1;
#endif
	  if(isMaster){	
	    //calculate mean and variance of energy at this coolness
	    A.CalculateLogEvidence(run, coolness, SumEnergy, SumEnergySq, samples - burnin);
	  } 
	} // end loop over coolnesses
#ifdef __ADMIXMAP__
      } 
      else { // evaluate energy for test individual only at all coolnesses simultaneously
	// call with argument AnnealedRun false - copies of test individual will be annealed anyway  
	Iterate(samples, burnin, A.GetCoolnesses(), NumAnnealedRuns, options, data, Log, SumEnergy, SumEnergySq, AISsumlogz, false,
		       loglikelihoodfile);
	// arrays of accumulated sums for energy and energy-squared have to be retrieved by function calls
	A.CalculateLogEvidence(getSumEnergy(), getSumEnergySq(), options.getNumAnnealedRuns());
	

      } // end evaluation of test individual
#endif
      if(isMaster)
	cout<< "\nIterations completed                       \n" << flush;
    
      // *************************** OUTPUT AT END ***********************************************************

      if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);	
      else Log.setDisplayMode(Quiet);

      Finalize(options, Log, data);
      
      A.PrintResults(Log, getDevianceAtPosteriorMean(&options, Log));
		
      if(options.getThermoIndicator()){
	Log << "\nAnnealed Importance Sampling estimates log marginal likelihood as " << AISsumlogz << "\n";
      }
		
      //Expected Outcome
      if(isMaster && options.getNumberOfOutcomes() > 0){
	Regression::FinishWritingEYAsRObject((options.getTotalSamples()-options.getBurnIn())/ options.getSampleEvery(), 
					     data.getOutcomeLabels());
      }
#ifdef __ADMIXMAP__
    }//end else
#endif
     if(isMaster)cout << "Output to files completed\n" << flush;
    
    // ******************* acceptance rates - output to screen and log ***************************
    if( options.getIndAdmixHierIndicator() ){
      PrintAcceptanceRates(options, Log);
    }
#ifdef PARALLEL
    t3 = MPI::Wtime()-t2;
#endif
}

void Model::Iterate(const int & samples, const int & burnin, const double* Coolnesses, unsigned coolness,
		    AdmixOptions & options, InputData & data, LogWriter& Log, 
		    double & SumEnergy, double & SumEnergySq, 
		    double& AISsumlogz, bool AnnealedRun, ofstream & loglikelihoodfile) {
  const bool isMaster = Comms::isMaster();
  //const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  //Accumulate Energy
  double Energy = 0.0;
  double AISz = 0.0;
  if(isMaster && !AnnealedRun) cout << endl;

  for( int iteration = 0; iteration <= samples; iteration++ ) {
    if( (isMaster || isWorker) && iteration > burnin) {
      //accumulate energy as minus loglikelihood, calculated using unnanealed genotype probs
      if( !options.getTestOneIndivIndicator() ) {
	Energy = IC->getEnergy(&options, R, AnnealedRun); // should store loglikelihood if not AnnealedRun
	if(isMaster){
	  SumEnergy += Energy;
	  SumEnergySq += Energy*Energy;
	  if(coolness>0)AISz += exp((Coolnesses[coolness]-Coolnesses[coolness-1])*(-Energy));
	  // write to file if not AnnealedRun
	  if(!AnnealedRun){
	      loglikelihoodfile << iteration<< "\t" << Energy <<endl;
	      if(options.getDisplayLevel()>2 && !(iteration%options.getSampleEvery()))cout << Energy << "\t";
	  }
	}
      } else {  
	IC->accumulateEnergyArrays(&options);
      }
    }
    //Write Iteration Number to screen
    if( isMaster && !AnnealedRun &&  !(iteration % options.getSampleEvery()) ) {
      WriteIterationNumber(iteration, (int)log10((double) samples+1 ), options.getDisplayLevel());
    }
    
    // if annealed run, anneal genotype probs - for testindiv only if testsingleindiv indicator set in IC
    if((isMaster || isWorker) && (AnnealedRun || options.getTestOneIndivIndicator())) 
      IC->annealGenotypeProbs(Loci.GetNumberOfChromosomes(), Coolnesses[coolness], Coolnesses); 

    //Sample Parameters    
    UpdateParameters(iteration, &options, Log, data.GetPopLabels(), Coolnesses[coolness], AnnealedRun);
    SubIterate(iteration, burnin, options, data, Log, SumEnergy, SumEnergySq, 
	       AnnealedRun);
	
  }// end loop over iterations
  //use Annealed Importance Sampling to calculate marginal likelihood
  if(coolness>0) AISsumlogz += log(AISz /= (double)(samples-burnin));

}

// void Model::UpdateParameters(int iteration, const AdmixOptions *options, const Genome *Loci, LogWriter& Log, 
// 		      const Vector_s& PopulationLabels, double coolness, bool anneal){
//   const bool isMaster = Comms::isMaster();
//   const bool isFreqSampler = Comms::isFreqSampler();
//   const bool isWorker = Comms::isWorker();

//   A.ResetAlleleCounts(options->getPopulations());

//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   if(isMaster || isWorker){
//     // ** update global sumintensities conditional on genotype probs and individual admixture proportions
//      if((options->getPopulations() > 1) && options->getIndAdmixHierIndicator() && !options->getHapMixModelIndicator() && 
//         (Loci->GetLengthOfGenome() + Loci->GetLengthOfXchrm() > 0.0))
//        L->UpdateGlobalSumIntensities(IC, (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1)); 
//     // leaves individuals with HMM probs bad, stored likelihood ok
//     // this function also sets locus correlations in Chromosomes
//   }
// //   if(options->getHapMixModelIndicator() && (iteration%2)){
// //     if(isMaster || isWorker){
// //       L->UpdateSumIntensitiesByRandomWalk(IC,  
// // 			      (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1) );
// //     }
// //   }
//   //find posterior modes of individual admixture at end of burn-in
//   //set Chib numerator
//   if(!anneal && iteration == options->getBurnIn() && (options->getChibIndicator() || strlen(options->getIndAdmixModeFilename()))) {
//     IC->FindPosteriorModes(options, R, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), &A, PopulationLabels);
//     if( options->getChibIndicator() ) {
//       IC->setChibNumerator(options, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), &A);
//     }
//   }

//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//   // Update individual-level parameters, sampling locus ancestry states
//   // then update jump indicators (+/- num arrivals if required for conjugate update of admixture or rho
//   if(isMaster || isWorker){
//   IC->SampleLocusAncestry(iteration, options, R, L->getpoptheta(), L->getalpha(), Scoretests.getAffectedsOnlyTest(), anneal);
//    }

//   if(isWorker || isFreqSampler) {
// #ifdef PARALLEL
//     MPE_Log_event(13, iteration, "sampleHapPairs");
// #endif
//     IC->SampleHapPairs(options, &A, Loci, anneal); // loops over individuals to sample hap pairs then increment allele counts
// #ifdef PARALLEL
//     MPE_Log_event(14, iteration, "sampledHapPairs");
// #endif
//   }
  
// #ifdef PARALLEL
//   if(isWorker || isFreqSampler){
//     A.SumAlleleCountsOverProcesses(options->getPopulations());
//   }
// #endif
  
//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   if( (isMaster || isWorker) && !anneal && iteration > options->getBurnIn() ){
//     //score tests
//     if( options->getScoreTestIndicator() ){
// #ifdef PARALLEL
//       MPE_Log_event(21, iteration, "ScoreTestUpdatestart");
// #endif
//       //TODO: broadcast dispersion in linear regression
//       Scoretests.Update(R);//score tests evaluated for first outcome var only
// #ifdef PARALLEL
//       MPE_Log_event(22, iteration, "ScoreTestUpdateEnd");
// #endif
//     }
//     if(options->getTestForResidualAllelicAssoc()){
// #ifdef PARALLEL
//       MPE_Log_event(21, iteration, "ResLDTestStart");
// #endif
//       Scoretests.UpdateScoresForResidualAllelicAssociation(A.GetAlleleFreqs());
// #ifdef PARALLEL
//       MPE_Log_event(22, iteration, "ResLDTestEnd");
// #endif
//     }
//   }
//       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
//   // sample individual admixture and sum-intensities 
//   IC->SampleParameters(iteration, options, R, L->getpoptheta(), L->getalpha(),  
// 		       L->getrhoalpha(), L->getrhobeta(), anneal);
//   // stored HMM likelihoods will now be bad if the sum-intensities are set at individual level
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
//   if( options->getChibIndicator() && !anneal && iteration >= options->getBurnIn() )
//     IC->updateChib(options, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), &A);      
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//   // update allele frequencies conditional on locus ancestry states
//   // TODO: this requires fixing to anneal allele freqs for historicallelefreq model
//   if( !options->getFixedAlleleFreqs()){
//     if(isFreqSampler){
// #ifdef PARALLEL
//       MPE_Log_event(7, iteration, "SampleFreqs");
// #endif
//       A.Update(IC, (iteration > options->getBurnIn() && !anneal), coolness);
// #ifdef PARALLEL
//     MPE_Log_event(8, iteration, "SampledFreqs");
// #endif
//     }
// #ifdef PARALLEL
//     if(isWorker || isFreqSampler ) { 
//       A.BroadcastAlleleFreqs();
//     }
// #endif
//   }//end if random allele freqs
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   // even for fixed allele freqs, must reset annealed genotype probs as unnannealed

//   if(isWorker && (!options->getFixedAlleleFreqs() || anneal)) {   
// #ifdef PARALLEL
//     MPE_Log_event(11, iteration, "setGenotypeProbs"); 
// #endif
//     IC->setGenotypeProbs(Loci, &A); // sets unannealed probs ready for getEnergy
//     IC->HMMIsBad(true); // update of allele freqs sets HMM probs and stored loglikelihoods as bad
// #ifdef PARALLEL
//     MPE_Log_event(12, iteration, "GenotypeProbsSet"); 
// #endif
//   } // update of allele freqs sets HMM probs and stored loglikelihoods as bad
  
//   // next update of stored loglikelihoods will be from getEnergy if not annealing run, from updateRhowithRW if globalrho, 
//   // or from update of individual-level parameters otherwise
  
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   if(isMaster || isWorker){
//     if(options->getHapMixModelIndicator()){
//       //L->UpdateGlobalTheta(iteration, IC);
//       //Hamiltonian Sampler, using sampled ancestry states. NB: requires accumulating of SumAncestry in IC
//       L->SampleHapMixLambda(IC->getSumAncestry(),  
// 			    (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1) );
//     }
    
//     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     else
//     if(!options->getHapMixModelIndicator())
//       //update population admixture Dirichlet parameters conditional on individual admixture
//       L->UpdatePopAdmixParams(iteration, IC, Log);
//   }
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
//   // ** update regression parameters (if regression model) conditional on individual admixture
//   if(isMaster){
//     bool condition = (!anneal && iteration > options->getBurnIn() && isMaster);
//     for(unsigned r = 0; r < R.size(); ++r){
//       R[r]->Update(condition, IC->getOutcome(r), coolness );
//       //output expected values of outcome variables to file every 'every' iterations after burnin
//       if(condition && !(iteration % options->getSampleEvery()) ) {
// 	R[r]->OutputExpectedY();
//       }
//     }
//   }
// #ifdef PARALLEL
//   if(isMaster || isWorker){    //broadcast regression parameters to workers
//     for(unsigned r = 0; r < R.size(); ++r){
//       Comms::BroadcastRegressionParameters(R[r]->getbeta(), R[r]->getNumCovariates());
//       //if(R[r]->getRegressionType()==Linear) Comms::BroadcastRegressionPrecision(R[r]->getlambda());
//     }
//   }
// #endif
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
// }

void Model::ResetStepSizeApproximators(int resetk){
  IC->resetStepSizeApproximators(resetk); 
  pA->resetStepSizeApproximator(resetk);
}

// double Model::getDevianceAtPosteriorMean(const AdmixOptions* const options, Genome* Loci, LogWriter& Log){
//   return IC->getDevianceAtPosteriorMean(options, R, Loci, Log, L->getSumLogRho(), Loci->GetNumberOfChromosomes(), &A);
// }


// void Model::InitialiseTests(AdmixOptions& options, const InputData& data, const Genome& Loci, LogWriter& Log){
//   const bool isMaster = Comms::isMaster();
//   //const bool isFreqSampler = Comms::isFreqSampler();
//   const bool isWorker = Comms::isWorker();

//   if( options.getScoreTestIndicator() && (isMaster || isWorker)){
//     Scoretests.Initialise(&options, IC, &Loci, data.GetPopLabels(), Log);
//   }
//   //if(isMaster){
//   if( options.getTestForDispersion() ){
//     DispTest.Initialise(&options, Log, Loci.GetNumberOfCompositeLoci());    
//   }
//   if(options.getStratificationTest())
//       StratTest.OpenOutputFile(options.getStratTestFilename(), Log);
//   if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
//     AlleleFreqTest.Initialise(&options, &Loci, Log );  
//   if( options.getHWTestIndicator() )
//     HWtest.Initialise(&options, Loci.GetTotalNumberOfLoci(), Log);
//   if(isMaster)
//     InitializeErgodicAvgFile(&options, Log, data.GetPopLabels(), data.getCovariateLabels());
//   //}
// }
//this function is here because three different objects have to write to avgstream
// void Model::InitializeErgodicAvgFile(const AdmixOptions* const options, LogWriter &Log,  
// 				     const Vector_s& PopulationLabels, const Vector_s& CovariateLabels){
//   Log.setDisplayMode(Quiet);
//   //Open ErgodicAverageFile  
//   if( strlen( options->getErgodicAverageFilename() ) ) {
//     avgstream.open( options->getErgodicAverageFilename(), ios::out );
//     if( ! avgstream ){
//       Log.setDisplayMode(On);
//       Log << "ERROR: Couldn't open Ergodic Average file\n";
//       //exit( 1 );
//     } else {
//       Log << "Writing ergodic averages of parameters to "
// 	  << options->getErgodicAverageFilename() << "\n\n";
//     }
    
//     // Header line of ergodicaveragefile
//     if( options->getIndAdmixHierIndicator() ){
//       if(options->getPopulations()>1){
// 	if(!options->getHapMixModelIndicator())
// 	  for( int i = 0; i < options->getPopulations(); i++ ){
// 	    avgstream << PopulationLabels[i] << "\t";
// 	  }
// 	if( options->isGlobalRho() ) avgstream << "sumIntensities\t";
// 	else avgstream << "sumIntensities.mean\t";
//       }
      
//       // dispersion parameters
//       if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
// 	for( int k = 0; k < options->getPopulations(); k++ ){
// 	  avgstream << "eta" << k << "\t";
// 	}
//       }
//     }//end if hierarchical model

//     //rate parameter of prior on frequency Dirichlet prior params
//     if(options->getHapMixModelIndicator() && !options->getFixedAlleleFreqs()){
//       avgstream << "FreqPriorRate\t"; 
//     }

//     // Regression parameters
//     if( options->getNumberOfOutcomes() > 0 ){
//       for(int r = 0; r < options->getNumberOfOutcomes(); ++r){
// 	if(IC->getOutcomeType(r)!=CoxData)
// 	avgstream << "intercept\t";

// 	//write covariate labels to header
// 	  copy(CovariateLabels.begin(), CovariateLabels.end(), ostream_iterator<string>(avgstream, "\t")); 

// 	if( IC->getOutcomeType(r)==Continuous )//linear regression
// 	  avgstream << "precision\t";
//       }
//     }

//     avgstream << "MeanDeviance\tVarDeviance\t";
//     if(options->getChibIndicator()){// chib calculation
//       avgstream << "LogPrior\tLogPosterior\tLogPosteriorAdmixture\tLogPosteriorSumIntensities\t"
// 		 << "LogPosteriorAlleleFreqs\tLogMarginalLikelihood";
//     }
//     avgstream << "\n";
//   } else {
//     Log << "No ergodicaveragefile given\n";
//   }
// }

void Model::OutputErgodicAvgDeviance(int samples, double & SumEnergy, double & SumEnergySq) {
  double EAvDeviance, EVarDeviance;
  EAvDeviance = 2.0*SumEnergy / (double) samples;//ergodic average of deviance
  EVarDeviance = 4.0 * SumEnergySq / (double)samples - EAvDeviance*EAvDeviance;//ergodic variance of deviance 
  avgstream << EAvDeviance << " "<< EVarDeviance <<" ";
}
