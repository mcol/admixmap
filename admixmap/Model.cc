#include "admixmap.h"

Model::Model(){

}

Model::~Model(){
  for(unsigned r = 0; r < R.size(); ++r)delete R[r];
  R.clear();
  delete IC;
  if(avgstream.is_open())avgstream.close();
}

void Model::Initialise(Genome& Loci, AdmixOptions& options, InputData& data,  LogWriter& Log){
  const bool isMaster = Comms::isMaster();
  //  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();
  
  A.Initialise(&options, &data, &Loci, Log); //checks allelefreq files, initialises allele freqs and finishes setting up Composite Loci
  
  IC = new IndividualCollection(&options, &data, &Loci);//NB call after A Initialise;//and before L and R Initialise
  if(isMaster || isWorker)IC->LoadData(&options, &data, (!options.getTestForAdmixtureAssociation() && options.getPopulations() > 1));    
  if(isWorker)IC->setGenotypeProbs(&Loci, &A); // sets unannealed probs
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

void Model::Iterate(const int & samples, const int & burnin, const double* Coolnesses, unsigned coolness,
		    AdmixOptions & options, InputData & data, const Genome & Loci, LogWriter& Log, 
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
    UpdateParameters(iteration, &options, &Loci, Log, data.GetPopLabels(), Coolnesses[coolness], AnnealedRun);
    SubIterate(iteration, burnin, options, data, Loci, Log, SumEnergy, SumEnergySq, 
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
  A.resetStepSizeApproximator(resetk);
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
