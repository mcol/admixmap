/* 
 *   ADMIXMAP
 *   AdmixMapModel.cc
 *   
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "AdmixMapModel.h"
#include "AdmixFilenames.h"
#include "bclib/LogWriter.h"

using namespace::bclib;

// AdmixMapModel::AdmixMapModel(){

// }

AdmixMapModel::~AdmixMapModel(){
    delete L;
}

void AdmixMapModel::Initialise(AdmixOptions& options, InputAdmixData& data,  LogWriter& Log){

  InitialiseGenome(Loci, options, data, Log);

  A.Initialise(&options, &data, &Loci, Log, options.getChibIndicator()); //checks allelefreq files, initialises allele freqs and finishes setting up Composite Loci
  pA = &A;//set pointer to AlleleFreqs object
  
  AdmixedIndividuals = new AdmixIndividualCollection(&options, &data, &Loci);//NB call after A Initialise;//and before L and R Initialise
  IC = (IndividualCollection*) AdmixedIndividuals;
  AdmixedIndividuals->LoadData(&options, &data);    
  AdmixedIndividuals->setGenotypeProbs(&Loci); // sets unannealed probs

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

  
  L = new PopAdmix(options, Loci);    
  L->Initialise(IC->getSize(), data.GetPopLabels(), Log);
  A.PrintPrior(data.GetPopLabels(), Log);
  
  if( options.getPopulations()>1 && (options.isGlobalRho())) {
    Loci.SetLocusCorrelation(L->getrho());
  }
  
  AdmixedIndividuals->Initialise(options, Loci, data.GetPopLabels(), Log);
  AdmixedIndividuals->DrawInitialAdmixture(L->getalpha());
  
  //initialise regression objects
  if (options.getNumberOfOutcomes()>0 ){
    InitialiseRegressionObjects(options, data, Log);
  }
}

void AdmixMapModel::TestIndivRun(Options& options, InputData& data, LogWriter& Log){
  double SumEnergy = 0.0;//cumulative sum of modified loglikelihood
  double SumEnergySq = 0.0;//cumulative sum of square of modified loglikelihood
  int samples = options.getTotalSamples();
  int burnin = options.getBurnIn();

  Start(options, data, Log);

  // call with argument AnnealedRun false - copies of test individual will be annealed anyway  
  Iterate(samples, burnin, _Annealer.GetCoolnesses(), 0, options, data, Log, 
	  SumEnergy, SumEnergySq, false);
  
  // arrays of accumulated sums for energy and energy-squared have to be retrieved by function calls
  _Annealer.CalculateLogEvidence(getSumEnergy(), getSumEnergySq(), options.getNumAnnealedRuns());

  Finish(options, data, Log);    
}

void AdmixMapModel::Iterate(const int & samples, const int & burnin, const double* Coolnesses, unsigned coolness,
			    Options & options, InputData & data, LogWriter& Log, 
			    double & SumEnergy, double & SumEnergySq,  bool AnnealedRun) {

  //Accumulate Energy
  double AISz = 0.0;
  if(!AnnealedRun) cout << endl;

  for( int iteration = 0; iteration <= samples; iteration++ ) {
    if( iteration > burnin) {
      if( options.getTestOneIndivIndicator() ) {
	AdmixedIndividuals->accumulateEnergyArrays(options);
      }
      AccumulateEnergy(Coolnesses, coolness, options, SumEnergy, SumEnergySq, AISz, AnnealedRun, iteration );
    }
  

    //Write Iteration Number to screen
    if( !AnnealedRun &&  !(iteration % options.getSampleEvery()) ) {
      WriteIterationNumber(iteration, (int)log10((double) samples+1 ), options.getDisplayLevel());
    }
    
    //Sample Parameters
    UpdateParameters(iteration, options, Log, data.GetHiddenStateLabels(), Coolnesses, Coolnesses[coolness], AnnealedRun);
    SubIterate(iteration, burnin, options, data, Log, SumEnergy, SumEnergySq,
	       AnnealedRun);
    
  }// end loop over iterations
  //use Annealed Importance Sampling to calculate marginal likelihood
  if(coolness>0) AISsumlogz += log(AISz /= (double)(samples-burnin));

}

void AdmixMapModel::UpdateParameters(int iteration, const Options& _options, LogWriter& Log, 
		      const Vector_s& PopulationLabels, const double* Coolnesses, double coolness, bool anneal){
    //cast Options pointer to AdmixOptions pointer for access to ADMIXMAP options
  //TODO: change pointer to reference
  const AdmixOptions& options = (const AdmixOptions&) _options;

  // if annealed run, anneal genotype probs - for testindiv only if testsingleindiv indicator set in IC
  if(anneal || options.getTestOneIndivIndicator())
    AdmixedIndividuals->annealGenotypeProbs(Loci.GetNumberOfChromosomes(), coolness, Coolnesses);


  A.ResetAlleleCounts(options.getPopulations());

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ** update global sumintensities conditional on genotype probs and individual admixture proportions
  if((options.getPopulations() > 1) && options.getIndAdmixHierIndicator() && 
     (Loci.GetLengthOfGenome() + Loci.GetLengthOfXchrm() > 0.0))
    L->UpdateGlobalSumIntensities(AdmixedIndividuals, (!anneal && iteration > options.getBurnIn() && options.getPopulations() > 1));
  // leaves individuals with HMM probs bad, stored likelihood ok
  // this function also sets locus correlations in Chromosomes
  

  //find posterior modes of individual admixture at end of burn-in
  //set Chib numerator
  if(!anneal && iteration == options.getBurnIn() && (options.getChibIndicator() || strlen(options.getIndAdmixModeFilename()))) {
    AdmixedIndividuals->FindPosteriorModes(options, R, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), &A, PopulationLabels);
    if( options.getChibIndicator() ) {
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
  AdmixedIndividuals->HMMUpdates(iteration, options, R, L->getpoptheta(), L->getalpha(), 
				   Scoretests.getAffectedsOnlyTest(), Scoretests.getAncestryAssocTest(), anneal);
  
  // loops over individuals to sample hap pairs then increment allele counts, skipping missing genotypes
  AdmixedIndividuals->SampleHapPairs(options, &A, &Loci, true, anneal, true);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if( !anneal && iteration > options.getBurnIn() ){
    //score tests
    if( options.getScoreTestIndicator() ){
      //TODO: broadcast dispersion in linear regression
      Scoretests.Update(R);//score tests evaluated for first outcome var only
    }
    if(options.getTestForResidualAllelicAssoc()){
      ResidualAllelicAssocScoreTest.Reset();
      ResidualAllelicAssocScoreTest.Update(A.GetAlleleFreqs(), options.getHapMixModelIndicator());
    }
  }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // sample individual admixture and sum-intensities 
  AdmixedIndividuals->SampleParameters(iteration, options, R, L->getpoptheta(), L->getalpha(),
		       L->getrhoalpha(), L->getrhobeta(), Scoretests.getAncestryAssocTest(),anneal);
  // stored HMM likelihoods will now be bad if the sum-intensities are set at individual level
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if( options.getChibIndicator() && !anneal && iteration >= options.getBurnIn() )
    AdmixedIndividuals->updateChib(options, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), &A);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // update allele frequencies conditional on locus ancestry states
  // TODO: this requires fixing to anneal allele freqs for historicallelefreq model
  if( !options.getFixedAlleleFreqs()){
    A.Update(IC, (iteration > options.getBurnIn() && !anneal), coolness);
  }//end if random allele freqs
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // even for fixed allele freqs, must reset annealed genotype probs as unnannealed

  if(!options.getFixedAlleleFreqs() || anneal) {   
    AdmixedIndividuals->setGenotypeProbs(&Loci); // sets unannealed probs ready for getEnergy
    AdmixedIndividuals->HMMIsBad(true); // update of allele freqs sets HMM probs and stored loglikelihoods as bad
  } // update of allele freqs sets HMM probs and stored loglikelihoods as bad
  
  // next update of stored loglikelihoods will be from getEnergy if not annealing run, from updateRhowithRW if globalrho, 
  // or from update of individual-level parameters otherwise
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //update population admixture Dirichlet parameters conditional on individual admixture
  L->UpdatePopAdmixParams(iteration, AdmixedIndividuals, Log);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // ** update regression parameters (if regression model) conditional on individual admixture
  bool condition = (!anneal && iteration > options.getBurnIn());
  for(unsigned r = 0; r < R.size(); ++r){
    R[r]->Update(condition, IC->getOutcome(r), coolness );
    //output expected values of outcome variables to file every 'every' iterations after burnin
    if(condition && !(iteration % options.getSampleEvery()) ) {
      R[r]->OutputExpectedY();
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
void AdmixMapModel::ResetStepSizeApproximators(int resetk){
  Model::ResetStepSizeApproximators(resetk);
 L->resetStepSizeApproximator(resetk);
}

void AdmixMapModel::SubIterate(int iteration, const int & burnin, Options& _options, InputData & data, 
			       LogWriter& Log, double& SumEnergy, double& SumEnergySq, 
			       bool AnnealedRun){
  //cast Options object to AdmixOptions for access to ADMIXMAP options
  AdmixOptions& options = (AdmixOptions&) _options;

    //Output parameters to file and to screen
    Log.setDisplayMode(Quiet);
     if(!AnnealedRun){    
      // output every 'getSampleEvery()' iterations
      if(!(iteration % options.getSampleEvery()))
	OutputParameters(iteration, &options, Log);
      
      // ** set merged haplotypes for allelic association score test 
      if( iteration == options.getBurnIn() ){

	if(options.getTestForHaplotypeAssociation())
	  Scoretests.MergeRareHaplotypes(L->getalpha0());

	if( options.getStratificationTest() )
	  StratTest.Initialize( &options, Loci, IC, Log);
      }
	    
      //Updates and Output after BurnIn     
      if( !AnnealedRun && iteration > burnin){
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
	  if( options.getScoreTestIndicator() )  
	    Scoretests.Output(data.GetHiddenStateLabels());
	  if(options.getTestForResidualAllelicAssoc())
	    ResidualAllelicAssocScoreTest.Output(data.getLocusLabels());
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
   
   if(options.getPopulations() > 1 ){
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
   A.OutputAlleleFreqSamplerAcceptanceRates(options.getResultsDir());
 }
}

void AdmixMapModel::Finalize(const Options& _options, LogWriter& Log, const InputData& data){
  const AdmixOptions& options = (const AdmixOptions&)_options;
  if(AdmixedIndividuals->getSize()==1){
    if( options.getChibIndicator()) {
      //MLEs of admixture & sumintensities used in Chib algorithm to estimate marginal likelihood
      AdmixedIndividuals->OutputChibResults(Log);
    }
  }
  else{
    AdmixedIndividuals->WritePosteriorMeans(options, data.GetHiddenStateLabels());
  }
  //FST
  if( strlen( options.getHistoricalAlleleFreqFilename()) ){
    A.OutputFST();
  }
  A.CloseOutputFile((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery(), 
		    data.GetHiddenStateLabels());

  if(options.getOutputAlleleFreq()){
//     ofstream klinfofile((options.getResultsDir() + "/KLinfo.txt").c_str());
//     A.WriteKLInfo((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery(), klinfofile);
//     klinfofile.close();

    A.WriteLocusInfo((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery(), options.getResultsDir(), data.GetHiddenStateLabels());
  }

  //stratification test
  if( options.getStratificationTest() ) StratTest.Output(Log);
  //dispersion test
  if( options.getTestForDispersion() )  DispTest.Output(options.getTotalSamples() - options.getBurnIn(), 
							Loci, data.GetHiddenStateLabels());
  //tests for mis-specified allele frequencies
  if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
    AlleleFreqTest.Output(options.getResultsDir(), &Loci, data.GetHiddenStateLabels(), Log); 
  //test for H-W eq
  if( options.getHWTestIndicator() )
    HWtest.Output( (options.getResultsDir() + "/" + HARDY_WEINBERG_TEST).c_str(), data.getLocusLabels(), Log); 

  if( options.getScoreTestIndicator() ) {
    //finish writing score test output as R objects
    //Scoretests.ROutput();
    //write final tables
    Scoretests.WriteFinalTables(data.GetHiddenStateLabels(), Log);
  }
  if(options.getTestForResidualAllelicAssoc())
    ResidualAllelicAssocScoreTest.WriteFinalTable((options.getResultsDir() + "/" + RESIDUAL_LD_TEST_FINAL).c_str(), data.getLocusLabels(), Log);

  //output to likelihood ratio file
  if(options.getTestForAffectedsOnly())
    Scoretests.OutputLikelihoodRatios(options.getResultsDir(), data.GetHiddenStateLabels());

  //print deviance at posterior mean, DIC
  _Annealer.PrintResults(Log, getDevianceAtPosteriorMean(options, Log));	

  //print results of annealed importance sampling
  if(options.getThermoIndicator()){
    Log << "\nAnnealed Importance Sampling estimates log marginal likelihood as " << AISsumlogz << "\n";
  }

}
void AdmixMapModel::InitialiseTests(Options& _options, const InputData& data, LogWriter& Log){

  //cast Options object to AdmixOptions for access to ADMIXMAP options
  AdmixOptions& options = (AdmixOptions&) _options;
  if( options.getScoreTestIndicator() ){
    Scoretests.Initialise(&options, IC, &Loci, data.GetHiddenStateLabels(), Log);
  }
  if(options.getTestForResidualAllelicAssoc())
    ResidualAllelicAssocScoreTest.Initialise((options.getResultsDir() + "/" + RESIDUAL_LD_TEST_PVALUES).c_str(), IC, &Loci);

  if( options.getTestForDispersion() ){
    DispTest.Initialise(options.getResultsDir(), Log, Loci.GetNumberOfCompositeLoci(), options.getPopulations());    
  }
  if(options.getStratificationTest())
      StratTest.OpenOutputFile(options.getResultsDir(), Log);
  if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
    AlleleFreqTest.Initialise(&options, &Loci, Log );  
  if( options.getHWTestIndicator() )
    HWtest.Initialise(Loci.GetTotalNumberOfLoci());
  InitializeErgodicAvgFile(&options, Log, data.GetHiddenStateLabels(), data.getCovariateLabels());


}
//this function is here because three different objects have to write to avgstream
void AdmixMapModel::InitializeErgodicAvgFile(const AdmixOptions* const options, LogWriter &Log,  
				     const Vector_s& PopulationLabels, const Vector_s& CovariateLabels){
  Log.setDisplayMode(bclib::Quiet);
  //Open ErgodicAverageFile  
  if( strlen( options->getErgodicAverageFilename() ) ) {
    avgstream.open( options->getErgodicAverageFilename(), ios::out );
    if( ! avgstream ){
      Log.setDisplayMode(bclib::On);
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


double AdmixMapModel::getDevianceAtPosteriorMean(const Options& options, LogWriter& Log){
  return AdmixedIndividuals->getDevianceAtPosteriorMean(options, R, &Loci, Log, L->getSumLogRho(), Loci.GetNumberOfChromosomes(), &A);
}
double* AdmixMapModel::getSumEnergy()const{
    return AdmixedIndividuals->getSumEnergy();
}
double* AdmixMapModel::getSumEnergySq()const{
    return AdmixedIndividuals->getSumEnergySq();
}

