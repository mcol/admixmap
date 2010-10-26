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

//=============================================================================
/// \file AdmixMapModel.cc
/// Implementation of the AdmixMapModel class.
//=============================================================================

#include "AdmixMapModel.h"
#include "AdmixFilenames.h"
#include "bclib/LogWriter.h"
#include "DispersionFreqs.h"
#include "CorrelatedFreqs.h"

#define DEBUG_ITER_TIMES    0
#if DEBUG_ITER_TIMES
    #include "CodeTimer.h"
#endif

using namespace::bclib;

AdmixMapModel::AdmixMapModel(){
  L = 0;
  A = 0;
}

AdmixMapModel::~AdmixMapModel(){
  delete L;
  delete A;
}

void AdmixMapModel::Initialise(AdmixOptions& options, InputAdmixData& data,  LogWriter& Log){

  InitialiseGenome(Loci, options, data, Log);

  //choose allele frequency model:
  //correlated frequencies
  if(options.getCorrelatedAlleleFreqs())
    A = new CorrelatedFreqs;
  //dispersion model ('historic allele freqs')
  else if(strlen(options.getHistoricalAlleleFreqFilename()))
    A = new DispersionFreqs;
  //normal random- or fixed-frequency model
  else
    A = new AdmixFreqs;

  A->Initialise(&options, &data, &Loci, Log, options.getChibIndicator()); //checks allelefreq files, initialises allele freqs and finishes setting up Composite Loci
  pA = A;//set pointer to AlleleFreqs object

  AdmixedIndividuals = new AdmixIndividualCollection( options, data, Loci ); //NB call after A Initialise;//and before L and R Initialise
  IC = AdmixedIndividuals;
  AdmixedIndividuals->LoadData( options, data );
  AdmixedIndividuals->setGenotypeProbs(&Loci); // sets unannealed probs

  if ( ! (data.isPedFile() || options.getUsePedForInd()) ) {
      const int numdiploid = AdmixedIndividuals->getNumDiploidIndividuals();
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

  int NumObservations = 0;
  for ( int i =  AdmixedIndividuals->getSize(); i-- != 0 ; ) {
    PedBase & el = AdmixedIndividuals->getElement(i);
    NumObservations += el.getNumObs();
  }

  L = new PopAdmix(options, Loci);
  L->Initialise(NumObservations, data.GetPopLabels(), Log);
  A->PrintPrior(data.GetPopLabels(), Log);

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

  Start(options, data, Log);

  // call with argument AnnealedRun false - copies of test individual will be annealed anyway
  Iterate(_Annealer.GetCoolnesses(), 0, options, data, Log,
	  SumEnergy, SumEnergySq, false);

  // arrays of accumulated sums for energy and energy-squared have to be retrieved by function calls
  _Annealer.CalculateLogEvidence(getSumEnergy(), getSumEnergySq(), options.getNumAnnealedRuns());

  Finish(options, data, Log);
}

void AdmixMapModel::Iterate(const double* Coolnesses, unsigned coolness,
			    Options & options, InputData & data, LogWriter& Log,
			    double & SumEnergy, double & SumEnergySq,  bool AnnealedRun) {

  const int samples = options.getTotalSamples();
  const int burnin  = options.getBurnIn();

  //Accumulate Energy
  double AISz = 0.0;
  if(!AnnealedRun) cout << endl;

  for( int iteration = 0; iteration <= samples; iteration++ ) {

  #if DEBUG_ITER_TIMES
      genepi::CodeTimer ct;
      fprintf( stderr, "\nIteration-%d: begin %s\n", iteration, ct.local_started().c_str() );
  #endif

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
    else if ( options.getDisplayLevel() > 3 )
      {
      putc( '.', stderr );
      fflush( stderr );
      }

    //Sample Parameters
    UpdateParameters(iteration, options, Log, data.GetHiddenStateLabels(),
                     Coolnesses, Coolnesses[coolness], AnnealedRun);
    SubIterate(iteration, options, data, Log, SumEnergy, SumEnergySq,
	       AnnealedRun);
    #if DEBUG_ITER_TIMES
	fprintf( stderr, "Iteration-%d: end %s\n", iteration, ct.local_elapsed().c_str() );
    #endif
  }// end loop over iterations
  //use Annealed Importance Sampling to calculate marginal likelihood
  if(coolness>0) AISsumlogz += log(AISz /= (double)(samples-burnin));

}

void AdmixMapModel::UpdateParameters(int iteration, const Options& _options, LogWriter& Log,
                                     const Vector_s& PopulationLabels,
                                     const double* Coolnesses, double coolness,
                                     bool anneal) {

  const AdmixOptions & options = dynamic_cast<const AdmixOptions &>( _options );
  const int Populations    = options.getPopulations();
  const int burnin         = options.getBurnIn();
  const bool ChibIndicator = options.getChibIndicator();
  const bool afterBurnIn   = iteration > burnin;

  // if annealed run, anneal genotype probs - for testindiv only if testsingleindiv indicator set in IC
  if ( anneal || options.getTestOneIndivIndicator() )
    AdmixedIndividuals->annealGenotypeProbs( Loci.GetNumberOfChromosomes(), coolness, Coolnesses );

  A->ResetAlleleCounts(Populations);

  // ** update global sumintensities conditional on genotype probs and individual admixture proportions
  if ( Populations > 1 && options.getIndAdmixHierIndicator() &&
       (Loci.GetLengthOfGenome() + Loci.GetLengthOfXchrm() > 0.0) )
    L->UpdateGlobalSumIntensities(*AdmixedIndividuals,
                                  (!anneal && afterBurnIn && Populations > 1));
  // leaves individuals with HMM probs bad, stored likelihood ok
  // this function also sets locus correlations in Chromosomes

  // update the odds ratios vector psi
  if ( Loci.isX_data() && Populations > 1 )
    L->UpdateOddsRatios(*AdmixedIndividuals, afterBurnIn);

  //find posterior modes of individual admixture at end of burn-in
  //set Chib numerator
  if ( !anneal && iteration == burnin &&
       (ChibIndicator || strlen(options.getIndAdmixModeFilename())) ) {
    AdmixedIndividuals->FindPosteriorModes(options, R, L->getalpha(),
                                           L->getrhoalpha(), L->getrhobeta(),
                                           A, PopulationLabels);
    if ( ChibIndicator )
      AdmixedIndividuals->setChibNumerator(options, L->getalpha(),
                                           L->getrhoalpha(), L->getrhobeta(), A);
  }

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
  if ( options.hasAnyAssociationTests() || options.getTestForAffectedsOnly() )
    AdmixedIndividuals->SampleHapPairs(options, A, &Loci, true, anneal, true);

  if ( !anneal && afterBurnIn ) {
    //score tests
    if ( options.getScoreTestIndicator() )
      //TODO: broadcast dispersion in linear regression
      Scoretests.Update(R);//score tests evaluated for first outcome var only

    if(options.getTestForResidualAllelicAssoc()){
      ResidualAllelicAssocScoreTest.Reset();
      ResidualAllelicAssocScoreTest.Update(A->GetAlleleFreqs(), options.getHapMixModelIndicator());
    }
  }

  // sample individual admixture and sum-intensities
  AdmixedIndividuals->SampleParameters(iteration, options, R, L->getpoptheta(), L->getalpha(),
		       L->getrhoalpha(), L->getrhobeta(), Scoretests.getAncestryAssocTest(),anneal);
  // stored HMM likelihoods will now be bad if the sum-intensities are set at individual level

  if ( ChibIndicator && !anneal && iteration >= burnin )
    AdmixedIndividuals->updateChib(options, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), A);

  // update allele frequencies conditional on locus ancestry states
  // TODO: this requires fixing to anneal allele freqs for historicallelefreq model
  if ( !options.getFixedAlleleFreqs() )
    A->Update(IC, !anneal && afterBurnIn, coolness);

  // even for fixed allele freqs, must reset annealed genotype probs as unnannealed
  if(!options.getFixedAlleleFreqs() || anneal) {
    AdmixedIndividuals->setGenotypeProbs(&Loci); // sets unannealed probs ready for getEnergy
    AdmixedIndividuals->HMMIsBad(true);
  } // update of allele freqs sets HMM probs and stored loglikelihoods as bad

  // next update of stored loglikelihoods will be from getEnergy if not annealing run, from updateRhowithRW if globalrho,
  // or from update of individual-level parameters otherwise

  //update population admixture Dirichlet parameters conditional on individual admixture
  L->UpdatePopAdmixParams(iteration, AdmixedIndividuals, Log);

  // ** update regression parameters (if regression model) conditional on individual admixture
  bool condition = (!anneal && afterBurnIn);
  for(unsigned r = 0; r < R.size(); ++r){
    R[r]->Update(condition, IC->getOutcome(r), coolness );
    //output expected values of outcome variables to file every 'every' iterations after burnin
    if(condition && !(iteration % options.getSampleEvery()) ) {
      R[r]->OutputExpectedY();
    }
  }
}

void AdmixMapModel::ResetStepSizeApproximators(int resetk){
  Model::ResetStepSizeApproximators(resetk);
 L->resetStepSizeApproximator(resetk);
}

void AdmixMapModel::SubIterate(int iteration, Options& _options,
                               InputData& data, LogWriter& Log,
                               double& SumEnergy, double& SumEnergySq, bool AnnealedRun) {

  //cast Options object to AdmixOptions for access to ADMIXMAP options
  AdmixOptions& options = (AdmixOptions&) _options;

  const int burnin = options.getBurnIn();

    //Output parameters to file and to screen
    Log.setDisplayMode(Quiet);
     if(!AnnealedRun){
      // output every 'getSampleEvery()' iterations
      if(!(iteration % options.getSampleEvery()))
	OutputParameters(iteration, &options, Log);

      // ** set merged haplotypes for allelic association score test
      if ( iteration == burnin ) {

	if(options.getTestForHaplotypeAssociation())
	  Scoretests.MergeRareHaplotypes(L->getalpha0().getVector_unsafe());

	if( options.getStratificationTest() )
	  StratTest.Initialize( &options, Loci, IC, Log);
      }

      //Updates and Output after BurnIn
      if( !AnnealedRun && iteration > burnin){
	//dispersion test
	if( options.getTestForDispersion() )DispTest.TestForDivergentAlleleFrequencies(A, IC);
	//stratification test
	if( options.getStratificationTest() )StratTest.calculate(IC, A->GetAlleleFreqs(), Loci.GetChrmAndLocus(),
								 options.getPopulations());
	//tests for mis-specified allelefreqs
	if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
	  AlleleFreqTest.Update(IC, A, &Loci);
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
	      A->OutputErgodicAvg(samples, &avgstream);//dispersion parameter in dispersion model, freq Dirichlet param prior rate in hapmixmodel
	    }
	    for(unsigned r = 0; r < R.size(); ++r)//regression params
	      R[r]->OutputErgodicAvg(samples,avgstream);

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
  bclib::Delimitedstdout ScreenWriter(' ');

  if ( options->getDisplayLevel() > 3 )
    cout << ' ';

  if(options->getIndAdmixHierIndicator()  ){
    //output population-level parameters only when there is a hierarchical model on indadmixture
    // ** pop admixture, sumintensities
    if(options->getPopulations() > 1){
      //write to file
      if( iteration > options->getBurnIn() ){
	L->OutputParams();
	//write to R object
	L->OutputParams(paramstream);
      }
      //write to screen
      if(options->getDisplayLevel()>2)
	L->OutputParams(ScreenWriter);
    }

    // ** dispersion parameter (if dispersion model)
    if( iteration > options->getBurnIn() ){
      //write to file
      A->OutputParams();
      //write to R object
      A->OutputParams(paramstream);
    }
    //write to screen
    if(options->getDisplayLevel()>2)
      A->OutputParams(ScreenWriter);
  }

  // ** regression parameters
  for(unsigned r = 0; r < R.size(); ++r){
    //output regression parameters to file if after burnin and to screen if displaylevel >2
    if( iteration > options->getBurnIn() ){
      //write to paramfile
      R[r]->OutputParams(options->getNumberOfOutcomes());
      //write to R object
      R[r]->OutputParams(paramstream);
    }
    //write to screen
    if(options->getDisplayLevel()>2){
      if(R.size()>1)
	cout << "\nRegression " << r +1 << " ";
      R[r]->OutputParams(ScreenWriter);
    }
  }

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
	A->OutputAlleleFreqs();
    }
  }
  // cout << endl;
  if( iteration > options->getBurnIn() )
    paramstream << bclib::newline;
}

void AdmixMapModel::PrintAcceptanceRates(const Options& _options, LogWriter& Log){
  const AdmixOptions& options = (const AdmixOptions&)_options;
 if( options.getIndAdmixHierIndicator() ){

   if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
   else Log.setDisplayMode(On);

   if(options.getPopulations() > 1 ){
     L->printAcceptanceRates(Log);
   }

   A->PrintAcceptanceRates(Log);

   A->OutputAlleleFreqSamplerAcceptanceRates(options.getResultsDir());
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
    if (Loci.isX_data() && options.isGlobalPsi())
      L->StoreOddsRatiosPosteriorMean(*AdmixedIndividuals);
    AdmixedIndividuals->WritePosteriorMeans(options, data.GetHiddenStateLabels(), &Loci);
  }
  //FST
//   if( strlen( options.getHistoricalAlleleFreqFilename()) ){
//     A->OutputFST();
//   }
  A->CloseOutputFile((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery(),
		    data.GetHiddenStateLabels());

  if(options.getOutputAlleleFreq()){
//     ofstream klinfofile((options.getResultsDir() + "/KLinfo.txt").c_str());
//     A->WriteKLInfo((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery(), klinfofile);
//     klinfofile.close();

    A->WriteLocusInfo((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery(), options.getResultsDir(), data.GetHiddenStateLabels());
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
    Scoretests.OutputLikelihoodRatios(data.GetHiddenStateLabels());

  //print deviance at posterior mean, DIC
  _Annealer.PrintResults(Log, getDevianceAtPosteriorMean(options, Log));

  //print results of annealed importance sampling
  if(options.getThermoIndicator()){
    Log << "\nAnnealed Importance Sampling estimates log marginal likelihood as " << AISsumlogz << "\n";
  }

  //if(options.outputParams())
  {
    WriteParamsAsRObjectDimensions(options, data);
  }

}

void AdmixMapModel::WriteParamsAsRObjectDimensions(const AdmixOptions& options, const InputData& data){
  //write dimensions and dimnames of R object holding sampled parameters

  const Vector_s& PopulationLabels = data.GetHiddenStateLabels();

  //population admixture and sumintensities
  vector<vector<string> > dimnames(1);
  //if(strlen(options.getParameterFilename()))
  {
    for( int i = 0; i < options.getPopulations(); i++ ) {
      stringstream ss;
      ss << "Dirichlet." << PopulationLabels[i];
      dimnames[0].push_back(ss.str());
    }
    //SumIntensities
    if( options.isGlobalRho() )  dimnames[0].push_back("sumIntensities");
    else  dimnames[0].push_back("sumIntensities.mean");

    // Odds ratios for the X chromosome
    if (Loci.isX_data()) {
      for (int i = 0; i < options.getPopulations(); i++) {
        stringstream ss;
        ss << "Psi." << PopulationLabels[i];
        dimnames[0].push_back(ss.str());
      }
    }
  }

  //allele freq precision
  //if(strlen(options.getEtaOutputFilename()))
  {
    //dispersion model
    if(strlen( options.getHistoricalAlleleFreqFilename() )){
      for( int k = 0; k < options.getPopulations(); k++ ){
	stringstream ss;
	ss << "eta." << PopulationLabels[k];
	dimnames[0].push_back(ss.str());
      }
    }
    //correlated freqs model
    else if(options.getCorrelatedAlleleFreqs()){
      dimnames[0].push_back("eta");
    }
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
  dims.push_back((options.getTotalSamples() - options.getBurnIn()) / options.getSampleEvery());

  paramstream.close(dims, dimnames);
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
  //InitializeErgodicAvgFile(&options, Log, data.GetHiddenStateLabels(), data.getCovariateLabels());


}
//this function is here because three different objects have to write to avgstream
void AdmixMapModel::InitializeErgodicAvgFile(const Options& _options, LogWriter &Log,
				     const Vector_s& PopulationLabels, const Vector_s& CovariateLabels){
  AdmixOptions& options = (AdmixOptions&) _options;

  Log.setDisplayMode(bclib::Quiet);
  //Open ErgodicAverageFile
  if( strlen( options.getErgodicAverageFilename() ) ) {
    avgstream.open( options.getErgodicAverageFilename(), ios::out );
    if( ! avgstream ){
      Log.setDisplayMode(bclib::On);
      Log << "ERROR: Couldn't open Ergodic Average file\n";
      //exit( 1 );
    } else {
      Log << "Writing ergodic averages of parameters to "
	  << options.getErgodicAverageFilename() << "\n\n";
    }

    // Header line of ergodicaveragefile
    if( options.getIndAdmixHierIndicator() ){
      if(options.getPopulations()>1){
	  for( int i = 0; i < options.getPopulations(); i++ ){
	    avgstream << PopulationLabels[i] << "\t";
	  }
	if( options.isGlobalRho() ) avgstream << "sumIntensities\t";
	else avgstream << "sumIntensities.mean\t";
      }

      // dispersion parameters
      if( strlen( options.getHistoricalAlleleFreqFilename() ) ){
	for( int k = 0; k < options.getPopulations(); k++ ){
	  avgstream << "eta" << k << "\t";
	}
      }
    }//end if hierarchical model

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
    if(options.getChibIndicator()){// chib calculation
      avgstream << "LogPrior\tLogPosterior\tLogPosteriorAdmixture\tLogPosteriorSumIntensities\t"
		 << "LogPosteriorAlleleFreqs\tLogMarginalLikelihood";
    }
    avgstream << "\n";
  } else {
    Log << "Not writing ergodic averages to file\n";
  }
}


double AdmixMapModel::getDevianceAtPosteriorMean(const Options& options, LogWriter& Log){
    return AdmixedIndividuals->getDevianceAtPosteriorMean(options, R, &Loci, Log, L->getSumLogRho(), Loci.GetNumberOfChromosomes(), A);
}

double* AdmixMapModel::getSumEnergy()const{
    return AdmixedIndividuals->getSumEnergy();
}

double* AdmixMapModel::getSumEnergySq()const{
    return AdmixedIndividuals->getSumEnergySq();
}
