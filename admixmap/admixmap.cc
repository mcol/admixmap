/** 
 *   ADMIXMAP
 *   admixmap.cc 
 *   Top-level source file
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "admixmap.h"
#include <fstream>
#include <dirent.h>//for OpenResultsDir
#include "Comms.h"

#define ADMIXMAP_VERSION 3.5
#define MAXNUMOPTIONS 50//maximum number of options specifiable.
#define MAXOPTIONLENGTH 1024//maximum number of characters in an option line (excluding spaces)

using namespace std;
double coolness = 1.0; // default

void MakeResultsDir(const char* dirname, bool verbose);
void InitializeErgodicAvgFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
			      LogWriter &Log, std::ofstream *avgstream, 
			      const Vector_s& PopLabels, const Vector_s& CovariateLabels);
void UpdateParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, vector<Regression *>&R, 
		      const AdmixOptions *options, const Genome *Loci, ScoreTests *S, LogWriter& Log, 
		      const Vector_s& PopulationLabels, double coolness, bool anneal);
void OutputParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, vector<Regression *>&R, 
		      const AdmixOptions *options, LogWriter& Log);
void WriteIterationNumber(const int iteration, const int width, int displayLevel);

void PrintCopyrightNotice(LogWriter & Log); 

void doIterations(const int & samples, const int & burnin, IndividualCollection *IC, Latent & L, AlleleFreqs  & A, 
		  vector<Regression *>&R, 
		  AdmixOptions & options,  const Genome  & Loci, LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
		  const double coolness, bool AnnealedRun, ofstream & loglikelihoodfile, 
		  ScoreTests & Scoretests, DispersionTest & DispTest, StratificationTest & StratTest, 
		  MisSpecAlleleFreqTest & AlleleFreqTest, HWTest & HWtest, ofstream & avgstream, InputData & data, 
		  const double* Coolnesses);

void OutputErgodicAvgDeviance(int samples, double & SumEnergy, double & SumEnergySq, std::ofstream *avgstream);

void PrintOptionsMessage();
void ThrowException(const string& msg, LogWriter & Log);

int main( int argc , char** argv ){
  if (argc < 2 || (argc ==2 && ( !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) ) ) {
    PrintOptionsMessage();
    exit(1); 
  } 
#ifdef PARALLEL
  MPI::Init(argc, argv);
  Comms::Initialise();

  double t1 = 0, t2 = 0, t3 = 0, t0 = MPI::Wtime();
  //define states for logging
  MPE_Init_log();
  if(Comms::isMaster()){
    MPE_Describe_state(1, 2, "Barrier", "red:vlines1");
    MPE_Describe_state(3, 4, "ReduceAncestry", "cyan:vlines3");
    MPE_Describe_state(5, 6, "ReduceAlleleCounts", "orange:gray2");
    MPE_Describe_state(7, 8, "SampleAlleleFreqs", "yellow:gray3");
    MPE_Describe_state(9, 10, "SampleRho", "LightBlue:gray3");
    MPE_Describe_state(11, 12, "SetGenotypeProbs", "LightGreen:hlines2");
    MPE_Describe_state(13, 14, "SampleHapPairs", "magenta:vlines2");
    MPE_Describe_state(15, 16, "SampleAncestry", "blue:vlines3");
    MPE_Describe_state(17, 18, "Bcastrho", "maroon:gray");
    MPE_Describe_state(19, 20, "BcastFreqs", "plum:hlines3");
    MPE_Describe_state(21, 22, "ScoreTests", "gray:hlines3");
  }
#endif

  // ******************* PRIMARY INITIALIZATION ********************************************************************************
  //read user options
  AdmixOptions options(argc, argv);

  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  if(isMaster){
    MakeResultsDir(options.getResultsDir().c_str(), (options.getDisplayLevel()>2));
  }
 
  //open logfile, start timer and print start message
  LogWriter Log(options.getLogFilename(), (bool)(options.getDisplayLevel()>1 && isMaster));
  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
  if(isMaster){
    //if(options.getDisplayLevel()>0 )
    PrintCopyrightNotice(Log);
    Log.StartMessage();
  }

  try{  
    Rand RNG;//allocate random number generator
    RNG.setSeed( options.getSeed() );  // set random number seed
  
    InputData data; //read data files and check (except allelefreq files)
    data.readData(&options, Log);//also sets 'numberofoutcomes' and 'populations' options

     //check user options
    if(options.checkOptions(Log, data.getNumberOfIndividuals())){
      cout << endl << endl;
      exit(1);
    }
  
    //print user options to args.txt; must be done after all options are set
    if(isMaster)options.PrintOptions();

    Genome Loci;
    Loci.Initialise(&data, options.getPopulations(), Log);//reads locusfile and creates CompositeLocus objects
    if(isFreqSampler){
      //print table of loci for R script to read
      string locustable = options.getResultsDir();
      locustable.append("/LocusTable.txt");
      Loci.PrintLocusTable(locustable.c_str(), data.getLocusMatrix().getCol(1));
      locustable.clear();
    }

    AlleleFreqs A(&Loci);
    if(isFreqSampler)//allele freq updater only
      A.Initialise(&options, &data, Log); //checks allelefreq files, initialises allele freqs and finishes setting up Composite Loci
    if(isFreqSampler || isWorker)A.AllocateAlleleCountArrays(options.getPopulations());
#ifdef PARALLEL
    //broadcast initial values of freqs
    if(!isMaster)A.BroadcastAlleleFreqs();
#endif

    IndividualCollection *IC = new IndividualCollection(&options, &data, &Loci);//NB call after A Initialise
    if(isMaster || isWorker)IC->LoadData(&options, &data);                             //and before L and R Initialise
    if(isWorker)IC->setGenotypeProbs(&Loci, &A); // sets unannealed probs

    Latent L( &options, &Loci);    
    if(isMaster || isWorker)L.Initialise(IC->getSize(), data.GetPopLabels(), Log);

    vector<Regression*> R;//vector of regression pointers
    if (options.getNumberOfOutcomes()>0 && (isMaster || isWorker)){
      if( isMaster ){
	Regression::OpenOutputFile(options.getNumberOfOutcomes(), options.getRegressionOutputFilename(), Log);  
	Regression::OpenExpectedYFile(options.getEYFilename(), Log);
      }
      for(int r = 0; r < options.getNumberOfOutcomes(); ++r){
	//determine regression type and allocate regression objects
	if( data.getOutcomeType(r)== Binary ) R.push_back( new LogisticRegression() );
	else if( data.getOutcomeType(r)== Continuous ) R.push_back( new LinearRegression());
	else if( data.getOutcomeType(r)== CoxData ) R.push_back(new CoxRegression());

	if(isMaster) {
	  if(R[r]->getRegressionType()==Cox)
	    R[r]->Initialise(r, options.getRegressionPriorPrecision(), IC->getCovariatesMatrix(),data.getCoxOutcomeVarMatrix(), Log);
	  else
	    R[r]->Initialise(r, options.getRegressionPriorPrecision(), IC->getCovariatesMatrix(), IC->getOutcomeMatrix(), Log);
	}
	else R[r]->Initialise(r, IC->GetNumCovariates());
	R[r]->InitializeOutputFile(data.getCovariateLabels(), options.getNumberOfOutcomes());
      }
    }
    
    if( options.getPopulations()>1 && (options.isGlobalRho() || options.getHapMixModelIndicator()) && isWorker) {
      Loci.SetLocusCorrelation(L.getrho());
      if(options.getHapMixModelIndicator())
	for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ ){
	  //set global state arrival probs in hapmixmodel
	  //TODO: can skip this if xonly analysis with no females
	  //NB: assumes always diploid in hapmixmodel
	  Loci.getChromosome(j)->SetHMMTheta(L.getGlobalTheta(), options.isRandomMatingModel(), true);
	  Loci.getChromosome(j)->SetStateArrivalProbs(options.isRandomMatingModel());
	}
    }
    
    if(isMaster || isWorker)
      IC->Initialise(&options, &Loci, data.GetPopLabels(), L.getalpha(), /*L.getrhoalpha(), L.getrhobeta(),*/ Log);
  
    data.Delete();

    //  ******** single individual, one population, fixed allele frequencies  ***************************
    if( IC->getSize() == 1 && options.getPopulations() == 1 && strlen(options.getAlleleFreqFilename()) )
      // nothing to do except calculate likelihood
      IC->getOnePopOneIndLogLikelihood(Log, data.GetPopLabels());
    else {
      // ******************* INITIALIZE TEST OBJECTS and ergodicaveragefile *******************************
      DispersionTest DispTest;
      StratificationTest StratTest(options.getStratTestFilename(), Log);
      ScoreTests Scoretests;
      MisSpecAlleleFreqTest AlleleFreqTest;
      HWTest HWtest;
      std::ofstream avgstream; //output to ErgodicAverageFile
      if( options.getScoreTestIndicator() && (isMaster || isWorker)){
	Scoretests.Initialise(&options, IC, &Loci, data.GetPopLabels(), Log);
      }
      //if(isMaster){
      if( options.getTestForDispersion() ){
	DispTest.Initialise(&options, Log, Loci.GetNumberOfCompositeLoci());    
      }
      //if( options.getStratificationTest() )
      //StratTest.Initialize( &options, Loci, IC, Log);
      if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
	AlleleFreqTest.Initialise(&options, &Loci, Log );  
      if( options.getHWTestIndicator() )
	HWtest.Initialise(&options, Loci.GetTotalNumberOfLoci(), Log);
      if(isMaster)
	InitializeErgodicAvgFile(&options, IC, Log, &avgstream, data.GetPopLabels(), data.getCovariateLabels());
      //}

      string s = options.getResultsDir()+"/loglikelihoodfile.txt";
      ofstream loglikelihoodfile(s.c_str());
    
      // ******************* initialise stuff for annealing ************************************************
      double IntervalRatio = 1.03; // size of increments of coolness increases geometrically
      int NumAnnealedRuns = options.getNumAnnealedRuns(); // number of annealed runs excluding last run at coolness of 1
      // set number of samples - 1 for annealing runs, use "samples" options otherwise. Overriden for final, unannealed run with "thermo" option
      int samples = options.getThermoIndicator()? 1 : options.getTotalSamples(); 
      int burnin = 1; // for annealing runs we only want burnin
    
      double SumEnergy = 0.0, SumEnergySq = 0.0, LogEvidence = 0.0;
      double MeanEnergy = 0.0, VarEnergy = 0.0;
      double LastMeanEnergy = 0.0;
    
      coolness = 1.0; // default
      bool AnnealedRun = false;
      std::ofstream annealstream;//for monitoring energy when annealing
    
      // set annealing schedule
      double *IntervalWidths = 0;
      double *Coolnesses = 0; // 
	IntervalWidths = new double[NumAnnealedRuns];
	Coolnesses = new double[NumAnnealedRuns + 1];
	Coolnesses[0] = 1.0;//default for unnannealed run

      if(NumAnnealedRuns > 0) {
	Coolnesses[0] = 0.0; // change this if you want annealing to start somewhere other than 0;
	// set initial increment of coolness so that geometric series of NumAnnealedRuns increments 
	// will sum to 1 - Coolnesses[0] after NumAnnealedRuns + 1 terms
	IntervalWidths[0] = (1.0 - Coolnesses[0]) * (1.0 - IntervalRatio) /(1.0 - pow(IntervalRatio, NumAnnealedRuns)); 
	Coolnesses[1] = Coolnesses[0] + IntervalWidths[0];
	if(NumAnnealedRuns > 1) {
	  for(int run=2; run < NumAnnealedRuns; ++run) {
	    IntervalWidths[run - 1] = IntervalWidths[run - 2] * IntervalRatio; // geometric increments in interval width
	    Coolnesses[run] = Coolnesses[run - 1] + IntervalWidths[run - 1];  
	  }
	}
	Coolnesses[NumAnnealedRuns] = 1.0;
      }
    
      if(options.getThermoIndicator()) { // set up output for thermodynamic integration
	if(isMaster){		
	  string s = options.getResultsDir()+"/annealmon.txt";
	  annealstream.open(s.c_str());
	  annealstream << "Coolness\tMeanEnergy\tVarEnergy\tlogEvidence" << endl;
	}
	samples = options.getTotalSamples();
	burnin = options.getBurnIn();
      }

      if( options.getTestOneIndivIndicator() )NumAnnealedRuns = 0;
      if(isMaster){
	if(NumAnnealedRuns > 0) {
	  Log << On << NumAnnealedRuns << " annealing runs of " << samples 
	      << " iteration(s) followed by final run of "; 
	}
	if(!options.getThermoIndicator()) {
	  Log << options.getTotalSamples();
	} else { 
	  Log << 2 * options.getTotalSamples();
	} 
	Log << " iterations at ";
	if( options.getTestOneIndivIndicator() ) {
	  Log << options.getNumAnnealedRuns()+1 
	      <<" coolnesses for test individual. Other individuals at ";
	}
	Log << "coolness of 1\n";
      }
    
      //Write initial values
      //     if(options.getIndAdmixHierIndicator()  ){
      //       if(options.getDisplayLevel()>2)Log.setDisplayMode(On);
      //       else Log.setDisplayMode(Quiet);
      //       //Log << "InitialParameterValues:\n"
      //       //OutputParameters(-1, IC, &L, &A, R, &options, Log);
      //       //Log << "\n";
      //     }
#ifdef PARALLEL
      t1 = MPI::Wtime()-t0;
#endif
      if(!options.getTestOneIndivIndicator()) {  
	for(int run=0; run <= NumAnnealedRuns; ++run) { //loop over coolnesses from 0 to 1
	  // should call a posterior mode-finding algorithm before last run at coolness of 1
	  //resets for start of each run
	  SumEnergy = 0.0;//cumulative sum of modified loglikelihood
	  SumEnergySq = 0.0;//cumulative sum of square of modified loglikelihood

	  if(run == NumAnnealedRuns) {
	    AnnealedRun = false;
	    coolness = 1.0;
	    if(options.getThermoIndicator()) {
	      samples *= 2 ; // last run is longer
	    }
	    burnin = options.getBurnIn();
	  } 
	  else {
	    AnnealedRun = true; 
	    coolness = Coolnesses[run];
	  }
	  if(NumAnnealedRuns > 0) {
	    cout <<"\rSampling at coolness of " << coolness << "         " << flush;
	    // reset approximation series in step size tuners
	    int resetk = NumAnnealedRuns; //   
	    if(samples < NumAnnealedRuns) {// samples=1 if annealing without thermo integration
	      resetk = 1 + run;
	    }
	    IC->resetStepSizeApproximators(resetk); 
	    A.resetStepSizeApproximator(resetk);
	    L.resetStepSizeApproximator(resetk);
	    
	  }
	  // accumulate scalars SumEnergy and SumEnergySq at this coolness
	  // array Coolnesses is not used unless TestOneIndivIndicator is true
	  doIterations(samples, burnin, IC, L, A, R, options, Loci, Log, SumEnergy, SumEnergySq, coolness, 
		       AnnealedRun, loglikelihoodfile, Scoretests, DispTest, StratTest, AlleleFreqTest, 
		       HWtest, avgstream, data, Coolnesses);
#ifdef PARALLEL
	  t2 = MPI::Wtime()-t1;
#endif
	  if(isMaster){	  
	    //calculate mean and variance of energy at this coolness
	    MeanEnergy = SumEnergy / ((double)samples - options.getBurnIn());
	    VarEnergy  = SumEnergySq / ((double)samples - options.getBurnIn()) - MeanEnergy * MeanEnergy;
	    if(options.getThermoIndicator()){// calculate thermodynamic integral
	      annealstream << coolness << "\t" << MeanEnergy << "\t" << VarEnergy;
	      if(run > 0) { // use trapezium rule to approximate integral
		LogEvidence -= 0.5*(LastMeanEnergy + MeanEnergy) * IntervalWidths[run];
	      } 
	      annealstream <<"\t"<< LogEvidence << endl; 
	      LastMeanEnergy = MeanEnergy;
	    }
	  } 
	} // end loop over coolnesses
      } 
      else { // evaluate energy for test individual only at all coolnesses simultaneously
	// call with argument AnnealedRun false - copies of test individual will be annealed anyway  
	doIterations(samples, burnin, IC, L, A, R, options, Loci, Log, SumEnergy, SumEnergySq, 1.0, false, 
		     loglikelihoodfile, Scoretests, DispTest, StratTest, AlleleFreqTest, HWtest, avgstream, data, Coolnesses);
	// arrays of accumulated sums for energy and energy-squared have to be retrieved by function calls
	double *MeanEner = IC->getSumEnergy(); 
	double *VarEner = IC->getSumEnergySq();
	
	double LastMeanEnergy = 0.0; 
	for(int ii = 0; ii < options.getNumAnnealedRuns()+1; ++ii) { // loop over coolnesses to evaluate integral
	  //calculate mean and variance of energy at each coolness
	  MeanEner[ii] /=  ((double)options.getTotalSamples() - options.getBurnIn());
	  VarEner[ii] = VarEner[ii] /  ((double)options.getTotalSamples() - options.getBurnIn()) - MeanEner[ii]*MeanEner[ii];
	  annealstream << Coolnesses[ii] << "\t" << MeanEner[ii] << "\t" << VarEner[ii];
	  // use trapezium rule to approximate integral
	  LogEvidence -= 0.5*(LastMeanEnergy + MeanEner[ii]) * IntervalWidths[ii]; 
	  annealstream <<"\t"<< LogEvidence << endl; 
	  LastMeanEnergy = MeanEner[ii];
	}
	MeanEnergy = MeanEner[options.getNumAnnealedRuns()];//mean at coolness of 1;
	VarEnergy = VarEner[options.getNumAnnealedRuns()];//var at   ""
      } // end evaluation of test individual

  
      delete[] IntervalWidths;
      delete[] Coolnesses;
      if(isMaster)
	cout<< "\nIterations completed                       \n" << flush;
    
      // *************************** OUTPUT AT END ***********************************************************

      if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);	
      else Log.setDisplayMode(Quiet);
      if( options.getChibIndicator()) {
	//IC->OutputChibEstimates(options.isRandomMatingModel(), Log, options.getPopulations());
	//MLEs of admixture & sumintensities used in Chib algorithm to estimate marginal likelihood
	if(IC->getSize()==1) IC->OutputChibResults(Log);
      }
		
      double Information = -LogEvidence - MeanEnergy;
      double MeanDeviance = 2.0 * MeanEnergy; 
      double VarDeviance = 4.0 * VarEnergy;
      Log << Quiet << "\nMeanDeviance(D_bar)\t" << MeanDeviance << "\n"
	  << "VarDeviance(V)\t" << VarDeviance << "\n"
	  << "PritchardStat(D_bar+0.25V)\t" << MeanDeviance + 0.25*VarDeviance << "\n";
      double D_hat = IC->getDevianceAtPosteriorMean(&options, R, &Loci, Log, L.getSumLogRho(), Loci.GetNumberOfChromosomes(), &A);
      double pD = MeanDeviance - D_hat;
      double DIC = MeanDeviance + pD;
      Log << Quiet << "DevianceAtPosteriorMean(D_hat)\t" << D_hat << "\n"
	  << "EffectiveNumParameters(pD)\t" << pD << "\n"
	  << "DevianceInformationCriterion\t" << DIC << "\n\n"; 
		
      if(options.getThermoIndicator()){
	Log << "thermodynamic integration for marginal likelihood yields:\n";
	Log << "LogEvidence " <<  LogEvidence << "\n"; 
	Log << "Information (negative entropy, measured in nats) " << Information << "\n";
      }
		
      //Expected Outcome
      if(isMaster && options.getNumberOfOutcomes() > 0){
	Regression::FinishWritingEYAsRObject((options.getTotalSamples()-options.getBurnIn())/ options.getSampleEvery(), 
					     data.getOutcomeLabels());
      }
      //FST
      if( strlen( options.getHistoricalAlleleFreqFilename()) && isFreqSampler  ){
	A.OutputFST();
      }
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

      if( options.getScoreTestIndicator() && isMaster ) {
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
      if(annealstream.is_open())annealstream.close();
      if(avgstream.is_open())avgstream.close();

    }//end else

    delete IC;//must call explicitly so IndAdmixOutputter destructor finishes writing to indadmixture.txt
    for(unsigned r = 0; r < R.size(); ++r)delete R[r];
    R.clear();

    if(isFreqSampler && !options.getHapMixModelIndicator())
      A.CloseOutputFile((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery(), data.GetPopLabels());
    if(isMaster)cout << "Output to files completed\n" << flush;
    
    // ******************* acceptance rates - output to screen and log ***************************
    if( options.getIndAdmixHierIndicator() ){
      if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
      else Log.setDisplayMode(On);
      if(options.getPopulations() > 1 && isMaster){
	L.printAcceptanceRates(Log);
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
    }
    A.OutputAlleleFreqSamplerAcceptanceRates((options.getResultsDir() + "/AlleleFreqSamplerAcceptanceRates.txt").c_str());
    if(options.getHapMixModelIndicator() && !options.getFixedAlleleFreqs() && isFreqSampler){
      Log << "Average expected Acceptance rate in allele frequency prior parameter sampler:\n" << A.getHapMixPriorSamplerAcceptanceRate()
	  << "\nwith average final step size of " << A.getHapMixPriorSamplerStepSize() << "\n";
    }
    if(isMaster){
      if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
      Log.ProcessingTime();
    }

#ifdef PARALLEL
    t3 = MPI::Wtime()-t2;
#endif
    //if(options.getDisplayLevel()>0){
#ifdef PARALLEL
    cout << "Rank " << MPI::COMM_WORLD.Get_rank() << " finished.\n Initialization: " << t1 <<"s Main loop: "<< t2 << "s Finalization:" << t3 << "s.\n";
#else
    cout << "Finished" << endl;
#endif

  } 
  catch (const string& msg) {//catch any stray error messages thrown upwards
    ThrowException(msg, Log);
  }
  catch (const char* msg) {//in case error messages thrown as char arrays instead of strings
    ThrowException(string(msg), Log);
  }
  catch (exception& e){
    ThrowException(e.what(), Log);
  }

#ifdef PARALLEL
  catch(MPI::Exception e){
    cout << "Error in process " << MPI::COMM_WORLD.Get_rank() << ": " << e.Get_error_code() << endl;
    MPI::COMM_WORLD.Abort(1);
  }
#endif
  catch(...){
    cout << "Unknown exception occurred. Contact the program authors for assistance" << endl;
    exit(1);
  }
#ifdef PARALLEL
  //MPI::COMM_WORLD.Barrier();
  Comms::Finalise();
  MPE_Finish_log("admixmap");
  MPI_Finalize();
#endif
  if(isMaster){//print line of *s
    cout <<setfill('*') << setw(80) << "*" <<endl;
  }
  putenv("ADMIXMAPCLEANEXIT=1");
  return 0;
} //end of main

void doIterations(const int & samples, const int & burnin, IndividualCollection *IC, Latent & L, AlleleFreqs & A, 
		  vector<Regression *>&R, AdmixOptions & options, 
		  const Genome & Loci, LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
		  double coolness, bool AnnealedRun, ofstream & loglikelihoodfile, 
		  ScoreTests & Scoretests, DispersionTest & DispTest, StratificationTest & StratTest, 
		  MisSpecAlleleFreqTest & AlleleFreqTest, HWTest & HWtest, ofstream & avgstream, InputData & data, const double* Coolnesses) {
  const bool isMaster = Comms::isMaster();
  //const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  double Energy = 0.0;
  if(isMaster && !AnnealedRun) cout << endl;
  for( int iteration = 0; iteration <= samples; iteration++ ) {
    if( (isMaster || isWorker) && iteration > burnin) {
      //accumulate energy as minus loglikelihood, calculated using unnanealed genotype probs
      if( !options.getTestOneIndivIndicator() ) {
	Energy = IC->getEnergy(&options, R, AnnealedRun); // should store loglikelihood if not AnnealedRun
	if(isMaster){
	  SumEnergy += Energy;
	  SumEnergySq += Energy*Energy;
	  // write to file if not AnnealedRun
	  if(!AnnealedRun)loglikelihoodfile << iteration<< "\t" << Energy <<endl;
	}
      } else {  
	IC->accumulateEnergyArrays(&options);
      }
    }
    if( isMaster && !AnnealedRun &&  !(iteration % options.getSampleEvery()) ) {
      WriteIterationNumber(iteration, (int)log10((double) samples+1 ), options.getDisplayLevel());
    }
    
    // if annealed run, anneal genotype probs - for testindiv only if testsingleindiv indicator set in IC
    if((isMaster || isWorker) && (AnnealedRun || options.getTestOneIndivIndicator())) 
      IC->annealGenotypeProbs(Loci.GetNumberOfChromosomes(), coolness, Coolnesses); 
    
    UpdateParameters(iteration, IC, &L, &A, R, &options, &Loci, &Scoretests, Log, data.GetPopLabels(), coolness, AnnealedRun);

    Log.setDisplayMode(Quiet);
 
    if(!AnnealedRun){    
      // output every 'getSampleEvery()' iterations
      if(!(iteration % options.getSampleEvery()) && isMaster)
	OutputParameters(iteration, IC, &L, &A, R, &options, Log);
      
      // ** set merged haplotypes for allelic association score test 
      if( (isMaster || isWorker) && iteration == options.getBurnIn() ){
#ifndef PARALLEL
	if(options.getTestForAllelicAssociation())
	  Scoretests.SetAllelicAssociationTest(L.getalpha0());
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
	      L.OutputErgodicAvg(samples,&avgstream);//pop admixture params, pop (mean) sumintensities
	      A.OutputErgodicAvg(samples, &avgstream);//dispersion parameter in dispersion model, freq Dirichlet param prior rate in hapmixmodel
	    }
	    for(unsigned r = 0; r < R.size(); ++r)//regression params
	      R[r]->OutputErgodicAvg(samples, &avgstream);

	    OutputErgodicAvgDeviance(samples, SumEnergy, SumEnergySq, &avgstream);
	    if(options.getChibIndicator()) IC->OutputErgodicChib(&avgstream, options.getFixedAlleleFreqs());
	    avgstream << endl;
	  }
	  //Score Test output
	  if( options.getScoreTestIndicator() )  Scoretests.Output(iteration-burnin, data.GetPopLabels(), data.getLocusLabels(), false);
	}//end "if every'*10" block
      }//end "if after BurnIn" block
    } // end "if not AnnealedRun" block
	
  }// end loop over iterations
}

//this function is here because three different objects have to write to avgstream
void InitializeErgodicAvgFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
			      LogWriter &Log, std::ofstream *avgstream, 
			      const Vector_s& PopulationLabels, const Vector_s& CovariateLabels){
  Log.setDisplayMode(Quiet);
  //Open ErgodicAverageFile  
  if( strlen( options->getErgodicAverageFilename() ) ) {
    avgstream->open( options->getErgodicAverageFilename(), ios::out );
    if( !*avgstream ){
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
	if(!options->getHapMixModelIndicator())
	  for( int i = 0; i < options->getPopulations(); i++ ){
	    *avgstream << PopulationLabels[i] << "\t";
	  }
	if( options->isGlobalRho() ) *avgstream << "sumIntensities\t";
	else *avgstream << "sumIntensities.mean\t";
      }
      
      // dispersion parameters
      if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
	for( int k = 0; k < options->getPopulations(); k++ ){
	  *avgstream << "eta" << k << "\t";
	}
      }
    }//end if hierarchical model

    //rate parameter of prior on frequency Dirichlet prior params
    if(options->getHapMixModelIndicator() && !options->getFixedAlleleFreqs()){
      *avgstream << "FreqPriorRate\t"; 
    }

    // Regression parameters
    if( options->getNumberOfOutcomes() > 0 ){
      for(int r = 0; r < options->getNumberOfOutcomes(); ++r){
	if(individuals->getOutcomeType(r)!=CoxData)
	*avgstream << "intercept\t";

	//write covariate labels to header
	  copy(CovariateLabels.begin(), CovariateLabels.end(), ostream_iterator<string>(*avgstream, "\t")); 

	if( individuals->getOutcomeType(r)==Continuous )//linear regression
	  *avgstream << "precision\t";
      }
    }

    *avgstream << "MeanDeviance\tVarDeviance\t";
    if(options->getChibIndicator()){// chib calculation
      *avgstream << "LogPrior\tLogPosterior\tLogPosteriorAdmixture\tLogPosteriorSumIntensities\t"
		 << "LogPosteriorAlleleFreqs\tLogMarginalLikelihood";
    }
    *avgstream << "\n";
  } else {
    Log << "No ergodicaveragefile given\n";
  }
}

void UpdateParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, vector<Regression *>&R, 
		      const AdmixOptions *options, const Genome *Loci, ScoreTests *S, LogWriter& Log, 
		      const Vector_s& PopulationLabels, double coolness, bool anneal){
  const bool isMaster = Comms::isMaster();
  const bool isFreqSampler = Comms::isFreqSampler();
  const bool isWorker = Comms::isWorker();

  A->ResetAlleleCounts(options->getPopulations());

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(isMaster || isWorker){
    // ** update global sumintensities conditional on genotype probs and individual admixture proportions
    if((options->getPopulations() > 1) && options->getIndAdmixHierIndicator() && !options->getHapMixModelIndicator() && 
       (Loci->GetLengthOfGenome() + Loci->GetLengthOfXchrm() > 0.0))
      L->UpdateGlobalSumIntensities(IC, (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1)); 
    // leaves individuals with HMM probs bad, stored likelihood ok
    // this function also sets locus correlations in Chromosomes
  }
//   if(options->getHapMixModelIndicator() && (iteration%2)){
//     if(isMaster || isWorker){
//       L->UpdateSumIntensitiesByRandomWalk(IC,  
// 			      (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1) );
//     }
//   }
  //find posterior modes of individual admixture at end of burn-in
  //set Chib numerator
  if(!anneal && iteration == options->getBurnIn() && (options->getChibIndicator() || strlen(options->getIndAdmixModeFilename()))) {
    IC->FindPosteriorModes(options, R, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), A, PopulationLabels);
    if( options->getChibIndicator() ) {
      IC->setChibNumerator(options, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), A);
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Update individual-level parameters, sampling locus ancestry states
  // then update jump indicators (+/- num arrivals if required for conjugate update of admixture or rho
  if(isMaster || isWorker){
    IC->SampleLocusAncestry(iteration, options, R, L->getpoptheta(), L->getalpha(), S->getAffectedsOnlyTest(), anneal);
  }

  if(isWorker || isFreqSampler) {
#ifdef PARALLEL
    MPE_Log_event(13, iteration, "sampleHapPairs");
#endif
    IC->SampleHapPairs(options, A, Loci, anneal); // loops over individuals to sample hap pairs then increment allele counts
#ifdef PARALLEL
    MPE_Log_event(14, iteration, "sampledHapPairs");
#endif
  }
  
#ifdef PARALLEL
  if(isWorker || isFreqSampler){
    A->SumAlleleCountsOverProcesses(options->getPopulations());
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
      S->Update(R);//score tests evaluated for first outcome var only
#ifdef PARALLEL
      MPE_Log_event(22, iteration, "ScoreTestUpdateEnd");
#endif
    }
    if(options->getTestForResidualAllelicAssoc()){
#ifdef PARALLEL
      MPE_Log_event(21, iteration, "ResLDTestStart");
#endif
      S->UpdateScoresForResidualAllelicAssociation(A->GetAlleleFreqs());
#ifdef PARALLEL
      MPE_Log_event(22, iteration, "ResLDTestEnd");
#endif
    }
  }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // sample individual admixture and sum-intensities 
  IC->SampleParameters(iteration, options, R, L->getpoptheta(), L->getalpha(),  
		       L->getrhoalpha(), L->getrhobeta(), anneal);
  // stored HMM likelihoods will now be bad if the sum-intensities are set at individual level
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if( options->getChibIndicator() && !anneal && iteration >= options->getBurnIn() )
    IC->updateChib(options, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), A);      
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // update allele frequencies conditional on locus ancestry states
  // TODO: this requires fixing to anneal allele freqs for historicallelefreq model
  if( !options->getFixedAlleleFreqs()){
    if(isFreqSampler){
#ifdef PARALLEL
      MPE_Log_event(7, iteration, "SampleFreqs");
#endif
      A->Update(IC, (iteration > options->getBurnIn() && !anneal), coolness);
#ifdef PARALLEL
    MPE_Log_event(8, iteration, "SampledFreqs");
#endif
    }
#ifdef PARALLEL
    if(isWorker || isFreqSampler ) { 
      A->BroadcastAlleleFreqs();
    }
#endif
  }//end if random allele freqs
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // even for fixed allele freqs, must reset annealed genotype probs as unnannealed

  if(isWorker && (!options->getFixedAlleleFreqs() || anneal)) {   
#ifdef PARALLEL
    MPE_Log_event(11, iteration, "setGenotypeProbs"); 
#endif
    IC->setGenotypeProbs(Loci, A); // sets unannealed probs ready for getEnergy
    IC->HMMIsBad(true); // update of allele freqs sets HMM probs and stored loglikelihoods as bad
#ifdef PARALLEL
    MPE_Log_event(12, iteration, "GenotypeProbsSet"); 
#endif
  } // update of allele freqs sets HMM probs and stored loglikelihoods as bad
  
  // next update of stored loglikelihoods will be from getEnergy if not annealing run, from updateRhowithRW if globalrho, 
  // or from update of individual-level parameters otherwise
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(isMaster || isWorker){
    if(options->getHapMixModelIndicator()){
      //L->UpdateGlobalTheta(iteration, IC);
      //conjugate sampler, using numbers of arrivals. NB: requires sampling of jump indicators in Individuals
      //L->SampleSumIntensities(IC->getSumNumArrivals(), IC->getSize(), 
      //Hamiltonian Sampler, using sampled ancestry states. NB: requires accumulating of SumAncestry in IC
      L->SampleSumIntensities(IC->getSumAncestry(),  
			      (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1) );
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else
    if(!options->getHapMixModelIndicator())
      //update population admixture Dirichlet parameters conditional on individual admixture
      L->UpdatePopAdmixParams(iteration, IC, Log);
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

void OutputParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, vector<Regression *>&R, 
		      const AdmixOptions *options, LogWriter& Log){
  // fix so that params can be output to console  
  Log.setDisplayMode(Quiet);
  if(options->getIndAdmixHierIndicator()  ){
    //output population-level parameters only when there is a hierarchical model on indadmixture
    // ** pop admixture, sumintensities
    if(options->getPopulations() > 1) L->OutputParams(iteration, Log);
    // ** dispersion parameter (if dispersion model)
    A->OutputEta(iteration, options, Log);
  }
  if(options->getHapMixModelIndicator() && (options->getDisplayLevel() > 2))cout << A->getHapMixPriorRate() << " " ;
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
    if(options->getOutputAlleleFreq()){
      if(options->getHapMixModelIndicator())
	A->OutputPriorParams();
      else
	A->OutputAlleleFreqs();
    }
  }
  // cout << endl;
}

void WriteIterationNumber(const int iteration, const int width, int displayLevel) {
  if( displayLevel > 2 ) { 
    cout << setiosflags( ios::fixed );
    cout.width(width );
    cout << "\n"<< iteration << " ";  
  }
  else if( displayLevel > 1 && iteration >0) { // display iteration counter only
    cout << "\rIterations so far: " << iteration;
  }
  cout.flush(); 
}


void OutputErgodicAvgDeviance(int samples, double & SumEnergy, double & SumEnergySq, std::ofstream *avgstream) {
  double EAvDeviance, EVarDeviance;
  EAvDeviance = 2.0*SumEnergy / (double) samples;//ergodic average of deviance
  EVarDeviance = 4.0 * SumEnergySq / (double)samples - EAvDeviance*EAvDeviance;//ergodic variance of deviance 
  *avgstream << EAvDeviance << " "<< EVarDeviance <<" ";
}
void PrintCopyrightNotice(LogWriter& Log){
  Log.setDisplayMode(On);
  cout << endl;
  Log << "-------------------------------------------------------\n"
      << "            ** ADMIXMAP (v" << ADMIXMAP_VERSION
#ifdef PARALLEL
      << " (Parallel) "
#endif
      << ") **\n"
      << "-------------------------------------------------------\n";
  Log.setDisplayMode(Quiet);
  cout << "Copyright(c) 2002-2006 " << endl
       << "David O'Donnell, Clive Hoggart and Paul McKeigue" << endl
       << "Send any comments or queries to david . odonnell@ucd.ie"<<endl
       << "-------------------------------------------------------"<<endl
       << "This program is free software distributed WITHOUT ANY WARRANTY " <<endl
       << "under the terms of the GNU General Public License. \nSee the file COPYING for details." <<endl
       << "-------------------------------------------------------" << endl;
}

void PrintOptionsMessage() {
  cout << "You must specify an options file\n"
       << "Consult the manual for a list of user options."
       << endl;
}

void MakeResultsDir(const char* dirname, bool verbose){
  if(strlen(dirname) && strcmp(dirname, ".") && strcmp(dirname, "..")){
    DIR *pdir;
    struct dirent *pent;
    
    string dirpath = "./";
    dirpath.append(dirname);
    pdir=opendir(dirpath.c_str()); //"." refers to the current dir
    if (!pdir){//dir does not exist
      cout << "Creating directory "<< dirpath <<endl;
      //this block is safer but only works in Windows
      //int status = mkdir(dirpath.c_str()/*, POSIX permissions go here*/);
      //if(status){
      //cout<< "Unable to create directory. Exiting." << endl;
      //exit(1);
      //}
      //KLUDGE: 'mkdir' not guaranteed to work on all systems; no error-checking
      //should be ok for normal use
      string cmd = "mkdir ";cmd.append(dirname);
      system(cmd.c_str());
    }
    else {
      cout << "Directory \"" << dirname << "\" exists. Contents will be deleted."<<endl;
      
      //list and delete contents of directory
      errno=0;
      while ((pent=readdir(pdir))){//read filenames
	if(strcmp(pent->d_name, ".") && strcmp(pent->d_name, "..")){//skip . and ..
	  string filepath = dirpath + "/"; 
	  filepath.append(pent->d_name);
	  if(verbose)
	    cout << "Deleting  " <<  filepath <<endl;
	  remove(filepath.c_str());//delete
	}
      }
      // cout << "errno = " << errno << endl << flush;
      if (errno && !errno==2){ // empty directory sets errno=2
	cerr << "readdir() failure; terminating";
	exit(1);
      }
      closedir(pdir);
      //rmdir(dirpath.c_str());
    }
  }
  else {
    cerr << "Invalid resultsdir. Exiting\n";
    exit(1);
  }
}

void ThrowException(const string& msg, LogWriter & Log){
#ifdef PARALLEL//print error message to screen as only master is allowed write with LogWriter
  Log << On << "rank " << MPI::COMM_WORLD.Get_rank() << ": " << msg << "\n Exiting...\n";
  Log.ProcessingTime();
  MPI::COMM_WORLD.Abort(1);
#else
  Log << On << "\n" << msg << "\n Exiting...\n";
  Log.ProcessingTime();
  exit(1);
#endif
}

