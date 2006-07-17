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
#ifdef PARALLEL
#include <mpi++.h>
#include <mpe.h>

MPI::Intracomm workers_and_master, workers_and_freqs;
#endif

#define MAXNUMOPTIONS 50//maximum number of options specifiable.
#define MAXOPTIONLENGTH 1024//maximum number of characters in an option line (excluding spaces)

using namespace std;
double coolness = 1.0; // default

int ReadArgsFromFile(char* filename, int* xargc, char **xargv);
void MakeResultsDir(const char* dirname, bool verbose);
void InitializeErgodicAvgFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
			      LogWriter &Log, std::ofstream *avgstream, const string* const PopulationLabels);
void UpdateParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, vector<Regression *>&R, const AdmixOptions *options, 
		      const Genome *Loci, ScoreTests *S, LogWriter& Log, const std::string* const PopulationLabels, 
		      double coolness, bool anneal);
void OutputParameters(int iteration, IndividualCollection *IC, Latent *L, AlleleFreqs *A, vector<Regression *>&R, const AdmixOptions *options, 
		      LogWriter& Log);
void WriteIterationNumber(const int iteration, const int width, int displayLevel);

void PrintCopyrightNotice(); 

void doIterations(const int & samples, const int & burnin, IndividualCollection *IC, Latent & L, AlleleFreqs  & A, 
		  vector<Regression *>&R, 
		  AdmixOptions & options,  const Genome  & Loci, LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
		  const double coolness, bool AnnealedRun, ofstream & loglikelihoodfile, 
		  ScoreTests & Scoretests, DispersionTest & DispTest, StratificationTest & StratTest, 
		  MisSpecAlleleFreqTest & AlleleFreqTest, HWTest & HWtest, ofstream & avgstream, InputData & data, 
		  const double* Coolnesses);

void OutputErgodicAvgDeviance(int samples, double & SumEnergy, double & SumEnergySq, std::ofstream *avgstream);

void PrintOptionsMessage();

int main( int argc , char** argv ){
#ifdef PARALLEL
  MPI::Init(argc, argv);
  const int rank = MPI::COMM_WORLD.Get_rank();

  //set up communicators
  int ranks[1];
  MPI::Group world_group = MPI::COMM_WORLD.Get_group();
  ranks[0] = 1;
  MPI::Group master_group = world_group.Excl(1, ranks);//exclude process (1) that will update freqs
  ranks[0] = 0;
  MPI::Group freqs_group = world_group.Excl(1, ranks);//exclude process 0 (master)
  workers_and_master = MPI::COMM_WORLD.Create(master_group);
  workers_and_freqs = MPI::COMM_WORLD.Create(freqs_group);
  world_group.Free();
  master_group.Free();
  freqs_group.Free();

  double t1 = 0, t2 = 0, t3 = 0, t0 = MPI::Wtime();
  //define states for logging
  MPE_Init_log();
  if(rank==0){
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
#else
  const int rank = -1;//cannot be 0 as process with rank 0 is excluded from some function calls
#endif

  int    xargc = argc;
  char **xargv = argv;    

   if (argc < 2) {
    PrintOptionsMessage();
    exit(1); 
  } else if (argc == 2) {     // using options text file        
    xargc = 1;//NB initialise to 1 to mimic argc (arg 0 is prog name), otherwise first option is ignored later
    xargv = new char*[MAXNUMOPTIONS];  
    ReadArgsFromFile(argv[1], &xargc, xargv);        
  }
  else {//fix broken arguments
    xargv = new char*[argc];
    xargv[1] = new char[strlen(argv[1])];
    strcpy(xargv[1], argv[1]);
    int ii = 2;
    for(int i = 2; i <argc; ++i)
      if( argv[i][0] != '-')  {//any 'args' not starting with '-' are appended to the previous line 
	strcat(xargv[ii-1], argv[i]);
	-- xargc;
      }
      else {
	  xargv[ii] = argv[i];
	++ii;
      }
  }
 
  // ******************* PRIMARY INITIALIZATION ********************************************************************************
  //read user options
  AdmixOptions options(xargc, xargv);
  if(rank<1){
    //if(options.getDisplayLevel()>0 )
    PrintCopyrightNotice();
	
    MakeResultsDir(options.getResultsDir().c_str(), (options.getDisplayLevel()>1));
  }
 
  //open logfile, start timer and print start message
  LogWriter Log(options.getLogFilename(), (bool)(options.getDisplayLevel()>1));
  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
  if(rank<1)Log.StartMessage();

  try{  
    Rand RNG;//allocate random number generator
    RNG.setSeed( options.getSeed() );  // set random number seed
  
    InputData data; //read data files and check (except allelefreq files)
    data.readData(&options, Log, rank);//also sets 'numberofoutcomes' and 'populations' options

     //check user options
    if(options.checkOptions(Log, data.getNumberOfIndividuals())){
      cout << endl << endl;
      exit(1);
    }
  
    //print user options to args.txt; must be done after all options are set
    if(rank<1)options.PrintOptions();
    Genome Loci;
    Loci.Initialise(&data, options.getPopulations(), Log, rank);//reads locusfile and creates CompositeLocus objects
    if(rank==1 || rank==-1){
      //print table of loci for R script to read
      string locustable = options.getResultsDir();
      locustable.append("/LocusTable.txt");
      Loci.PrintLocusTable(locustable.c_str(), data.getLocusMatrix().getCol(1));
      locustable.clear();
    }

    AlleleFreqs A(&Loci);
    if(rank ==-1 || rank ==1)//allele freq updater only
      A.Initialise(&options, &data, Log); //checks allelefreq files, initialises allele freqs and finishes setting up Composite Loci
    if(rank!=0)A.AllocateAlleleCountArrays(options.getPopulations());
#ifdef PARALLEL
    //broadcast initial values of freqs
    if(rank!=0)A.BroadcastAlleleFreqs(workers_and_freqs);

#endif

    IndividualCollection *IC = new IndividualCollection(&options, &data, &Loci);//NB call after A Initialise
    if(rank!=1)IC->LoadData(&options, &data);                             //and before L and R Initialise
    if(options.getNumberOfOutcomes()>0 && rank<1)IC->OpenExpectedYFile(options.getResidualFilename(), Log );
    if(rank!=0 && rank!=1)IC->setGenotypeProbs(&Loci, &A); // sets unannealed probs

    Latent L( &options, &Loci);    
    if(rank!=1)L.Initialise(IC->getSize(), data.GetPopLabels(), Log);

    vector<Regression*> R;//vector of regression pointers
    if (options.getNumberOfOutcomes()>0 && rank !=1){
      for(int r = 0; r < options.getNumberOfOutcomes(); ++r){
	//determine regression type and allocate regression objects
	if( IC->getOutcomeType(r)== Binary ) R.push_back( new LogisticRegression() );
	else if( IC->getOutcomeType(r)== Continuous ) R.push_back( new LinearRegression());

	if(rank < 1) R[r]->Initialise(r, options.getRegressionPriorPrecision(), IC, Log);
	else R[r]->Initialise(r, IC);
      }
      if(rank<1)Regression::OpenOutputFile(&options, IC, data.GetPopLabels(), Log);  
    }
    
    if( (options.isGlobalRho() || options.getHapMixModelIndicator()) && (rank>1 || rank==-1)) {
      Loci.SetLocusCorrelation(L.getrho());
      if(options.getHapMixModelIndicator())
	for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ )
	  //set global state arrival probs in hapmixmodel
	  //TODO: can skip this if xonly analysis with no females
	  Loci.getChromosome(j)->SetStateArrivalProbs(L.getGlobalTheta(), options.isRandomMatingModel(), true);
    }
    
    if(rank !=1)IC->Initialise(&options, &Loci, data.GetPopLabels(), L.getalpha(), /*L.getrhoalpha(), L.getrhobeta(),*/ Log);
  
    //set expected Outcome
    if(rank < 1)
      for(int r = 0; r < options.getNumberOfOutcomes(); ++r)
	R[r]->SetExpectedY(IC);

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
      if( options.getScoreTestIndicator() && rank!=1){
#ifdef PARALLEL
	Scoretests.SetComm(&workers_and_master, &(data.getLocusLabels()));
#endif
	Scoretests.Initialise(&options, IC, &Loci, data.GetPopLabels(), Log);
      }
      //if(rank==0){
      if( options.getTestForDispersion() ){
	DispTest.Initialise(&options, Log, Loci.GetNumberOfCompositeLoci());    
      }
      //if( options.getStratificationTest() )
      //StratTest.Initialize( &options, Loci, IC, Log);
      if( options.getTestForMisspecifiedAlleleFreqs() || options.getTestForMisspecifiedAlleleFreqs2())
	AlleleFreqTest.Initialise(&options, &Loci, Log );  
      if( options.getHWTestIndicator() )
	HWtest.Initialise(&options, Loci.GetTotalNumberOfLoci(), Log);
      if(rank<1)
	InitializeErgodicAvgFile(&options, IC, Log, &avgstream,data.GetPopLabels());
      //}

      string s = options.getResultsDir()+"/loglikelihoodfile.txt";
      ofstream loglikelihoodfile(s.c_str());
    
      // ******************* initialise stuff for annealing ************************************************
      double IntervalRatio = 1.03; // size of increments of coolness increases geometrically
      int NumAnnealedRuns = options.getNumAnnealedRuns(); // number of annealed runs excluding last run at coolness of 1
      int samples = 1; // default for annealing runs - overridden for final unannealed run or if AnnealIndicator = true 
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
      IntervalWidths = new double[NumAnnealedRuns + 1];
      Coolnesses = new double[NumAnnealedRuns + 1];
      Coolnesses[0] = 0.0;
      if(NumAnnealedRuns > 0) {
	// initial increment of coolness from 0 is set so that geometric series of increments will sum to 1 
	// after NumAnnealedRuns additional terms
	IntervalWidths[1] = (IntervalRatio - 1.0) /(pow(IntervalRatio, NumAnnealedRuns) - 1.0 ); 
	Coolnesses[1] = IntervalWidths[1];
	if(NumAnnealedRuns > 1) {
	  for(int run=2; run < NumAnnealedRuns+1; ++run) {
	    IntervalWidths[run] = IntervalWidths[run - 1] * IntervalRatio; // geometric increments in interval width
	    Coolnesses[run] = Coolnesses[run - 1] + IntervalWidths[run];  
	  }
	}
      }
      Coolnesses[NumAnnealedRuns] = 1.0;
    
      if(options.getThermoIndicator()) { // set up output for thermodynamic integration
	if(rank<1){		
	  string s = options.getResultsDir()+"/annealmon.txt";
	  annealstream.open(s.c_str());
	  annealstream << "Coolness\tMeanEnergy\tVarEnergy\tlogEvidence" << endl;
	}
	samples = options.getTotalSamples();
	burnin = options.getBurnIn();
      }
      Log.setDisplayMode(On);
      if( options.getTestOneIndivIndicator() )NumAnnealedRuns = 0;
      if(rank<1){
	if(NumAnnealedRuns > 0) {
	  Log << NumAnnealedRuns << " annealing runs of " << samples 
	      << " iteration(s) followed by final run of "; 
	}
	Log << options.getTotalSamples() << " iterations at ";
	if( options.getTestOneIndivIndicator() )Log << options.getNumAnnealedRuns()+1 
						    <<" coolnesses for test individual. Other individuals at ";
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
	for(int run=0; run < NumAnnealedRuns + 1; ++run) { //loop over coolnesses from 0 to 1
	  // should call a posterior mode-finding algorithm before last run at coolness of 1
	  //resets for start of each run
	  SumEnergy = 0.0;//cumulative sum of modified loglikelihood
	  SumEnergySq = 0.0;//cumulative sum of square of modified loglikelihood

	  if(run == NumAnnealedRuns) {
	    AnnealedRun = false;
	    samples = options.getTotalSamples();//redundant?
	    burnin = options.getBurnIn();
	  } else AnnealedRun = true; 
	  coolness = Coolnesses[run];
	  if(NumAnnealedRuns > 0) {
	    cout <<"\rSampling at coolness of " << coolness << "         " << flush;
	    IC->resetStepSizeApproximators(NumAnnealedRuns); // reset k <= NumAnnealedRuns in step size tuners
	  }
	  // accumulate scalars SumEnergy and SumEnergySq at this coolness
	  // array Coolnesses is not used unless TestOneIndivIndicator is true
	  doIterations(samples, burnin, IC, L, A, R, options, Loci, Log, SumEnergy, SumEnergySq, coolness, AnnealedRun, 
		       loglikelihoodfile, Scoretests, DispTest, StratTest, AlleleFreqTest, HWtest, avgstream, data, Coolnesses);
#ifdef PARALLEL
	  t2 = MPI::Wtime()-t1;
#endif
	  if(rank<1){	  
	    //calculate mean and variance of energy at this coolness
	    MeanEnergy = SumEnergy / ((double)options.getTotalSamples() - options.getBurnIn());
	    VarEnergy  = SumEnergySq / ((double)options.getTotalSamples() - options.getBurnIn()) - MeanEnergy * MeanEnergy;
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
      if(rank<1)
	cout<< "\nIterations completed                       \n" << flush;
    
      // *************************** OUTPUT AT END ***********************************************************

      if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);	
      else Log.setDisplayMode(On);
      if( options.getChibIndicator()) {
	//IC->OutputChibEstimates(options.isRandomMatingModel(), Log, options.getPopulations());
	//MLEs of admixture & sumintensities used in Chib algorithm to estimate marginal likelihood
	if(IC->getSize()==1) IC->OutputChibResults(Log);
      }
		
      double Information = -LogEvidence - MeanEnergy;
      double MeanDeviance = 2.0 * MeanEnergy; 
      double VarDeviance = 4.0 * VarEnergy;
      Log << "\nMeanDeviance(D_bar)\t" << MeanDeviance << "\n"
	  << "VarDeviance(V)\t" << VarDeviance << "\n"
	  << "PritchardStat(D_bar+0.25V)\t" << MeanDeviance + 0.25*VarDeviance << "\n";
      double D_hat = IC->getDevianceAtPosteriorMean(&options, R, &Loci, Log, L.getSumLogRho(), Loci.GetNumberOfChromosomes(), &A);
      double pD = MeanDeviance - D_hat;
      double DIC = MeanDeviance + pD;
      Log << "DevianceAtPosteriorMean(D_hat)\t" << D_hat << "\n"
	  << "EffectiveNumParameters(pD)\t" << pD << "\n"
	  << "DevianceInformationCriterion\t" << DIC << "\n\n"; 
		
      if(options.getThermoIndicator()){
	Log << "thermodynamic integration for marginal likelihood yields:\n";
	Log << "LogEvidence " <<  LogEvidence << "\n"; 
	Log << "Information (negative entropy, measured in nats) " << Information << "\n";
      }
		
      //Residuals
      if(rank <1 && options.getNumberOfOutcomes() > 0)
	IC->//OutputResiduals(options.getResidualFilename(), data.getOutcomeLabels(), options.getTotalSamples()-options.getBurnIn());
	  FinishWritingEYAsRObject((options.getTotalSamples()-options.getBurnIn())/ options.getSampleEvery(), data.getOutcomeLabels());
      //FST
      if( strlen( options.getHistoricalAlleleFreqFilename()) && (rank==-1 || rank==1)  ){
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

      if( options.getScoreTestIndicator() && (rank <1) ) {
	//finish writing score test output as R objects
	Scoretests.ROutput();
	//write final tables
	Scoretests.Output(options.getTotalSamples() - options.getBurnIn(), data.GetPopLabels(), true);
      }
    
      //output to likelihood ratio file
      if(options.getTestForAffectedsOnly())
	Individual::OutputLikRatios(options.getLikRatioFilename(), options.getTotalSamples()-options.getBurnIn(), data.GetPopLabels());
		
      if(annealstream.is_open())annealstream.close();
      if(avgstream.is_open())avgstream.close();

    }//end else

    delete IC;//must call explicitly so IndAdmixOutputter destructor finishes writing to indadmixture.txt
    for(int r = 0; r < options.getNumberOfOutcomes(); ++r)delete R[r];
    R.clear();

    if((rank==-1 || rank==1) && !options.getHapMixModelIndicator())
      A.CloseOutputFile((options.getTotalSamples() - options.getBurnIn())/options.getSampleEvery(), data.GetPopLabels());
    if(rank<1)cout << "Output to files completed\n" << flush;
    
    // ******************* acceptance rates - output to screen and log ***************************
    if( options.getIndAdmixHierIndicator() ){
      if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
      else Log.setDisplayMode(On);
      if(options.getPopulations() > 1 && (rank<1)){
	L.printAcceptanceRates(Log);
      }
      if(options.getCorrelatedAlleleFreqs()){
	Log<< "Expected acceptance rates in sampler for allele frequency proportion parameters: \n";
	for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci(); i++){
	  if(Loci(i)->GetNumberOfStates()>2)
	    Log << A.getAlphaSamplerAcceptanceRate(i) << " ";
	}
	Log<< "Expected acceptance rate in sampler for allele frequency dispersion parameter: \n";
	Log << A.getEtaRWSamplerAcceptanceRate(0)
	    << "\nwith final step sizes of \n";
	for(unsigned int i = 0; i < Loci.GetNumberOfCompositeLoci(); i++){
	  if(Loci(i)->GetNumberOfStates()>2)
	    Log << A.getAlphaSamplerStepsize(i) << " " ;
	}
	Log <<  A.getEtaRWSamplerStepsize(0) << "\n" ;
      }
    
#if ETASAMPLER ==1
      if( strlen( options.getHistoricalAlleleFreqFilename() )){
	Log << "Expected acceptance rates in allele frequency dispersion parameter samplers:\n ";
	for(int k = 0; k < options.getPopulations(); ++k){Log << A.getEtaRWSamplerAcceptanceRate(k)<< " " ;}
	Log << "\nwith final step sizes of ";
	for(int k = 0; k < options.getPopulations(); ++k){Log <<  A.getEtaRWSamplerStepsize(k) << " ";}
	Log << "\n";
      }
#endif
    }
    if(options.getHapMixModelIndicator() && (rank == -1 || rank== 1)){
      Log << "Average expected Acceptance rate in allele frequency prior parameter sampler:\n" << A.getHapMixPriorSamplerAcceptanceRate()
	  << "\nwith average final step size of " << A.getHapMixPriorSamplerStepSize() << "\n";
    }
    if(rank<1){
      if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
      Log.ProcessingTime();
    }

#ifdef PARALLEL
    t3 = MPI::Wtime()-t2;
#endif
    //if(options.getDisplayLevel()>0){
#ifdef PARALLEL
    cout << "Rank " << rank << " finished.\n Initialization: " << t1 <<"s Main loop: "<< t2 << "s Finalization:" << t3 << "s.\n";
#else
    cout << "Finished" << endl;
#endif
    //}
    if(argc ==2){//options file supplied so need to delete xargv
      for(int i = 1; i < xargc; ++i)delete[] xargv[i];
      delete[] xargv;
    }
    else
      if(xargv != argv){
	delete[] argv[1];
	delete[] xargv;
      }
    //end try block	
  } catch (string msg) {//catch any stray error messages thrown upwards
    Log.setDisplayMode(On);
    Log << "\n" << msg << "\n Exiting...\n";
    Log.ProcessingTime();
#ifdef PARALLEL//print error message to screen as only master is allowed write with LogWriter
    if(rank>0)cerr << "rank " << rank << ": " << msg << endl;
    MPI::COMM_WORLD.Abort(1);
#else
    exit(1);
#endif
  }
  catch (char *msg){//in case error messages thrown as char arrays instead of strings
    throw string(msg);
  }

#ifdef PARALLEL
  catch(MPI::Exception e){
    cout << "Error in process " << rank << ": " << e.Get_error_code() << endl;
    MPI::COMM_WORLD.Abort(1);
  }
#endif
  catch(...){
    cout << "Unknown exception occurred. Contact the program authors for assistance" << endl;
    exit(1);
  }
#ifdef PARALLEL
  //MPI::COMM_WORLD.Barrier();
  MPE_Finish_log("admixmap");

  if(rank!=1)workers_and_master.Free();
  if(rank!=0)workers_and_freqs.Free();
  MPI_Finalize();
#endif
  if(rank<1){//print line of *s
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
#ifdef PARALLEL
  const int rank = MPI::COMM_WORLD.Get_rank();
#else
  const int rank = 0;
#endif

  double Energy = 0.0;
  if(rank==0 && !AnnealedRun) cout << endl;
  for( int iteration = 0; iteration <= samples; iteration++ ) {
    if(rank!=1 && iteration > burnin) {
      //accumulate energy as minus loglikelihood, calculated using unnanealed genotype probs
      if( !options.getTestOneIndivIndicator() ) {
	Energy = IC->getEnergy(&options, R, AnnealedRun); // should store loglikelihood if not AnnealedRun
	if(rank==0){
	  SumEnergy += Energy;
	  SumEnergySq += Energy*Energy;
	  // write to file if not AnnealedRun
	  if(!AnnealedRun)loglikelihoodfile << iteration<< "\t" << Energy <<endl;
	}
      } else {  
	IC->accumulateEnergyArrays(&options);
      }
    }
    if( rank==0 && !AnnealedRun &&  !(iteration % options.getSampleEvery()) ) {
      WriteIterationNumber(iteration, (int)log10((double) samples+1 ), options.getDisplayLevel());
    }
    
    // if annealed run, anneal genotype probs - for testindiv only if testsingleindiv indicator set in IC
    if((rank!=1) && (AnnealedRun || options.getTestOneIndivIndicator())) IC->annealGenotypeProbs(Loci.GetNumberOfChromosomes(), coolness, Coolnesses); 
    
    UpdateParameters(iteration, IC, &L, &A, R, &options, &Loci, &Scoretests, Log, data.GetPopLabels(), coolness, AnnealedRun);

    Log.setDisplayMode(Quiet);
 
    if(!AnnealedRun){    
      // output every 'getSampleEvery()' iterations
      if(!(iteration % options.getSampleEvery()) && rank==0)
	OutputParameters(iteration, IC, &L, &A, R, &options, Log);
      
      // ** set merged haplotypes for allelic association score test 
      if( rank!=1 && iteration == options.getBurnIn() ){
#ifndef PARALLEL
	if(options.getTestForAllelicAssociation())
	  Scoretests.SetAllelicAssociationTest(L.getalpha0());
#endif
	if( options.getStratificationTest() )
	  StratTest.Initialize( &options, Loci, IC, Log, rank);
      }
	    
      //Updates and Output after BurnIn     
      if( !AnnealedRun && iteration > burnin && rank<1){
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
	      A.OutputErgodicAvg(samples, &avgstream);//dispersion parameter in dispersion model
	    }
	    for(int r = 0; r < options.getNumberOfOutcomes(); ++r)//regression params
	      R[r]->OutputErgodicAvg(samples, &avgstream);

	    OutputErgodicAvgDeviance(samples, SumEnergy, SumEnergySq, &avgstream);
	    if(options.getChibIndicator()) IC->OutputErgodicChib(&avgstream, options.getFixedAlleleFreqs());
	    avgstream << endl;
	  }
	  //Score Test output
	  if( options.getScoreTestIndicator() )  Scoretests.Output(iteration-burnin, data.GetPopLabels(), false);
	}//end "if every'*10" block
      }//end "if after BurnIn" block
    } // end "if not AnnealedRun" block
	
  }// end loop over iterations
}

int ReadArgsFromFile(char* filename, int* xargc, char **xargv){
  ifstream fin(filename);
  if (0 == filename || 0 == strlen(filename)) return 1;
  if (!fin.is_open()) {
    string msg = "Cannot open file \"";
    msg += filename;
    msg += "\". Aborting";
    cerr << msg << endl;
    exit(1);
  } 

  std::string str;
  //read in line from file
  while (getline(fin,str,'\n')){// ## apparent memory leak 

    if( str.find_first_of("#") < str.length() ) str.erase( str.find_first_of("#") );//ignore #comments
    if(str.find_first_not_of(" \t\n\r") < str.length() )//skip blank lines. 
      {   
	str.erase(0, str.find_first_not_of(" \t\n\r") );//trim leading whitespace
	//trim remaining whitespace
	str.erase( str.find_last_not_of(" \t\n\r") + 1 );//trailing whitespace
	if( str.find_first_of(" \t\n\r") <= str.length() ){//check for any whitespace left
	  string::size_type eq = str.find("="), pos = str.find_first_of(" \t\n\r");
	  if( pos < eq )//check for space before '='
	    str.erase( pos, eq - pos );//remove space before '='
	  //str.erase( str.find_first_of(" \t\n\r"),str.find_last_of(" \t\n") - str.find_first_of(" \t\n\r") +1 );//after '='
	  eq = str.find("=");
	  pos = str.find_first_of(" \t\n\r", eq);//position of first space after the = 
	  str.erase( pos, str.find_first_not_of(" \t\n", pos) - pos );//remove space after '='
	}
	//add line to xargv
	xargv[*xargc]=new char[MAXOPTIONLENGTH];
	strcpy(xargv[*xargc],"--");
	strcat(xargv[*xargc],str.c_str());
	++(*xargc);
      }
    str.clear();
  }
  fin.close();
  return 0;
}

//this function is here because three different objects have to write to avgstream
void InitializeErgodicAvgFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
			      LogWriter &Log, std::ofstream *avgstream, const std::string* const PopulationLabels){
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
    if(options->getHapMixModelIndicator()){
      *avgstream << "FreqPriorRate\t"; 
    }

    // Regression parameters
    if( options->getNumberOfOutcomes() > 0 ){
      for(int r = 0; r < individuals->getNumberOfOutcomeVars(); ++r){
	*avgstream << "intercept\t";
	if(strlen(options->getCovariatesFilename()) > 0){//if covariatesfile specified
	  for( int i = 0; i < individuals->GetNumberOfInputCovariates(); i++ ){
	    *avgstream << individuals->getCovariateLabels(i) << "\t";//write covariate labels to header
	  }
	}
	if( !options->getHapMixModelIndicator() && !options->getTestForAdmixtureAssociation() ){
	  for( int k = 1; k < options->getPopulations(); k++ ){
	    *avgstream << PopulationLabels[k] << "\t";//write population labels (admixture covariates) to header
	  }
	}
	if( individuals->getOutcomeType(r)==0 )//linear regression
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
		      const std::string* const PopulationLabels, double coolness, bool anneal){
#ifdef PARALLEL
  const int rank = MPI::COMM_WORLD.Get_rank();
#else 
  const int rank = -1;
#endif
  A->ResetAlleleCounts(options->getPopulations());

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(rank<1){
    // ** update global sumintensities conditional on genotype probs and individual admixture proportions
    if((options->getPopulations() > 1) && options->getIndAdmixHierIndicator() && !options->getHapMixModelIndicator() && 
       (Loci->GetLengthOfGenome() + Loci->GetLengthOfXchrm() > 0.0))
      L->UpdateGlobalSumIntensities(IC, (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1)); 
    // leaves individuals with HMM probs bad, stored likelihood ok
    // this function also sets locus correlations in Chromosomes
  }

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
  if(rank!=1)IC->SampleLocusAncestry(iteration, options, R, L->getpoptheta(), L->getalpha(), anneal);

  if(rank!=0) {
#ifdef PARALLEL
    MPE_Log_event(13, iteration, "sampleHapPairs");
#endif
    IC->SampleHapPairs(options, A, Loci, anneal);
#ifdef PARALLEL
    MPE_Log_event(14, iteration, "sampledHapPairs");
#endif
  }
  
#ifdef PARALLEL
  if(rank>0){
    A->SumAlleleCountsOverProcesses(workers_and_freqs, options->getPopulations());
  }
#endif
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if( rank!=1 && !anneal && iteration > options->getBurnIn() ){
    //score tests
    if( options->getScoreTestIndicator() ){
      double dispersion = 1.0;
      if(rank != 1 && (options->getNumberOfOutcomes()>0) ) dispersion = R[0]->getDispersion();
#ifdef PARALLEL
      MPE_Log_event(21, iteration, "ScoreTestUpdatestart");
#endif
      S->Update(dispersion);//score tests evaluated for first outcome var only
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
  // this function samples individual admixture with conjugate update on odd-numbered iterations
  // samples admixture of test individuals at every iteration  
  IC->SampleParameters(iteration, options, R, L->getpoptheta(), L->getalpha(),  
		       L->getrhoalpha(), L->getrhobeta(), anneal);
  // stored HMM likelihoods will now be bad if the sum-intensities are set at individual level
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if( options->getChibIndicator() && !anneal && iteration >= options->getBurnIn() )
    IC->updateChib(options, L->getalpha(), L->getrhoalpha(), L->getrhobeta(), A);      

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // update allele frequencies conditional on locus ancestry states
  // TODO: this requires fixing to anneal allele freqs for historicallelefreq model
  if( (rank==-1 || rank==1) && A->IsRandom()){//may be a tautology, non freq updaters have RandomAlleleFreqs=false
#ifdef PARALLEL
    MPE_Log_event(7, iteration, "SampleFreqs");
#endif
    //bool thermoSampler = (anneal && options->getThermoIndicator() && !options->getTestOneIndivIndicator());
    A->Update(IC, (iteration > options->getBurnIn() && !anneal), coolness);

#ifdef PARALLEL
    MPE_Log_event(8, iteration, "SampledFreqs");
#endif
  }
  
  // even for fixed allele freqs, must reset annealed genotype probs as unnannealed
  if(rank!=0 ) { 
#ifdef PARALLEL
    if(!options->getFixedAlleleFreqs()) { 
      A->BroadcastAlleleFreqs(workers_and_freqs);
    }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
    if(rank!=1 && (!options->getFixedAlleleFreqs() || anneal)) {   
#ifdef PARALLEL
    MPE_Log_event(11, iteration, "setGenotypeProbs"); 
#endif
      IC->setGenotypeProbs(Loci, A); // sets unannealed probs ready for getEnergy
      IC->HMMIsBad(true); // update of allele freqs sets HMM probs and stored loglikelihoods as bad
#ifdef PARALLEL
      MPE_Log_event(12, iteration, "GenotypeProbsSet"); 
#endif
    } // update of allele freqs sets HMM probs and stored loglikelihoods as bad
  }   
  // next update of stored loglikelihoods will be from getEnergy if not annealing run, from updateRhowithRW if globalrho, 
  // or from update of individual-level parameters otherwise
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(options->getHapMixModelIndicator()){
    if(rank!=1){
#ifdef PARALLEL
      MPE_Log_event(1, iteration, "BarrierStart");
      workers_and_master.Barrier();//force master to wait until workers have finished
      MPE_Log_event(2, iteration, "BarrierEnd");
      
#endif
      //L->UpdateGlobalTheta(iteration, IC);
      //conjugate sampler, using numbers of arrivals. NB: requires sampling of jump indicators in Individuals
      //L->SampleSumIntensities(IC->getSumNumArrivals(), IC->getSize(), 
      //Hamiltonian Sampler, using sampled ancestry states. NB: requires accumulating of SumAncestry in IC
      L->SampleSumIntensities(IC->getSumAncestry(),  
			      (!anneal && iteration > options->getBurnIn() && options->getPopulations() > 1)
#ifdef PARALLEL
			      , workers_and_master
#endif
			      );
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  else
    //update population admixture Dirichlet parameters conditional on individual admixture
    L->UpdatePopAdmixParams(iteration, IC, Log);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // ** update regression parameters (if regression model) conditional on individual admixture
  if(rank != 1){
    bool condition = (!anneal && iteration > options->getBurnIn() && (rank <1));
    for(int r = 0; r < options->getNumberOfOutcomes(); ++r){
      R[r]->Update(condition, IC, coolness 
#ifdef PARALLEL
		  , workers_and_master
#endif
);
      //IC->UpdateSumResiduals();
      //output expected values of outcome variables to file every 'every' iterations after burnin
      if(condition && !(iteration % options->getSampleEvery()) ) IC->OutputExpectedY(r);
    }
  }
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
  for(int r = 0; r < options->getNumberOfOutcomes(); ++r)
    R[r]->Output(iteration, options, Log);
  
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
void PrintCopyrightNotice(){
  cout << "\n-------------------------------------------------------" << endl;
  cout << "            ** ADMIXMAP (v" << ADMIXMAP_VERSION << ") **" << endl;
  cout << "-------------------------------------------------------" << endl;
  cout << "Copyright(c) 2002-2006 " <<endl;
  cout << "David O'Donnell, Clive Hoggart and Paul McKeigue"<<endl;
  cout << "Send any comments or queries to david . odonnell@ucd.ie"<<endl;
  cout << "-------------------------------------------------------"<<endl;
  cout << "This program is free software distributed WITHOUT ANY WARRANTY " <<endl;
  cout << "under the terms of the GNU General Public License. \nSee the file COPYING for details." <<endl;
  cout << "-------------------------------------------------------" << endl;
}

void PrintOptionsMessage() {
  cout << "You must specify an options file or list of arguments on command line\n"
       << "Usage:\n"
       << "1. (not recommended) admixmap --[optionname]=[value] ...\n"
       << "2. admixmap [optionfile], where optionfile is a text file containg a list of user options\n"
       << "3. use a Perl script to call the program with command-line arguments. \nSee sample script supplied with this program.\n"
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
