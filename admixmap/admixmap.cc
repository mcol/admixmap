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

#define ADMIXMAP_VERSION 3.5
#define HAPMIXMAP_VERSION 0.2
#define MAXNUMOPTIONS 50//maximum number of options specifiable.
#define MAXOPTIONLENGTH 1024//maximum number of characters in an option line (excluding spaces)

using namespace std;
double coolness = 1.0; // default



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
  //  const bool isWorker = Comms::isWorker();

  if(isMaster){
    MakeResultsDir(options.getResultsDir().c_str(), false/*(options.getDisplayLevel()>2)*/);
  }
 
  //open logfile, start timer and print start message
  LogWriter Log(options.getLogFilename(), (bool)(options.getDisplayLevel()>1 && isMaster));
  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
  if(isMaster){
    //if(options.getDisplayLevel()>0 )
    if(options.getHapMixModelIndicator())
      PrintHAPMIXMAPCopyrightNotice(Log);
    else
      PrintADMIXMAPCopyrightNotice(Log);
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

    Model* M;
    if(options.getHapMixModelIndicator())M = new HapMixModel;
    else M = new AdmixMapModel;

    Genome Loci;
    Loci.Initialise(&data, options.getPopulations(), Log);//reads locusfile and creates CompositeLocus objects
    if(isFreqSampler){
      //print table of loci for R script to read
      string locustable = options.getResultsDir();
      locustable.append("/LocusTable.txt");
      Loci.PrintLocusTable(locustable.c_str(), data.getLocusMatrix().getCol(1));
      locustable.clear();
    }

    M->Initialise(Loci, options, data, Log);

    data.Delete();
 
    //  ******** single individual, one population, fixed allele frequencies  ***************************
    if( M->getNumIndividuals() == 1 && options.getPopulations() == 1 && strlen(options.getAlleleFreqFilename()) )
      // nothing to do except calculate likelihood
      M->getOnePopOneIndLogLikelihood(Log, data.GetPopLabels());
    else {
      // ******************* Initialize test objects and ergodicaveragefile *******************************
      M->InitialiseTests(options, data, Loci, Log);

      //open file to output loglikelihood
      string s = options.getResultsDir()+"/loglikelihoodfile.txt";
      ofstream loglikelihoodfile(s.c_str());
    
      // ******************* Set annealing schedule ************************************************
      double IntervalRatio = 1.03; // size of increments of coolness increases geometrically
      int NumAnnealedRuns = options.getNumAnnealedRuns(); // number of annealed runs excluding last run at coolness of 1
      double SumEnergy = 0.0, SumEnergySq = 0.0, LogEvidence = 0.0;
      double MeanEnergy = 0.0, VarEnergy = 0.0;
      double LastMeanEnergy = 0.0;
      double AISsumlogz = 0.0; //for computing marginal likelihood by Annealed Importance Sampling
    
      coolness = 1.0; // default
      bool AnnealedRun = false;
      std::ofstream annealstream;//for monitoring energy when annealing

      // set number of samples : 1 for annealing runs, "samples" option otherwise. Overriden for final, unannealed run with "thermo" option
      int samples, burnin;
      if(options.getThermoIndicator()) { // set up output for thermodynamic integration
	if(isMaster){		
	  string s = options.getResultsDir()+"/annealmon.txt";
	  annealstream.open(s.c_str());
	  annealstream << "Coolness\tMeanEnergy\tVarEnergy\tlogEvidence" << endl;
	}
	samples = options.getTotalSamples();
	burnin = options.getBurnIn();
      }
      else if(NumAnnealedRuns > 0){
	samples = 1;
	burnin = 1;
      }
      else{//default case: no thermo, no annealing
	samples = options.getTotalSamples();
	burnin = options.getBurnIn();
      }

      if( options.getTestOneIndivIndicator() )NumAnnealedRuns = 0;
      if(isMaster){
	if(NumAnnealedRuns > 0) {
	  Log << On << NumAnnealedRuns << " annealing runs of " << samples 
	      << " iteration(s) followed by final run of "; 
	}
	if(options.getThermoIndicator()) {
	  Log << 2*samples;  //last run is twice as long with thermo option
	} else { 
	  Log << options.getTotalSamples();
	} 
	Log << " iterations at ";
	if( options.getTestOneIndivIndicator() ) {
	  Log << options.getNumAnnealedRuns()+1 
	      <<" coolnesses for test individual. Other individuals at ";
	}
	Log << "coolness of 1\n";
      }
    
      // set annealing schedule
      double *IntervalWidths = 0;
      double *Coolnesses = 0; // 
      IntervalWidths = new double[NumAnnealedRuns+1];
      Coolnesses = new double[NumAnnealedRuns + 1];
      Coolnesses[0] = 1.0;//default for unnannealed run

      if(NumAnnealedRuns > 0) {
	Coolnesses[0] = 0.0; // change this if you want annealing to start somewhere other than 0;
	// set initial increment of coolness so that geometric series of NumAnnealedRuns increments 
	// will sum to 1 - Coolnesses[0] after NumAnnealedRuns + 1 terms
	IntervalWidths[0] = (1.0 - Coolnesses[0]) * (1.0 - IntervalRatio) /(1.0 - pow(IntervalRatio, NumAnnealedRuns)); 
	//Coolnesses[1] = Coolnesses[0] + IntervalWidths[0];
	if(NumAnnealedRuns > 1) {
	  for(int run=1; run <= NumAnnealedRuns; ++run) {
	    IntervalWidths[run] = IntervalWidths[run - 1] * IntervalRatio; // geometric increments in interval width
	    Coolnesses[run] = Coolnesses[run - 1] + IntervalWidths[run - 1];  
	  }
	}
	Coolnesses[NumAnnealedRuns] = 1.0;
      }
    
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
	    else{
	      samples = options.getTotalSamples();
	    }
	    burnin = options.getBurnIn();
	  } 
	  else {
	    AnnealedRun = true; 
	    coolness = Coolnesses[run];
	  }
	  if(NumAnnealedRuns > 0) {
	    cout <<"\rSampling at coolness of " << coolness << "       "<< flush;
	    // reset approximation series in step size tuners
	    int resetk = NumAnnealedRuns; //   
	    if(samples < NumAnnealedRuns) {// samples=1 if annealing without thermo integration
	      resetk = 1 + run;
	    }
	    M->ResetStepSizeApproximators(resetk); 
	    
	  }
	  // accumulate scalars SumEnergy and SumEnergySq at this coolness
	  // array Coolnesses is not used unless TestOneIndivIndicator is true
	  M->Iterate(samples, burnin, Coolnesses, run, options, data, Loci, Log, SumEnergy, SumEnergySq, AISsumlogz,
			 AnnealedRun, loglikelihoodfile);

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
	    cout << "\t\tMeanEnergy = " << MeanEnergy << "        " << flush;
	  } 
	} // end loop over coolnesses
      } 
      else { // evaluate energy for test individual only at all coolnesses simultaneously
	// call with argument AnnealedRun false - copies of test individual will be annealed anyway  
	M->Iterate(samples, burnin, Coolnesses, NumAnnealedRuns, options, data, Loci, Log, SumEnergy, SumEnergySq, AISsumlogz, false, 
		       loglikelihoodfile);
	// arrays of accumulated sums for energy and energy-squared have to be retrieved by function calls
	double *MeanEner = M->getSumEnergy(); 
	double *VarEner = M->getSumEnergySq();
	
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

      M->Finalize(options, Log, data, Loci);
      
      double Information = -LogEvidence - MeanEnergy;
      double MeanDeviance = 2.0 * MeanEnergy; 
      double VarDeviance = 4.0 * VarEnergy;
      Log << Quiet << "\nMeanDeviance(D_bar)\t" << MeanDeviance << "\n"
	  << "VarDeviance(V)\t" << VarDeviance << "\n"
	  << "PritchardStat(D_bar+0.25V)\t" << MeanDeviance + 0.25*VarDeviance << "\n";
      double D_hat = M->getDevianceAtPosteriorMean(&options, &Loci, Log);
      double pD = MeanDeviance - D_hat;
      double DIC = MeanDeviance + pD;
      Log << Quiet << "DevianceAtPosteriorMean(D_hat)\t" << D_hat << "\n"
	  << "EffectiveNumParameters(pD)\t" << pD << "\n"
	  << "DevianceInformationCriterion\t" << DIC << "\n\n"; 
		
      if(options.getThermoIndicator()){
	Log << "thermodynamic integration for marginal likelihood yields:\n";
	Log << "LogEvidence " <<  LogEvidence << "\n"; 
	Log << "Information (negative entropy, measured in nats) " << Information << "\n";

	Log << "\nAnnealed Importance Sampling estimates log marginal likelihood as " << AISsumlogz << "\n";
      }
		
      //Expected Outcome
      if(isMaster && options.getNumberOfOutcomes() > 0){
	Regression::FinishWritingEYAsRObject((options.getTotalSamples()-options.getBurnIn())/ options.getSampleEvery(), 
					     data.getOutcomeLabels());
      }

      if(annealstream.is_open())annealstream.close();

    }//end else

     if(isMaster)cout << "Output to files completed\n" << flush;
    
    // ******************* acceptance rates - output to screen and log ***************************
    if( options.getIndAdmixHierIndicator() ){
      M->PrintAcceptanceRates(options, Loci, Log);
    }
    delete M;
    if(isMaster){
      if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
      Log.ProcessingTime();
    }

#ifdef PARALLEL
    t3 = MPI::Wtime()-t2;
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
  cout << "Rank " << MPI::COMM_WORLD.Get_rank() << " finished.\n Initialization: " << t1 <<"s Main loop: "<< t2 << "s Finalization:" << t3 << "s.\n";
  MPI_Finalize();
#else
  cout << "Finished" << endl;
#endif

  if(isMaster){//print line of *s
    cout <<setfill('*') << setw(80) << "*" <<endl;
  }
  putenv("ADMIXMAPCLEANEXIT=1");
  return 0;
} //end of main

void PrintADMIXMAPCopyrightNotice(LogWriter& Log){
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
       << "-------------------------------------------------------"<<endl
       << "This program is free software distributed WITHOUT ANY WARRANTY " <<endl
       << "under the terms of the GNU General Public License. \nSee the file COPYING for details." <<endl
       << "-------------------------------------------------------" << endl;
}

void PrintHAPMIXMAPCopyrightNotice(LogWriter& Log){
  Log.setDisplayMode(On);
  cout << endl;
  Log << "-------------------------------------------------------\n"
      << "            ** HAPMIXMAP (v" << HAPMIXMAP_VERSION
#ifdef PARALLEL
      << " (Parallel) "
#endif
      << ") **\n"
      << "-------------------------------------------------------\n";
  Log.setDisplayMode(Quiet);
  cout << "Copyright(c) 2006 " << endl
       << "David O'Donnell and Paul McKeigue" << endl
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



