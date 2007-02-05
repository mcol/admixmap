/** 
 *   ADMIXMAP
 *   admixmap.cc 
 *   Top-level source file
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "AdmixMapModel.h"
#include <fstream>

#define ADMIXMAP_VERSION 3
#define SUBVERSION 5.3

using namespace std;

int main( int argc , char** argv ){
  //-v flag prints version number and exits
  if(argc==2 & !strcmp(argv[1], "-v")){
    LogWriter LW;
    PrintCopyrightNotice(LW);
    exit(0);
  }
  //-h flag or --help print a usage messsage
  else if (argc < 2 || (argc ==2 && ( !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) ) ) {
    PrintOptionsMessage();
    exit(1); 
  } 
#ifdef PARALLEL
  MPI::Init(argc, argv);
  Comms::Initialise();

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
  //const bool isFreqSampler = Comms::isFreqSampler();
  //  const bool isWorker = Comms::isWorker();

  //create results directory, or if it exists, deletes the contents
  if(isMaster){
    MakeResultsDir(options.getResultsDir().c_str(), false/*(options.getDisplayLevel()>2)*/);
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
      Log << On << "\nProgram aborted due to bad options.\n";
      exit(1);
    }
  
    //print user options to args.txt; must be done after all options are set
    if(isMaster)options.PrintOptions();

    AdmixMapModel M;

    M.Initialise(options, data, Log);

    data.Delete();

   //  ******** single individual, one population, fixed allele frequencies  ***************************
    if( M.getNumIndividuals() == 1 && options.getPopulations() == 1 && strlen(options.getAlleleFreqFilename()) )
      // nothing to do except calculate likelihood
      M.getOnePopOneIndLogLikelihood(Log, data.GetPopLabels());
    else {
      M.Run(options, data, Log);
    }

    if(isMaster){
      if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
      Log.ProcessingTime();
    }

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
  cout << "Rank " << MPI::COMM_WORLD.Get_rank() << " finished.\n";

  //MPI::COMM_WORLD.Barrier();
  Comms::Finalise();
  MPE_Finish_log("admixmap");
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

void PrintCopyrightNotice(LogWriter& Log){
  Log.setDisplayMode(On);
  cout << endl;
  Log << "-------------------------------------------------------\n"
      << "            ** ADMIXMAP (v" << ADMIXMAP_VERSION << "." << SUBVERSION
#ifdef PARALLEL
      << " (Parallel) "
#endif
      << ") **\n"
      << "-------------------------------------------------------\n";
  Log.setDisplayMode(Quiet);
  cout << "Copyright(c) 2002-2007 " << endl
       << "David O'Donnell, Clive Hoggart and Paul McKeigue" << endl
       << "-------------------------------------------------------"<<endl
       << "This program is free software distributed WITHOUT ANY WARRANTY " <<endl
       << "under the terms of the GNU General Public License. \nSee the file COPYING for details." <<endl
       << "-------------------------------------------------------" << endl;
}



