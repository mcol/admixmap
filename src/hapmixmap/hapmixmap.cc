/** 
 *   HAPMIXMAP
 *   hapmixmap.cc
 *   Top-level source file
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "HapMixModel.h"
#include "EventLogger.hh"
#include <fstream>

#define HAPMIXMAP_VERSION 0
#define SUBVERSION 9

using namespace std;

int main( int argc , char** argv ){

  HapMixOptions options(argc, argv);

  //print version number and copyright info, if requested, and exit
  if(options.getFlag("version")){
    LogWriter LW;
    PrintCopyrightNotice(LW);
    exit(1);
  }

  //if no options specified or help requested, print help message and list of options, then exit
  if(!options.hasOptions() || options.getFlag("help")){
    LogWriter LW;
    //PrintCopyrightNotice(LW);
    PrintUsage("hapmixmap");
    options.PrintAllOptions(cout);
    exit(1);
  }

  if(!options.CheckRequiredOptions() || !options.SetOptions())
    exit(1);

#ifdef PARALLEL
  MPI::Init(argc, argv);
  Comms::Initialise();

#endif
  //define states for logging
  EventLogger::Initialise();
  if(Comms::isMaster()){
    EventLogger::DefineEvent(1, 2, "Barrier", "red:vlines1");
    EventLogger::DefineEvent(3, 4, "ReduceAncestry", "cyan:vlines3");
    EventLogger::DefineEvent(5, 6, "ReduceAlleleCounts", "orange:gray2");
    EventLogger::DefineEvent(7, 8, "SampleAlleleFreqs", "yellow:gray3");
    EventLogger::DefineEvent(9, 10, "SampleLambda", "LightBlue:gray3");
    EventLogger::DefineEvent(11, 12, "SetGenotypeProbs", "LightGreen:hlines2");
    EventLogger::DefineEvent(13, 14, "SampleHapPairs", "magenta:vlines2");
    EventLogger::DefineEvent(15, 16, "SampleAncestry", "blue:vlines3");
    EventLogger::DefineEvent(17, 18, "BcastLambda", "maroon:gray");
    EventLogger::DefineEvent(19, 20, "BcastFreqs", "plum:hlines3");
    EventLogger::DefineEvent(21, 22, "ScoreTests", "gray:hlines3");
    EventLogger::DefineEvent(23, 24, "UpdateAlleleCounts", "DarkGreen:hlines4");
  }

  const bool isMaster = Comms::isMaster();
  //const bool isFreqSampler = Comms::isFreqSampler();
  //  const bool isWorker = Comms::isWorker();

  //create results directory, or if it exists, deletes the contents
  if(isMaster){
    MakeResultsDir(options.getResultsDir().c_str(), false/*(options.getDisplayLevel()>2)*/, options.getDeleteOldResultsIndicator());
  }
 
  //open logfile, start timer and print start message
  LogWriter Log(options.getLogFilename(), (bool)(options.getDisplayLevel()>1 && isMaster));
  if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);

  if(isMaster){
    //if(options.getDisplayLevel()>0 )
   PrintCopyrightNotice(Log);
   if(options.getFlag("checkmode"))
     Log << On << "  *** Check Mode Active *** \n"
	 << "-------------------------------------------------------\n";
   else{
     if(options.doPrintBuildInfo())PrintBuildInfo(Log);
     Log.StartMessage();
   }
  }

  try{  
    Rand RNG;//allocate random number generator
    RNG.setSeed( options.getSeed() );  // set random number seed
  
    //read data files and check (except allelefreq files)
    //also sets 'numberofoutcomes' and 'states' options
    InputHapMixData data(&options, Log);

     //check user options
    if(options.checkOptions(Log, data.getNumberOfIndividuals())){
      Log << On << "\nProgram aborted due to bad options. See logfile for details\n";
      exit(1);
    }
  
    //print user options to args.txt; must be done after all options are set
    if(isMaster)options.PrintUserOptions("args.txt");

    //end of program, in checkmode
    if(options.getFlag("checkmode")){
      Log << On << "-------------------------------------------------------\n"
	  <<  "  *** Everything looks good ***\n\n  *** Check Mode Complete ***\n" 
	  << "-------------------------------------------------------\n";
      exit(0);
    }

    HapMixModel M;
    M.Initialise(options, data, Log);

    data.Delete();
 
    M.Run(options, data, Log, options.getNumAnnealedRuns());

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

  EventLogger::Finalise("hapmixmap");
  //print run times to screen and log
  if(isMaster){
    if(options.getDisplayLevel()==0)Log.setDisplayMode(Off);
    else  Log << "\n";
    Log.ProcessingTime();
  }
    
#ifdef PARALLEL
  cout << "Rank " << MPI::COMM_WORLD.Get_rank() << " finished.\n";
  //MPI::COMM_WORLD.Barrier();
  Comms::Finalise();
  MPI_Finalize();
#else
  cout << "Finished" << endl;
#endif

    cout << "-------------------------------------------------------" <<endl;

  putenv("HAPMIXMAPCLEANEXIT=1");
  return 0;
} //end of main

void PrintCopyrightNotice(LogWriter& Log){
  Log.setDisplayMode(On);
  cout << endl;
  Log << "-------------------------------------------------------\n"
      << "            ** HAPMIXMAP (v" << HAPMIXMAP_VERSION << "." << SUBVERSION
#ifdef PARALLEL
      << " (Parallel) "
#endif
      << ") **\n"
      << "-------------------------------------------------------\n";
  Log.setDisplayMode(Quiet);
  cout << "Copyright(c) 2006, 2007 " << endl
       << "David O'Donnell and Paul McKeigue" << endl
       << "-------------------------------------------------------"<< endl
       << "This program is free software distributed WITHOUT ANY WARRANTY" << endl
       << "under the terms of the GNU General Public License."<< endl
       << "See the file COPYING for details." << endl
       << "-------------------------------------------------------" << endl;
}

