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
#include <fstream>

#define HAPMIXMAP_VERSION 1
#define SUBVERSION 0

using namespace std;

int main( int argc , char** argv ){

  HapMixOptions options(argc, argv);

  //print version number and copyright info, if requested, and exit
  if(options.getFlag("version")){
    bclib::LogWriter LW;
    PrintCopyrightNotice(LW);
    exit(1);
  }

  //if no options specified or help requested, print help message and list of options, then exit
  if(!options.hasOptions() || options.getFlag("help")){
    bclib::LogWriter LW;
    //PrintCopyrightNotice(LW);
    PrintUsage("hapmixmap");
    options.PrintAllOptions(cout);
    exit(1);
  }

  if(!options.CheckRequiredOptions() || !options.SetOptions())
    exit(1);

  //create results directory, or if it exists, deletes the contents
  CreateDirectory(options.getResultsDir().c_str(), options.getDeleteOldResultsIndicator());

  //open logfile, start timer and print start message
  bclib::LogWriter Log(options.getLogFilename(), (bool)(options.getDisplayLevel()>1));
  if(options.getDisplayLevel()==0)Log.setDisplayMode(bclib::Off);

  //if(options.getDisplayLevel()>0 )
  PrintCopyrightNotice(Log);
  if(options.getFlag("checkmode"))
    Log << bclib::On << "  *** Check Mode Active *** \n"
	<< "-------------------------------------------------------\n";
  else{
    if(options.getFlag("printbuildinfo"))PrintBuildInfo(Log);
    Log.StartMessage();
  }

  try{  
    bclib::Rand RNG;//allocate random number generator
    RNG.setSeed( options.getSeed() );  // set random number seed
  
    //read data files and check (except allelefreq files)
    //also sets 'numberofoutcomes' and 'states' options
    InputHapMixData data(&options, Log);

     //check user options
    if(options.checkOptions(Log, data.getNumberOfIndividuals())){
      Log << bclib::On << "\nProgram aborted due to bad options. See logfile for details\n";
      exit(1);
    }
  
    //print user options to args.txt; must be done after all options are set
    options.PrintUserOptions("args.txt");

    //end of program, in checkmode
    if(options.getFlag("checkmode")){
      Log << bclib::On << "-------------------------------------------------------\n"
	  <<  "  *** Everything looks good ***\n\n  *** Check Mode Complete ***\n" 
	  << "-------------------------------------------------------\n";
      exit(0);
    }

    //create directory to write final parameter values, unless it is the same as resultsdir
    const string& finalvaluedir = options.getFinalValueDir();
    if(finalvaluedir.compare(options.getResultsDir()))
      CreateDirectory(finalvaluedir.c_str(), false);


    HapMixModel M;
    M.Initialise(options, data, Log);

    data.Delete();
 
    M.Run(options, data, Log);

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

  catch(...){
    cout << "Unknown exception occurred. Contact the program authors for assistance" << endl;
    exit(1);
  }

  //print run times to screen and log
  if(options.getDisplayLevel()==0)Log.setDisplayMode(bclib::Off);
  else  Log << "\n";
  Log.ProcessingTime();
    
  cout << "Finished" << endl
       << "-------------------------------------------------------" << endl;

  putenv("HAPMIXMAPCLEANEXIT=1");
  return 0;
} //end of main

void PrintCopyrightNotice(bclib::LogWriter& Log){
  Log.setDisplayMode(bclib::On);
  cout << endl;
  Log << "-------------------------------------------------------\n"
      << "            ** HAPMIXMAP (v" << HAPMIXMAP_VERSION << "." << SUBVERSION << ") **\n"
      << "-------------------------------------------------------\n";
  Log.setDisplayMode(bclib::Quiet);
  cout << "Copyright(c) 2006, 2007 " << endl
       << "David O'Donnell and Paul McKeigue" << endl
       << "-------------------------------------------------------"<< endl
       << "This program is free software distributed WITHOUT ANY WARRANTY" << endl
       << "under the terms of the GNU General Public License."<< endl
       << "See the file COPYING for details." << endl
       << "-------------------------------------------------------" << endl;
}

