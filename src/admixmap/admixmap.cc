/** 
 *   ADMIXMAP
 *   \file admixmap.cc
 *   Top-level source file, containing the main() program for admixmap
 */
/*
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "AdmixMapModel.h"
#include "config.h" // USE_GENOTYPE_PARSER
#if USE_GENOTYPE_PARSER
    #include "DataValidError.h"
#endif
#include <fstream>

#define ADMIXMAP_VERSION "3"
#define SUBVERSION	 "8."

extern const char SVN_VERSION [];

using namespace std;

int main( int argc , char** argv ){
  AdmixOptions options(argc, argv);

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
    PrintUsage("admixmap");
    //options.PrintAllOptions(cout);
    exit(1);
  }
  if(!options.CheckRequiredOptions() || !options.SetOptions())
    exit(1);

  // ******************* PRIMARY INITIALIZATION ********************************************************************************

  //create results directory, or if it exists, deletes the contents
  CreateDirectory(options.getResultsDir().c_str());
 
  //open logfile, start timer and print start message
  bclib::LogWriter Log(options.getLogFilename(), (bool)(options.getDisplayLevel()>1));
  if(options.getDisplayLevel()==0)Log.setDisplayMode(bclib::Off);
 
  //if(options.getDisplayLevel()>0 )
  PrintCopyrightNotice(Log);
  if(options.getFlag("printbuildinfo"))PrintBuildInfo(Log);
  Log.StartMessage();
 

  try{  
    bclib::Rand RNG;//allocate random number generator
    RNG.setSeed( options.getSeed() );  // set random number seed
  
    //read data files and check
    //also sets 'numberofoutcomes' and 'populations' options
    InputAdmixData data(&options, Log); 

     //check user options
    if(options.checkOptions(Log, data.getNumberOfIndividuals())){
      Log << bclib::On << "\nProgram aborted due to bad options. See logfile for details\n";
      exit(1);
    }
  
    //print user options to args.txt; must be done after all options are set
    options.PrintUserOptions("args.txt");

    AdmixMapModel M;

    M.Initialise(options, data, Log);

    data.Delete();

   //  ******** single individual, one population, fixed allele frequencies  ***************************
    if( M.getNumIndividuals() == 1 && options.getPopulations() == 1 && strlen(options.getAlleleFreqFilename()) )
      // nothing to do except calculate likelihood
      M.getOnePopOneIndLogLikelihood(Log, data.GetPopLabels());
    else {
      if(options.getTestOneIndivIndicator()) { 
	M.TestIndivRun(options, data, Log);
      } 
      else
	M.Run(options, data, Log);
    }

    if(options.getDisplayLevel()==0)Log.setDisplayMode(bclib::Off);
    Log.ProcessingTime();

  }

#if USE_GENOTYPE_PARSER
  // Catch data-validation errors and format nicely:
  catch ( genepi::DataValidError & e )
    {
    // options->getProgramName()
    cerr << '\n' << argv[0] << ": ***** input data validation error *****"
	<< (genepi::DataValidError::getEmacsStyleFmt() ? "\n  " : " ")
	<< e.what() << "\n\n";
    return 1;
    }
#endif

  catch (const string& msg) {//catch any stray error messages thrown upwards
    ThrowException(msg, Log);
  }
  catch (const char* msg) {//in case error messages thrown as char arrays instead of strings
    ThrowException(string(msg), Log);
  }
  catch (const exception& e){
    ThrowException(e.what(), Log);
  }

  catch(...){
    cout << "Unknown exception occurred. Contact the program authors for assistance" << endl;
    exit(1);
  }
  cout << "Finished" << endl
       << setfill('*') << setw(80) << "*" <<endl;

  putenv((char*)"ADMIXMAPCLEANEXIT=1");
  return 0;
} //end of main

void PrintCopyrightNotice(bclib::LogWriter& Log){
  Log.setDisplayMode(bclib::On);
  cout << endl;
  Log << "-------------------------------------------------------\n"
      << "            ** ADMIXMAP (v" ADMIXMAP_VERSION "." SUBVERSION << SVN_VERSION << ") **\n"
      << "-------------------------------------------------------\n";
  Log.setDisplayMode(bclib::Quiet);
  cout << "Copyright(c) 2002-2007 " << endl
       << "David O'Donnell, Clive Hoggart and Paul McKeigue" << endl
       << "-------------------------------------------------------"<<endl
       << "This program is free software distributed WITHOUT ANY WARRANTY " <<endl
       << "under the terms of the GNU General Public License. \nSee the file COPYING for details." <<endl
       << "-------------------------------------------------------" << endl;
}



