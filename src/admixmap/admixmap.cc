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
using namespace genepi;



#define DEBUG_PRINT_EPROBS  0 // **** PRINT EMISSION PROBABILITIES DEBUG-CODE ****
#define max_print_states    12



#if DEBUG_PRINT_EPROBS

    #include "HiddenStateSpace.h"

    //-----------------------------------------------------------------------------
    /// Output a summary of a pedigree
    //-----------------------------------------------------------------------------

    ostream & ped_sum( ostream & os, const Pedigree & ped )
	{
	os << "pedigree \""
		    << ped.getId() << "\" (" << ped.getNMembers() << " members, "
		    << ped.getNFounders() << " founders, "
		    << ped.getNNonFndrs() << " non-founders)";
	return os;
	}

#endif



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

    // "check user options"
    // DDF: well, actually finish constructing them also.  They can't be used
    // prior to this call.
    if(options.checkOptions(Log, data.getNumberOfIndividuals())){
      Log << bclib::On << "\nProgram aborted due to bad options. See logfile for details\n";
      exit(1);
    }

    // DDF: What can one say?  This shouldn't be required, but circular
    // dependencies (AdmixOptions/InputAdmixData) do that to you.
    data.finishConstructing( options );

    //print user options to args.txt; must be done after all options are set
    options.PrintUserOptions("args.txt");

    AdmixMapModel M;

    M.Initialise(options, data, Log);

    // ---- DDF NOTE: ----
    // The old (pre-pedigree) method was to use the InputAdmixData (data) as a
    // temporary holding/conversion object which read the input-data from the
    // files, and used it to construct all of the necessary objects in the
    // AdmixMapModel (M).  The data object was then deleted (presumably just to
    // free memory) prior to beginning to run the model.  The Pedigree objects
    // violate the assumptions of this method in that they are more-or-less
    // constructed in-place when the files are parsed (being rather large,
    // complex, and unwieldy to copy), so either the model object must be passed
    // into the input-data object (presumably the best idea), which requires
    // that the method signatures of InputAdmixData (and InputData and
    // InputHapmixData) change, and that the model be constructed prior to the
    // input-data object; or that (currently the case) the input-data object
    // hold the "master" copies of these objects (Pedigree, SimpleLocusArray,
    // etc.), which requires that the input-data object remain in existence for
    // the lifetime of the model object (ugly).  A third option would be to copy
    // (take-ownership or use reference-counted objects would be preferable) the
    // objects from the input-data object to the model object.
    #if ! USE_GENOTYPE_PARSER
	data.Delete();
    #endif


    #if DEBUG_PRINT_EPROBS // **** PRINT EMISSION PROBABILITIES DEBUG-CODE ****

	const bool   print_eprobs	= true;
	const bool   print_ssize	= true;

	const cvector<Pedigree> & peds = const_cast<const InputAdmixData &>( data ) .getPeds();
	const SimpleLocusArray &  loci = data.getSimpleLoci();

	for ( cvector<Pedigree>::const_iterator iter = peds.begin() ; iter != peds.end() ; ++iter )
	    {

	    const Pedigree & ped = *iter;

	    cout << "\nPrinting hidden states and emission-probabilities for " << ped.getId() << ":\n";

	    if ( ped.getNFounders() > 1 )
		{
		cout << "\n  For this pedigree, AVs are listed in founder organism order:";

		for ( Pedigree::Iterator it = ped.getFirstFounder() ; it != ped.getEndFounder() ; ++it )
			cout << ' ' << (*it)->getOrgId();
		}

	    if ( ped.getNNonFndrs() > 1 )
		{
		cout << "\n"
			"  IVs are listed in non-founder organism order:";
		for ( Pedigree::Iterator it = ped.getFirstNonFndr() ; it != ped.getEndNonFndr() ; ++it )
			cout << ' ' << (*it)->getOrgId();
		}

	    if ( (ped.getNNonFndrs() > 1) || (ped.getNFounders() > 1) )
		cout << '\n';

	    for ( SLocIdxType sLocIdx = 0 ; sLocIdx < loci.size() ; ++sLocIdx )
		{
		const HiddenStateSpace & space = ped.getStateProbs( sLocIdx );

		if ( print_eprobs || print_ssize )
		    cout << (print_eprobs ? "\n" : "  ") << "Simple-locus \"" << loci[ sLocIdx ].getName() << "\":";

		size_t n_states = 0;

		for ( HiddenStateSpace::Iterator it(space); it; ++it )
		    {
		    ++n_states;
		    if ( print_eprobs && ((max_print_states == 0) || (n_states <= size_t(max_print_states))) )
			{
			const HiddenStateSpace::State & st = *it;
			cout << "\n  " << st.av << ' ' << st.iv << "  " << st.emProb;
			}
		    }

		if ( print_eprobs )
		    cout << "\n  ";

		if ( print_eprobs || print_ssize )
		    cout << ' ' << space.getNStates() << " states, of which "
			<< n_states << " have non-zero emission probability.\n";

		if ( print_eprobs )
		    cout << "\n\n";
		}

	    }

    #endif // DEBUG_PRINT_EPROBS



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



