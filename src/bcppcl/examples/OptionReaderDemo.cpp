/*
 * Here is a sample of how to use OptionReader to
 * parse comand line argumnets and an ptions file
 *
 * Create  sample.txt as follows
 *
 *      # sample options file
 *      # this is a comment 
 *      zip  
 *      size  : 42 
 *      title : This is a test title.
 *
 * Run the sample with any combination of the options
 *
 *      a.out -c --zip -s 20 --name foo.jpg argtest1 argtest2 
 */

#include <iostream>
#include "bcppcl/OptionReader.h"

void example( int argc, char* argv[] );
void PrintUsage();

using namespace::std;
using namespace::bcppcl;

int main( int argc, char* argv[] ){

  example( argc, argv ); 
  return 0 ;
}

void example( int argc, char* argv[] ){

  /* 1. CREATE AN OBJECT */
  OptionReader opt;

  /* 2. SET PREFERENCES  */
  opt.setVerbose(true); /* print warnings about unknown options */

  /* 3. SET THE OPTIONS STRINGS/CHARACTERS */
  int size = 0;
  string Name;

  opt.addFlag(  'h', "help" );   /* a flag (takes no argument), supporting long and short form */ 
  opt.addFlag( 'c' );            /* a flag (takes no argument), supporting only short form */
  opt.addOption( 's', "size", intOption, &size); /* an option (takes an argument), supporting long and short form */
  opt.addOption(  "name", stringOption, &Name, true ); /* a required option (takes an argument), supporting only long form */

  /* 4A. PROCESS THE COMMANDLINE AND RESOURCE FILE */
  /* read options from a  option/resource file with ':' separated options or flags, one per line */
  //bool allClear = opt.ReadArgsFromFile( argv[0] );  

  /* go through the command line and get the options  */
  //bool allClear = opt.ReadCommandLineArgs( argc, argv );
  /* supply an extra parameter to provide a prefix for option file argument. Then you can specify options both on the command-line and in file*/
  //bool allClear = opt.ReadCommandLineArgs( argc, argv, "-f" );

  /* 4B. PROCESS BOTH COMMANDLINE AND RESOURCE FILE */
  // autodetect - a single argument is assumed to be a filename, otherwise command-line options are read.
  // option values are then set, for convenience
  bool allClear = opt.ReadUserOptions(argc, argv);

  /* print usage if no options or if requested*/
  if( (! opt.hasOptions()) || opt.getFlag( "help" ) || opt.getFlag( 'h' )) { 
    //PrintUsage();
    opt.PrintAllOptions(cout);
    return;
  }

  /* 5A. CHECK REQUIRED OPTIONS AND SET THE VALUES */
  //use this if you used ReadCommandLineArgs or ReadArgsFromFile above
  //if(!opt.CheckRequiredOptions() | !opt.SetOptions() )
  //return;

  /* 5B. CHECK REQUIRED OPTIONS  */
  //use this if you used ReadUserOptions above (values have already been set)
  if(!opt.CheckRequiredOptions() || !allClear )
    return;


  /* 6. PRINT THE VALUES */
  cout << "size = " << size << endl ;
  cout << "name = " << Name << endl ;
  if( opt.getFlag( 'c' ) )  
    cout << "c = flag set " << endl ;
  cout << endl ;
  
  /* 7. CLEAN UP */
  opt.clear();
  
}

 /*  PRINT THE USAGE/HELP   */
void PrintUsage(){

  cout 
    << "Usage: "  << endl
    <<" -h  --help  		Prints this help " << endl
    << " -s  --size 42 	        Image Size " << endl 
    << " -c   			convert Image " << endl
    << "     --name image.jpg	Image Name " << endl;
}
