/*
 *   Copyright (c) 2006, 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

//=============================================================================
/// \file misc.cc
/// Miscellaneoud functions.
//=============================================================================

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "Model.h"
#include <dirent.h>//for OpenResultsDir
#ifdef HAVE_CONFIG_H 
#include "config.h"
#endif

/**
   Prints a message to screen in case of user not specifying any arguments
*/
void PrintUsage(const char* ProgName) {
  cout << "Usage: " << endl
       << "   " << ProgName << " <options-file>" << endl
       << "or " << ProgName << " -f<optionsfile> [extra options]" << endl
    //       << "or " << ProgName << "<options> (deprecated)" << endl
       << endl << "Consult the manual for details of user options." << endl;
}

/**
   Opens a directory for program output or if one exists, empties contents.
*/
void CreateDirectory(const char* dirname,  bool DeleteExistingFiles){
  if(strlen(dirname) && strcmp(dirname, ".") && strcmp(dirname, "..")){
    DIR *pdir;
    struct dirent *pent;
    
    string dirpath = "./";
    dirpath.append(dirname);
    
#ifdef WIN32
    //replace '/' with '\'
    string::size_type i = 0;
    while(i != string::npos){
      i = dirpath.find("/", i+1);
      dirpath[i] = '\\';
    }
    string cmd = "mkdir ";
#else
    string cmd = "mkdir -p ";
#endif
    
    pdir=opendir(dirpath.c_str()); 
    if (!pdir) {
      //dir does not exist
      cout << "Creating directory "<< dirpath << " with system command " << cmd <<endl;

      //KLUDGE: 'mkdir' not guaranteed to work on all systems; no error-checking
      //should be ok for normal use
      // -p = make recursively and no warning if dir exists

      cmd.append(dirpath);
      if ( system(cmd.c_str()) ) {;}; // "if" is to suppress compiler warning
      if(!opendir(dirpath.c_str())){
	cerr << "Unable to create directory. Exiting." << endl;
	exit(EXIT_FAILURE);
      }
      
    }
    else if(DeleteExistingFiles){
      cout << "Directory \"" << dirname << "\" exists. Contents will be deleted."<<endl;
      
      //list and delete contents of directory
      errno=0;
      while ((pent=readdir(pdir))){//read filenames
	if(strcmp(pent->d_name, ".") && strcmp(pent->d_name, "..")//skip . and ..
	   && strncmp(pent->d_name, "state", 5)){//skip any files starting 
	  string filepath = dirpath + "/"; 
	  filepath.append(pent->d_name);
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
    cerr << "Invalid directory name: " << dirname << " Exiting\n";
    exit(1);
  }
}



void ThrowException( const string & msg, bclib::LogWriter & Log )
    {
    Log << bclib::On << "\n\n"
	">>>FATAL ERROR: " << msg << "\n>>>ABORTING!\n\n";
    Log.ProcessingTime();
    exit(1);
    }



/**
   Writes an iteration number to screen, with newline for verbose output or as a counter for reduced output.
*/
void WriteIterationNumber(const int iteration, const int width, int displayLevel) {
  if( displayLevel > 2 ) {
    cout << '\n'
	 << std::setiosflags( ios::fixed ) //causes memory leak/overwrite error
	 << setw( width )
	 << iteration << ' ';
  }
  else if ( displayLevel > 1 && iteration > 0 ) { // display iteration counter only
    cout << "\rIterations so far: " << iteration;
  }

  cout.flush();
}


///print build info (compiler, flags, host type etc) to Log & screen
void PrintBuildInfo(bclib::LogWriter& Log){
  Log << bclib::On
#ifdef HAVE_CONFIG_H 
        << "Build Info:\n"
#ifdef _compiler_name
	<< "Compiler       = " << _compiler_name << "\n"
#endif
#ifdef _compiler_flags
	<< "Compiler Flags = " << _compiler_flags << "\n"
#endif
#ifdef _host
	<< "Host           = " << _host << "\n"
#endif
#ifdef _platform
	<< "Target         = " << _platform << "\n"
#endif
#if defined(_OPENMP)
	<< "OpenMP support = " << "enabled" << "\n"
#else
	<< "OpenMP support = " << "disabled" << "\n"
#endif
#ifdef _config_date
	<< "Config date    = " << _config_date << "\n"
#endif
        << "-------------------------------------------------------\n"
#endif
      << bclib::Quiet;
}


