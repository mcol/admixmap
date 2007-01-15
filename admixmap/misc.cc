/** 
 *   misc.cc 
 *   miscellaneous functions
 *   Copyright (c) 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "admixmap.h"
#include "parallel.h"//for PARALLEL
#include <dirent.h>//for OpenResultsDir

/**
   Prints a message to screen in case of user not specifying any arguments
*/
void PrintOptionsMessage() {
  cout << "You must specify an options file\n"
       << "Consult the manual for a list of user options."
       << endl;
}

/**
   Opens a directory for program output
*/
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

/**
   Writes an iteration number to screen, with newline for verbose output or as a counter for reduced output.
*/
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



