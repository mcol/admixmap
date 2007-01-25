/*
 *   LogWriter.cc 
 *   Class for writing to a log file
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "LogWriter.h"
#include "parallel.h"
#ifdef PARALLEL
#include <mpi++.h>
#endif

using namespace::std;

struct nullstream:
  std::ostream {
  nullstream(): std::ios(0), std::ostream(0) {}
};

LogWriter::LogWriter(){
  toscreen = Off;
  verbose = false;
}

LogWriter::LogWriter(const char *LogFilename, const bool isverbose){
    rank = 0;
#ifdef PARALLEL 
    rank = MPI::COMM_WORLD.Get_rank();
#endif
    if(rank==0){ 
      LogFileStream.open(LogFilename, ios::out );
      if(!LogFileStream.is_open()){
	cerr << "ERROR: unable to open logfile"<<endl;
	exit(1);
      }
    }
    verbose = isverbose;
    toscreen = On;
}

LogWriter::~LogWriter(){
  if(LogFileStream.is_open())LogFileStream.close();
}

/**
 * Set whether the messages should be displayed.
 * Accepted values: On, Quiet, Off.
 */
void LogWriter::setDisplayMode(DisplayMode d){
  toscreen = d;
}

/**
 * Overloaded stream insertion operators,
 * write to log and screen, unless DisplayMode
 * switched off
 */
LogWriter& LogWriter::operator<<(const int message){
  LogFileStream << message;
  if(toscreen==On || (verbose && toscreen==Quiet)){
    cout << message;
  }
  return *this;
}
LogWriter& LogWriter::operator<<(const unsigned message){
  LogFileStream << message;
  if(toscreen==On || (verbose && toscreen==Quiet)){
    cout << message;
  }
  return *this;
}
LogWriter& LogWriter::operator<<(const long message){
  LogFileStream << message;
  if(toscreen==On || (verbose && toscreen==Quiet)){
    cout << message;
  }
  return *this;
}
LogWriter& LogWriter::operator<<(const double message){
  LogFileStream << message;
  if(toscreen==On || (verbose && toscreen==Quiet)){
    cout << message;
  }
  return *this;
}
LogWriter& LogWriter::operator<<(const string message){
  LogFileStream << message << flush;
  if(toscreen==On || (verbose && toscreen==Quiet)){
    cout << message << flush;
  }
  return *this;
}
LogWriter& LogWriter::operator<<(const char* message){
  LogFileStream << message << flush;
  if(toscreen==On || (verbose && toscreen==Quiet)){
    cout << message << flush;
  }
  return *this;
}

LogWriter& LogWriter::operator<<(DisplayMode d){
  toscreen = d;
  return *this;
}

void LogWriter::width(const unsigned w){
  if(rank==0)LogFileStream.width(w);
}
void LogWriter::setPrecision(int p){
  if(rank==0){
    LogFileStream<<setprecision(p);
    cout<<setprecision(p);
  }
}

void LogWriter::StartMessage(){
  if(rank==0){
    //start timer
    //#ifdef PARALLEL
    //StartTime = MPI::Wtime();
    //#else
    StartTime = time(0);
    //#endif
    tm timer = *localtime( &StartTime );
    
    toscreen = On;
    
#ifdef PARALLEL
    *this << "Running on " << MPI::COMM_WORLD.Get_size() << " processors\n";
#endif
    *this << "Program started at "
	  << timer.tm_hour << ":" << (timer.tm_min < 10 ? "0" : "")  << timer.tm_min << "." << (timer.tm_sec < 10 ? "0" : "") 
	  << timer.tm_sec << " " << timer.tm_mday << "/" << timer.tm_mon+1 << "/" << 1900+timer.tm_year << "\n";
  }
}

void LogWriter::ProcessingTime()
{
  if(rank==0){
    long EndTime = time(0);
    tm timer;
    timer = *localtime( &EndTime );
    
    toscreen = On;
    *this << "Program finished at " << timer.tm_hour << ":" << (timer.tm_min < 10 ? "0" : "")
	  << timer.tm_min << "." << (timer.tm_sec < 10 ? "0" : "")  << timer.tm_sec << " "
	  << timer.tm_mday << "/" << timer.tm_mon+1 << "/" << 1900+timer.tm_year << "\n";
    
    double realtime = difftime(EndTime, StartTime);
    *this << "Elapsed seconds = " << (int)realtime << "\n";
    *this << "Elapsed time = ";
    if(realtime > 3600.0){
      int hours = (int)(realtime/3600);
      *this << hours << "h, ";
      realtime -= (double)(hours*3600);
    }
    //if(realtime > 60.0){
    int mins = (int)(realtime/60);
    *this <<  mins << "m, ";
    realtime -= (double)(mins*60);
    //}
    
    *this << (int)realtime << "s\n";
  }
}

