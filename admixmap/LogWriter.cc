/** 
 *   ADMIXMAP
 *   LogWriter.cc 
 *   Class responsible for writing to the log file
 *   Copyright (c) 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#include "LogWriter.h"
#include "admixmap.h"

using namespace::std;

LogWriter::LogWriter(){
  toscreen = On;
}

LogWriter::LogWriter(const char *LogFilename, const bool isverbose){
  LogFileStream.open(LogFilename, ios::out );
  if(!LogFileStream.is_open()){
    cerr << "ERROR: unable to open logfile"<<endl;
    exit(1);
  }
  verbose = isverbose;
  toscreen = On;
}

LogWriter::~LogWriter(){
  if(LogFileStream.is_open())LogFileStream.close();
}

void LogWriter::setDisplayMode(DisplayMode d){
  toscreen = d;
}

// ** Overloaded stream insertion operators, write to log and screen, unless DisplayMode switched off
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
  LogFileStream << message;
  if(toscreen==On || (verbose && toscreen==Quiet)){
    cout << message;
  }
  return *this;
}
LogWriter& LogWriter::operator<<(const char* message){
  LogFileStream << message;
  if(toscreen==On || (verbose && toscreen==Quiet)){
    cout << message;
  }
  return *this;
}

void LogWriter::width(const unsigned w){
  LogFileStream.width(w);
}
void LogWriter::setPrecision(int p){
  LogFileStream<<setprecision(p);
  cout<<setprecision(p);
}

void LogWriter::StartMessage(){
  //start timer
  StartTime = time(0);
  tm timer = *localtime( &StartTime );

  toscreen = On;
  LogFileStream << "-----------------------------------------------" << endl;
  LogFileStream << "            ** ADMIXMAP (v" << ADMIXMAP_VERSION << ") **" << endl;
  LogFileStream << "-----------------------------------------------" << endl;
  *this << "Program started at "
	<< timer.tm_hour << ":" << (timer.tm_min < 10 ? "0" : "")  << timer.tm_min << "." << (timer.tm_sec < 10 ? "0" : "") 
	<< timer.tm_sec << " " << timer.tm_mday << "/" << timer.tm_mon+1 << "/" << 1900+timer.tm_year << "\n\n";
}

void LogWriter::ProcessingTime()
{
  long EndTime = time(0);
  tm timer;
  timer = *localtime( &EndTime );

  toscreen = On;
  *this << "\nProgram finished at " << timer.tm_hour << ":" << (timer.tm_min < 10 ? "0" : "")
	<< timer.tm_min << "." << (timer.tm_sec < 10 ? "0" : "")  << timer.tm_sec << " "
	<< timer.tm_mday << "/" << timer.tm_mon+1 << "/" << 1900+timer.tm_year << "\n";

  double realtime = difftime(EndTime, StartTime);
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

