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
}

LogWriter::LogWriter(const char *LogFilename, const bool useCout){
  LogFileStream.open(LogFilename, ios::out );
  if(!LogFileStream.is_open()){
    cerr << "ERROR: unable to open logfile"<<endl;
    exit(1);
  }
  useCOUTOption = useCout;
}

LogWriter::~LogWriter(){
  if(LogFileStream.is_open())LogFileStream.close();
}

void LogWriter::logmsg (const bool display, const string message)
{
  LogFileStream << message;
  if(useCOUTOption || display){
    cout << message;
  }
}

void LogWriter::logmsg (const bool display, const char * message)
{
  LogFileStream << message;
  if(useCOUTOption || display){
    cout << message;
  }
}

void LogWriter::logmsg (const bool display, const int number)
{
  LogFileStream << number;
  if(useCOUTOption  || display){
    cout << number;
  }
}

void LogWriter::logmsg (const bool display, const unsigned number)
{
  LogFileStream << number;
  if(useCOUTOption  || display){
    cout << number;
  }
}

void LogWriter::logmsg (const bool display, const long number)
{
  LogFileStream << number;
  if(useCOUTOption || display ){
    cout << number;
  }
}

void LogWriter::logmsg (const bool display, const double number)
{
  LogFileStream << number;
  if(useCOUTOption  || display){
    cout << number;
  }
}

void LogWriter::write(const char* message){
  LogFileStream<<message<<" ";
}
void LogWriter::write(const string message){
  LogFileStream<<message<<" ";
}
void LogWriter::write(const int number){
  LogFileStream<<number<<" ";
}
void LogWriter::write(const long number){
  LogFileStream<<number<<" ";
}
void LogWriter::write(const double number){
  LogFileStream<<number<<" ";
}
void LogWriter::write(const double number, const unsigned prec){
  LogFileStream<<setprecision(prec)<<number<<" ";
  LogFileStream<<setprecision(6);//restore default
}
void LogWriter::write(const double *array, const size_t dim){
  if(array)//to avoid seg faults with unallocated arrays
    for(size_t i = 0; i < dim;++i)LogFileStream<<array[i]<<" ";
}

void LogWriter::width(const unsigned w){
  LogFileStream.width(w);
}

void LogWriter::StartMessage(tm *timer){
  LogFileStream << "-----------------------------------------------" << endl;
  LogFileStream << "            ** ADMIXMAP (v" << ADMIXMAP_VERSION << ") **" << endl;
  LogFileStream << "-----------------------------------------------" << endl;
  logmsg(true,"Program started at ");
  logmsg(true,timer->tm_hour);
  logmsg(true,":");
  logmsg(true,timer->tm_min < 10 ? "0" : "" );
  logmsg(true,timer->tm_min);
  logmsg(true,".");
  logmsg(true,timer->tm_sec < 10 ? "0" : "" );
  logmsg(true,timer->tm_sec);
  logmsg(true," ");
  logmsg(true,timer->tm_mday);
  logmsg(true,"/");
  logmsg(true,timer->tm_mon+1);
  logmsg(true,"/");
  logmsg(true,1900+timer->tm_year);
  logmsg(true,"\n\n");
}
void LogWriter::Reset(const int iteration, const int width){
  if( !useCOUTOption || iteration == 0 )
    //output params to log on first iteration and every other when coutindicator = 0
    {
      LogFileStream << setiosflags( ios::fixed );
      LogFileStream.width( width );
      LogFileStream << iteration << " ";
    }

}

