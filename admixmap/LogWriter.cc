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

void LogWriter::Initialise(std::ofstream *LogFileStream, int useCout){
  //LogFileStreamPtr.open(LogFilename(), ios::out );
  LogFileStreamPtr = LogFileStream;
  useCOUTOption=useCout;
}

void
LogWriter::logmsg (bool useCOUT, std::string message)
{
  *LogFileStreamPtr << message;
  if(useCOUTOption || useCOUT){
    cout << message;
  }
}

void
LogWriter::logmsg (bool useCOUT, const char * message)
{
  *LogFileStreamPtr << message;
  if(useCOUTOption || useCOUT){
    cout << message;
  }
}

void
LogWriter::logmsg (bool useCOUT, int number)
{
  *LogFileStreamPtr << number;
  if(useCOUTOption  || useCOUT){
    cout << number;
  }
}

void
LogWriter::logmsg (bool useCOUT, unsigned int number)
{
  *LogFileStreamPtr << number;
  if(useCOUTOption  || useCOUT){
    cout << number;
  }
}

void
LogWriter::logmsg (bool useCOUT, long number)
{
  *LogFileStreamPtr << number;
  if(useCOUTOption || useCOUT ){
    cout << number;
  }
}

void
LogWriter::logmsg (bool useCOUT, double number)
{
  *LogFileStreamPtr << number;
  if(useCOUTOption  || useCOUT){
    cout << number;
  }
}

void LogWriter::write(const char* message){
  *LogFileStreamPtr<<message;
}
void LogWriter::write(std::string message){
  *LogFileStreamPtr<<message;
}
void LogWriter::write(int number){
  *LogFileStreamPtr<<number;
}
void LogWriter::write(long number){
  *LogFileStreamPtr<<number;
}
void LogWriter::write(double number){
  *LogFileStreamPtr<<number;
}

void LogWriter::StartMessage(tm *timer){
  *LogFileStreamPtr << "-----------------------------------------------" << endl;
  *LogFileStreamPtr << "            ** ADMIXMAP (v" << ADMIXMAP_VERSION << ") **" << endl;
  *LogFileStreamPtr << "-----------------------------------------------" << endl;
  logmsg(false,"Program started at ");
  logmsg(false,timer->tm_hour);
  logmsg(false,":");
  logmsg(false,timer->tm_min < 10 ? "0" : "" );
  logmsg(false,timer->tm_min);
  logmsg(false,".");
  logmsg(false,timer->tm_sec < 10 ? "0" : "" );
  logmsg(false,timer->tm_sec);
  logmsg(false," ");
  logmsg(false,timer->tm_mday);
  logmsg(false,"/");
  logmsg(false,timer->tm_mon+1);
  logmsg(false,"/");
  logmsg(false,1900+timer->tm_year);
  logmsg(false,"\n\n");
}


