// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   LogWriter.h
 *   Header file for LogWriter class
 *   all writing to the logfile is done via this class
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
#ifndef LOGWRITER_H
#define LOGWRITER_H 1

#include <fstream>
#include <iostream>
#include <iomanip>

enum DisplayMode {Off, Quiet, On};
//set DisplayMode to Off to output only to logfile and not screen
//                   Quiet        to logfile and if displaylevel>1 to screen
//                   On            to logfile and to screen
//EXAMPLES:
//(1) to write to logfile only (logging messages):
//setDisplayMode(Off); Log<<message;
//(2) to write to log and screen (important messages):
//setDisplayMode(On); Log<<message;
// (3) to write to log and let user determine whether to write to screen (unimportant information)
//setDisplayMode(Quiet); Log<<message;

class LogWriter
{
public:
  LogWriter();
  LogWriter(const char *LogFilename, const bool isverbose);
  ~LogWriter();

  void setDisplayMode(DisplayMode);

  LogWriter& operator<<(const int);
  LogWriter& operator<<(const unsigned);
  LogWriter& operator<<(const long);
  LogWriter& operator<<(const double);
  LogWriter& operator<<(const std::string);
  LogWriter& operator<<(const char*);

  void width(const unsigned w);//calls width function for LogFileStream
  void setPrecision(int);

  void StartMessage();//prints startup message
  void ProcessingTime();//prints finish time and time elapsed

private:
  std::ofstream LogFileStream;
  bool verbose;
  long StartTime;
  DisplayMode toscreen;
};

#endif /* !defined LATENT_H */
