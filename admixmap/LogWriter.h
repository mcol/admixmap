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

class LogWriter
{
public:
  LogWriter();
  LogWriter(const char *LogFilename, const bool useCout);
  ~LogWriter();

  void Reset(const int iteration, const int);

  //logmsg functions write to logfile
  //also write to screen unless cout = 0 and first arg is false
  void logmsg(const bool , const std::string message);

  void logmsg(const bool , const char * message);

  void logmsg(const bool , const int number);

  void logmsg(const bool , const unsigned number);

  void logmsg(const bool , const long number);

  void logmsg(const bool, const double number);

  //write functions write only to log with a space at end
  void write(const char*);
  void write (const std::string message);
  void write(const int);
  void write(const long);
  void write(const double);
  void write(const double number, const unsigned prec);
  void write(const double *array, const size_t dim);

  void width(const unsigned w);//calls width function for LogFileStream

  void StartMessage(tm *);//prints startup message

private:
  std::ofstream LogFileStream;
  bool useCOUTOption;
};

#endif /* !defined LATENT_H */
