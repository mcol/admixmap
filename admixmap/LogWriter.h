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

  void Reset(int iteration, bool, int);

  //logmsg functions write to logfile
  //also write to screen unless cout = 0 and first arg is false
  void logmsg(bool , std::string message);

  void logmsg(bool , const char * message);

  void logmsg(bool , int number);

  void logmsg(bool , unsigned int number);

  void logmsg(bool , long number);

  void logmsg(bool, double number);

  //write functions write only to log with a space at end
  void write(const char*);
  void write (std::string message);
  void write(int);
  void write(long);
  void write(double);
  void write(double number, unsigned prec);
  void write(double *array, size_t dim);

  void width(unsigned w);//calls width function for LogFileStream

  void StartMessage(tm *);//prints startup message

private:
  std::ofstream LogFileStream;
  bool useCOUTOption;
};

#endif /* !defined LATENT_H */
