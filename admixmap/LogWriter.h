// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   LogWriter.h
 *   Header file for LogWriter class
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
  void Initialise(std::ofstream *LogFileStream, int useCout);

  void logmsg(bool useCOUT, std::string message);

  void logmsg(bool useCOUT, const char * message);

  void logmsg(bool useCOUT, int number);

  void logmsg(bool useCOUT, unsigned int number);

  void logmsg(bool useCOUT, long number);

  void logmsg(bool useCOUT, double number);

  void write(const char*);
  void write (std::string message);
  void write(int);
  void write(long);
  void write(double);

  void StartMessage(tm *);

private:
  std::ofstream *LogFileStreamPtr;
  int useCOUTOption;
};

#endif /* !defined LATENT_H */
