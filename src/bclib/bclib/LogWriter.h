// *-*-C++-*-*
/* 
 *   LogWriter.h
 *   Class for writing to a log file
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef LOGWRITER_H
#define LOGWRITER_H 1

#include "FileWriter.h"

BEGIN_BCLIB_NAMESPACE

/**
   enum for three different display modes.
   set DisplayMode to Off to output only to logfile and not screen
   Quiet        to logfile and if displaylevel>1 to screen
   On            to logfile and to screen
   EXAMPLES:
   (1) to write to logfile only (logging messages):
   setDisplayMode(Off); Log<<message;
   (2) to write to log and screen (important messages):
   setDisplayMode(On); Log<<message;
   (3) to write to log and let user determine whether to write to screen (unimportant information)
   setDisplayMode(Quiet); Log<<message;
*/
enum DisplayMode {Off, Quiet, On};

///Class to write to a logfile.
///Use insertion operator after setting appropriate display mode.
//Currently does not accept manipulators like flush, endl;
class LogWriter : public FileWriter
{
public:
  ///constructor - supply filename and indicate behaviour for quiet mode
  LogWriter(const char *LogFilename, const bool isverbose);
  void open(const char *LogFilename, const bool isverbose = false);
  ///default constructor, makes verbose output, no logfile
  LogWriter();
  ~LogWriter();

  void setDisplayMode(DisplayMode);

  LogWriter& operator<<(const int);
  LogWriter& operator<<(const unsigned);
  LogWriter& operator<<(const long);
  LogWriter& operator<<(const double);
  LogWriter& operator<<(const std::string);
  LogWriter& operator<<(const char*);
  LogWriter& operator<<(DisplayMode d);
  ///calls width function for LogFileStream
  void width(const unsigned w);
  //set number of digits in output
  void setPrecision(int);
  ///prints startup message
  void StartMessage();
  ///prints finish time and time elapsed
  void ProcessingTime();

private:
  bool verbose;///< determines if output goes to screen in quiet mode
  time_t StartTime;
  DisplayMode toscreen;
  int rank;

};
END_BCLIB_NAMESPACE
#endif 
