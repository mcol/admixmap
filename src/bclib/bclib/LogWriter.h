// *-*-C++-*-*
/* 
 *   LogWriter.h
 *   Class for writing to a log file
 */
/*
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

#include "bclib/bclib.h"
#include <fstream>
#include <string>

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



/**
 * enum for three different display modes.
 * <BR>
 * EXAMPLES:
 *  -# to write to logfile only (logging messages):
 *	<PRE>
 *	setDisplayMode(Off); Log<<message;
 *	</PRE>
 *  -# to write to log and screen (important messages):
 *	<PRE>
 *	setDisplayMode(On); Log<<message;
 *	</PRE>
 *  -# to write to log and let user determine whether to write to screen (unimportant information)
 *	<PRE>
 *	setDisplayMode(Quiet); Log<<message;
 *	</PRE>
 */
enum DisplayMode
    {
    Off	    , ///< output only to logfile and not screen
    Quiet   , ///< output to logfile and if displaylevel>1 to screen
    On	      ///< output to logfile and to screen
    };

///Class to write to a logfile.
///Use insertion operator after setting appropriate display mode.
//Currently does not accept manipulators like flush, endl;
class LogWriter 
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
  LogWriter& operator<<(unsigned int);
  LogWriter& operator<<(const long);
  LogWriter& operator<<(unsigned long);
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
  std::ofstream file;
  bool verbose;///< determines if output goes to screen in quiet mode
  time_t StartTime;
  DisplayMode toscreen;
  int rank;

};

/** @} */

END_BCLIB_NAMESPACE

#endif 
