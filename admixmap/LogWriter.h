// *-*-C++-*-*
#ifndef LOGWRITER_H
#define LOGWRITER_H 1

#include <fstream>
#include <iostream>
#include <iomanip>

class LogWriter
{
private:
  std::ofstream *LogFileStreamPtr;
  int useCOUTOption;

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
  void write(int);
  void write(long);
  void write(double);

  void StartMessage(int, int, tm *);

  void StartOutput(int iteration, int TotalSamples);
  void EndOutput(int iteration);

  void Flush();
};

#endif /* !defined LATENT_H */
