//CLASS RESPONSIBLE FOR WRITING TO THE LOG FILE
#include "LogWriter.h"
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
void LogWriter::write(int number){
  *LogFileStreamPtr<<number;
}
void LogWriter::write(long number){
  *LogFileStreamPtr<<number;
}
void LogWriter::write(double number){
  *LogFileStreamPtr<<number;
}

void LogWriter::StartMessage(int TotalSamples, int BurnIn,tm *timer){
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
  logmsg(false,"\n");
  logmsg(false,"Total number of iterations: ");
  logmsg(false,TotalSamples);
  logmsg(false,"\n");
  logmsg(false,"Burnin: ");
  logmsg(false,BurnIn);
  logmsg(false,"\n");
   
}

void LogWriter::StartOutput(int iteration, int TotalSamples){
  //setup output streams
  if( !useCOUTOption || iteration == 0 )//do we really want output pars to log?
    {
      (*LogFileStreamPtr) << setiosflags( ios::fixed );
      LogFileStreamPtr->width( (int)( log10((double)TotalSamples)+1 ) );
      (*LogFileStreamPtr) << iteration << " ";
    }

  if( useCOUTOption )
    {
      cout << setiosflags( ios::fixed );
      (cout).width( (int)( log10((double)TotalSamples)+1 ) );
      cout << iteration << " ";
    }
}
void LogWriter::EndOutput(int iteration){
  //Print newline in output
  if( !useCOUTOption || iteration == 0 )
    {    *LogFileStreamPtr << endl;
    }
  if( useCOUTOption )
    {      cout << endl;
    }
}
void LogWriter::Flush(){
  *LogFileStreamPtr<<flush;
}
