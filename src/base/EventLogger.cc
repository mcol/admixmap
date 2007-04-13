#include "EventLogger.hh"
#include <sstream>
#include <cassert>
#ifndef _HAVE_MPE_
//define macros to replace MPE function calls
#define MPE_Init_log()
#define MPE_Finish_log(X)
#define MPE_Describe_state(M, N, NAME, DISPLAY)
#define MPE_Log_event(N, M, S)
#endif

unsigned EventLogger::NumEvents = 0;
bool EventLogger::isInitialised = false;

void EventLogger::Initialise(){
  MPE_Init_log();
  isInitialised = true;
}

void EventLogger::Finalise(const char* logname){
  assert(logname);
  if(isInitialised)MPE_Finish_log(logname);
}

void EventLogger::DefineEvent(int begin, int end, const char* eventName, const char* LineColour){
  assert( (end = begin+1) && eventName && LineColour);
  if(!isInitialised)Initialise();
  ++NumEvents;
  MPE_Describe_state(begin, end, eventName, LineColour);
}

void EventLogger::LogEvent(unsigned num, int ref, const char* message){
  assert(message && ref>=0);
  if(num > 2*NumEvents) {
    std::stringstream ss;
    ss <<"EventLogger Error: invalid event number (" << num << ")\n";
    throw(ss.str());
  }
  MPE_Log_event(num, ref, message);
}
