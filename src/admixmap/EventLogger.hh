
///class to implement MPE event logging

#ifdef _HAVE_MPE_
#include <mpe.h>
#endif

class EventLogger{
public:

  static void Initialise();
  static void Finalise(const char* logname);
  static void DefineEvent(int, int, const char*, const char*);
  static void LogEvent(unsigned, int, const char* );

private:
  //not allowed to instantiate an object of this class
  EventLogger();
  ~EventLogger();

  static bool isInitialised;
  static unsigned NumEvents;

};
