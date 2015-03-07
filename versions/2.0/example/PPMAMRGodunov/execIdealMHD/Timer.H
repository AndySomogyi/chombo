#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _TIMER_H_
#define _TIMER_H_

//#define CATFISH  // to compile stand-alone Timer tester

#ifdef CH_AIX
#define FOUR_COUNTERS
#endif

#include <cstdio>

#ifdef CH_MPI
#include "mpi.h"
#endif

#ifndef CATFISH
//#include "Counter.H"     // temp holding spot for this temp hack
#include "parstream.H"   // for pout()'s
#endif

#define DTIMESAMR    200
#define DTIMELS      202
#define DTIMEFORT    204
#define DTIMECOMM    206
#define DTIMEDEFINES 210
#define DTIMEALLGSRB 212

// if TIMER is not defined, then we just have a stub class.
#ifdef TIMER

extern "C" {
  // these are for gettimeofday() wall-clock timing
#include <unistd.h>
#include <sys/time.h>

  // obviously, if you don't have PAPI installed on your system
  // (which requires kernel modification) then the counters wont
  // work, so don't compile with -DPAPI
#ifdef PAPI
#include <papi.h>
#endif
}

#ifndef TIMER_COUNTER
#define TIMER_COUNTER 0
#endif

#include <list>
#include <string>
using namespace std;

class Timer {

public:

  Timer();                                  // unmanaged
  Timer(const string&, const int);          // root parent-only
  Timer(const string&, Timer&, const int);  // parent/child
  Timer(const string&, Timer&);             // child only
  Timer(const string&, Timer&, const int, bool);  // diagnostic

  //Timer(const string&, const int);  // counter

  ~Timer();

  void setup();
  void start();
  void stop();
  void clear();

  inline double getTimeStampWC();
  double wc_time() { return m_accumulated_WCtime; }

#ifdef PAPI
  long long int papi_counter0() { return m_accumulated_counter0; }
  long long int papi_counter1() { return m_accumulated_counter1; }
  double total_papi_counter0()  { return m_totalPapiCounter0; }
  double total_papi_counter1()  { return m_totalPapiCounter1; }
#ifdef FOUR_COUNTERS
  long long int papi_counter2() { return m_accumulated_counter2; }
  long long int papi_counter3() { return m_accumulated_counter3; }
  double total_papi_counter2()  { return m_totalPapiCounter2; }
  double total_papi_counter3()  { return m_totalPapiCounter3; }
#endif
#endif

  string        Name()    { return m_name;     }
  long long int Count()   { return m_count;    }
  int           tableID() { return m_tableID;  }

  static list<Timer*> *TimerList; // must be initialized somewhere

  static void TimerSummary(const int bert, const int ernie) {
    TimerSummary(); // backwards compatibilty...
  }

  static void TimerInit(const int rank);
  static void TimerSummary(void);

private:

  // internally called...
  static void TimerSummary_(void);
  static void writeParentTables(FILE *out, const double TimerCost);
  static void writeDiagnosticTables(FILE *out, const double TimerCost);

  //bool timer_on;    // State of timer, either on(true) or off(false)
  int m_tableID;  // the table ID -- where parent goes in the summary

  double m_accumulated_WCtime, m_last_WCtime_stamp;

  // wall-clock timer data
#ifndef CH_MPI
  struct timeval tv;   //  Values from call to gettimeofday
  struct timezone tz;
#endif

  int m_ID;
  string m_name;
  bool m_diagnostic;
  Timer& m_Parent;

  double m_avgWC, m_minWC, m_maxWC, m_avgCount;

#ifdef PAPI
  long long int m_accumulated_counter0;
  long long int m_accumulated_counter1;
  long long int m_previous_counter0;
  long long int m_previous_counter1;
  double m_totalPapiCounter0;
  double m_totalPapiCounter1;
#ifdef FOUR_COUNTERS
  long long int m_values[4];
  long long int m_accumulated_counter2;
  long long int m_accumulated_counter3;
  long long int m_previous_counter2;
  long long int m_previous_counter3;
  double m_totalPapiCounter2;
  double m_totalPapiCounter3;
#else
  long long int m_values[2];
#endif

#endif

  long long int m_count;
  long long int m_totalCount;
};

#else

// stub Timer class
#include <string>
#include <iostream>
using namespace std;

class Timer {

public:

  Timer(){}
  Timer(const string&, const int){}
  Timer(const string&, Timer&, const int){}
  Timer(const string&, Timer&){}
  Timer(const string&, Timer&, const int, bool){}

  ~Timer(){}

  void setup(){}
  void start(){}
  void stop(){}
  void clear(){}

  inline double getTimeStampWC(){return 0.0;}
  double wc_time() {return 0.0;}

#ifdef PAPI
  long long int papi_counter0() {return 0;}
  long long int papi_counter1() {return 0;}
  double total_papi_counter0()  {return 0.0;}
  double total_papi_counter1()  {return 0.0;}
#ifdef FOUR_COUNTERS
  long long int papi_counter2() {return 0;}
  long long int papi_counter3() {return 0;}
  double total_papi_counter2()  {return 0.0;}
  double total_papi_counter3()  {return 0.0;}
#endif
#endif

  string        Name()    { return NULL;   }
  long long int Count()   { return 0;   }
  int           tableID() { return 0;  }

  static void TimerSummary(const int bert, const int ernie) {
    TimerSummary(); // backwards compatibilty...
  }

  static void TimerInit(const int rank){}
  static void TimerSummary(void){}

private:

  // internally called...
  static void TimerSummary_(void){}
  static void writeParentTables(FILE *out, const double TimerCost){}
  static void writeDiagnosticTables(FILE *out, const double TimerCost){}
};

#endif // TIMER

#endif // TIMER_H
