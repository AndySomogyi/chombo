#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifdef TIMER

#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <cstdio>
#include "SPACE.H"
using namespace std;

#include "Timer.H"

#ifndef CATFISH
#include "parstream.H"
#endif

#ifdef MEMORY_USAGE
#include "memusage.H"
static char stuff[180];
#endif

// Must initialize the static list defined in Timer.hh
list<Timer*> *Timer::TimerList = NULL;

static int ID_counter=0;
static double NtotalStartStops = 0.0;

#define TIME_UNIT_FACTOR 1.0

#ifdef PAPI
static int papiEventID=1;
static float CPU_MHZ = 0.1;

// These are just helper functions that call PAPI functions
//   Only used here and only for making code easier to read.
static int PAPI_initSingleEvent(const int event, int *tag);
static int PAPI_initDualEvent(const int event1, const int event2, int *tag);
static void PAPI_initMultiEvent(const int event1, const int event2,
                                const int event3, const int event4, int *tag);
static int PAPI_checkEvent(const int event);

static int PAPI_init(void);
static void PAPI_hwinfo(const int rank, float &mhz);
static void PAPI_initMultiplex(void);
static void PAPI_cleanup(int *i);
double computeDerivedCounter(double pc0, double pc1, double pc2);
#endif

// These are helper functions to write the formatted output
//  summary.  Nothing special.
static void writeDoubleLineSeparator(FILE *out);
static void writeSingleLineSeparator(FILE *out);
static void writeTableHeader(FILE *out);
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const int number_procs);

#ifdef PAPI
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const double counter0,
                            const double counter1,
                            const int number_procs);
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const double counter0,
                            const double counter1,
                            const double counter2,
                            const double counter3,
                            const int number_procs);
#endif
static void writeTableTotals(FILE *out,
                             const int table_number,
                             const double table_count_sum,
                             const double parent_avg,
                             const double table_avg_sum,
                             const double table_min_sum,
                             const double table_max_sum,
                             const double table_pc0_sum,
                             const double table_pc1_sum,
                             const double table_pc2_sum,
                             const double table_pc3_sum,
                             const double TimerCost);

// Timer construction for the root of a tree.
Timer::Timer(const string& name, const int tableid):
  m_name(name),  m_Parent(*this)  {

  //printf(" Root Timer: %s  ID_counter=%d tableid=%d\n", Name().c_str(), ID_counter, tableid);
  setup();
  m_tableID = tableid;

  // Make a list of pointers of all Parent Timers
  if(TimerList==NULL) TimerList = new list<Timer*>();
  TimerList->push_back(this);
}

// Non-root Managed parent/child Timer construction
Timer::Timer(const string& name, Timer& parent, const int tableid):
  m_name(name),  m_Parent(parent) {

  //printf(" Timer: %s  ID_counter=%d tableid=%d\n", Name().c_str(), ID_counter, tableid);
  setup();
  m_tableID = tableid;

  if(TimerList==NULL) TimerList = new list<Timer*>();
  TimerList->push_back(this);
}

// child-only Timer
Timer::Timer(const string& name, Timer& parent):
  m_name(name),  m_Parent(parent) {

  //printf(" child Timer: %s  ID_counter=%d\n", Name().c_str(), ID_counter);
  setup();

  if(TimerList==NULL) TimerList = new list<Timer*>();
  TimerList->push_back(this);
}

// diagnostic
Timer::Timer(const string& name, Timer& parent, const int tableid, const bool diag):
  m_name(name),  m_Parent(parent) {

  //printf(" Diagnostic Timer: %s  ID_counter=%d tableid=%d\n", Name().c_str(), ID_counter, tableid);
  setup();
  m_diagnostic = true;
  m_tableID = tableid;

  if(TimerList==NULL) TimerList = new list<Timer*>();
  TimerList->push_back(this);
}

// // Counter Timer
// Timer::Timer(const string& name, const int table):
//   timer_name(name),  Parent(parent) {
//
//   // printf(" Timer: %s  ID_counter=%d\n", Name().c_str(), ID_counter);
//   setup();
//   diagnostic_table = table;
//
//   if(TimerList==NULL) TimerList = new list<Timer*>();
//   TimerList->push_back(this);
// }

// Unmanaged Timer construction
Timer::Timer(): m_Parent(*this) {
  setup();
}

void Timer::setup() {

  //printf(" Timer setup: %s  ID: %d\n", Name().c_str(), ID_counter);

  // set this Timer to curent id
  m_ID = ID_counter;
  // increment static integer.
  ID_counter++;

  m_count = 0;
  m_tableID = -1;
  m_diagnostic = false;
  //timer_on = false;

  m_accumulated_WCtime = 0.0;
  m_last_WCtime_stamp  = 0.0;

#ifdef PAPI
  m_accumulated_counter0 = 0;
  m_accumulated_counter1 = 0;

  m_values[0] = 0;
  m_values[1] = 0;

#ifdef FOUR_COUNTERS
  m_accumulated_counter2 = 0;
  m_accumulated_counter3 = 0;

  m_values[2] = 0;
  m_values[3] = 0;
#endif
#endif
}

void Timer::TimerInit(const int rank) {

#ifdef PAPI

  // this only needs to be called once.
  // by all ranks.
  // just after MPI_init, before any start/stop calls.
  // and before any PAPI calls.  not sure where to put it.
  PAPI_init();
  PAPI_hwinfo(rank, CPU_MHZ);

#ifdef FOUR_COUNTERS
  PAPI_initMultiEvent(PAPI_TOT_CYC, PAPI_FP_INS,
                      PAPI_FMA_INS, PAPI_TLB_TL,
                      &papiEventID);
  //PAPI_FP_INS + PAPI_FMA_INS - (PM_FPU_LD_ST_ISSUES - PM_FPU_LD)

  //PAPI_initMultiEvent(PAPI_LD_INS, PAPI_L1_LDM,
  //                  PAPI_L1_STM, PAPI_L1_ICM,
  //                  &papiEventID);
  //PAPI_initMultiEvent(PAPI_TOT_CYC,
  //                  PAPI_SR_INS, PAPI_L2_LDM,
  //                PAPI_MEM_SCY,
  //                &papiEventID);

#else
  if(TIMER_COUNTER == 0) {
    PAPI_initDualEvent(PAPI_TOT_CYC, PAPI_FP_INS, &papiEventID);
    //PAPI_initDualEvent(PAPI_TOT_INS, PAPI_FP_INS, &papiEventID);
  } else if(TIMER_COUNTER == 1) {
    //PAPI_initDualEvent(PAPI_L1_TCM,  PAPI_L2_TCM, &papiEventID);
    PAPI_initDualEvent(PAPI_L1_DCM,  PAPI_L2_DCM, &papiEventID);
  } else {
    PAPI_initDualEvent(PAPI_TOT_INS, PAPI_BR_INS, &papiEventID);
  }
  //PAPI_initSingleEvent(PAPI_TOT_CYC, &myEventSet);
#endif

  PAPI_start(papiEventID);
#endif

}

// Destructor for Managed and Unmanaged Timers.
Timer::~Timer() {
  //PAPI_stop(papiEventID, m_values);
  //  cout << " Timer::~Timer() " << Name() << endl;
}

void Timer::start(void) {

  ++m_count;
  ++NtotalStartStops;
  m_last_WCtime_stamp = getTimeStampWC();

#ifdef PAPI

#ifdef NDEBUG
  PAPI_read(papiEventID, m_values);
#else
  CH_assert(PAPI_read(papiEventID, m_values) == PAPI_OK);
#endif

  m_previous_counter0 = m_values[0];
  m_previous_counter1 = m_values[1];
#ifdef FOUR_COUNTERS
  m_previous_counter2 = m_values[2];
  m_previous_counter3 = m_values[3];
#endif

#endif

  //printf(" strt Timer: %20s  v0=%20.10e pc0=%20.10e v1=%20.10e pc1=%20.10e\n",
  //     Name().c_str(),
  //     (double)m_values[0], (double)previous_counter0,
  //     (double)m_values[1], (double)previous_counter1);
  //if (m_values[0] < previous_counter0) {
  //printf(" in strt, m_values[0] < previous_counter0\n");
  //}
#ifdef MEMORY_USAGE
  sprintf(stuff, "%30s start %10d  mem=%-10.3f\n",
          Name().c_str(), (long)Count(), get_memory_usage());
  pout() << stuff;
#endif
}

void Timer::stop(void) {

  m_accumulated_WCtime += getTimeStampWC() - m_last_WCtime_stamp;

#ifdef PAPI

#ifdef NDEBUG
  PAPI_read(papiEventID, m_values);
#else
  CH_assert(PAPI_read(papiEventID, m_values) == PAPI_OK);
#endif

  m_accumulated_counter0 += m_values[0] - m_previous_counter0;
  m_accumulated_counter1 += m_values[1] - m_previous_counter1;
#ifdef FOUR_COUNTERS
  m_accumulated_counter2 += m_values[2] - m_previous_counter2;
  m_accumulated_counter3 += m_values[3] - m_previous_counter3;
#endif

#endif //PAPI

  //printf(" stop Timer: %20s  papiEventID=%2d v0=%20.10e pc0=%20.10e v1=%20.10e pc1=%20.10e\n",
  //     Name().c_str(), papiEventID,
  //     (double)m_values[0], (double)m_previous_counter0,
  //     (double)m_values[1], (double)m_previous_counter1);

  //if (m_values[0] < m_previous_counter0) {
  //printf(" in stop, m_values[0] < m_previous_counter0\n");
  //}
  //if (m_values[1] < m_previous_counter1) {
  //  printf(" in stop, m_values[1] < m_previous_counter1\n");
  //}

  // CH_assert(m_values[0] > 0);

  //#ifndef NDEBUG
  //if (m_values[0] < m_previous_counter0 || m_values[1] < m_previous_counter1) {
  //printf(" stop Timer: %20s  papiEventID=%2d v0=%20.10e pc0=%20.10e v1=%20.10e pc1=%20.10e\n",
  //       Name().c_str(), papiEventID,
  //       (double)m_values[0], (double)m_previous_counter0,
  //       (double)m_values[1], (double)m_previous_counter1);
  //}
  // CH_assert(m_values[0] >= m_previous_counter0);
  // CH_assert(m_values[1] >= m_previous_counter1);
  //#endif

#ifdef MEMORY_USAGE
  sprintf(stuff, "%30s stop  %10d  mem=%-10.3f\n",
          Name().c_str(), (long)Count(), get_memory_usage());
  pout() << stuff;
#endif

}

void Timer::clear(void) {

  //printf("clear() Timer\n");

  m_accumulated_WCtime = 0.;
  m_count = 0;

#ifdef PAPI
  m_accumulated_counter0 = 0;
  m_accumulated_counter1 = 0;
  m_previous_counter0 = 0;
  m_previous_counter1 = 0;
#ifdef FOUR_COUNTERS
  m_accumulated_counter2 = 0;
  m_accumulated_counter3 = 0;
  m_previous_counter2 = 0;
  m_previous_counter3 = 0;
#endif

#endif

}

inline double Timer::getTimeStampWC(){

#ifdef CH_MPI
  return( MPI_Wtime() );
#else
  gettimeofday(&tv, &tz);
  //printf( "%d seconds, %d microseconds\n", tv.tv_sec, tv.tv_usec);
  //cout << tvbuf->tv_sec << endl;
  //return(  ((double)((tvbuf->tv_sec)*1000000 + tvbuf->tv_usec))*1.0e-6 );

  // this doesn't always work cuz of potential integer overflow.
  //return( (tv.tv_sec * 1000000 + tv.tv_usec) * 1.0e-6 );

  // this is what MPICH uses for MPI_Wtime on machines
  // which have gettimeofday.
  return((double)tv.tv_sec + 0.000001 * (double)tv.tv_usec);
#endif
}

// This is the Summary for the Managed Timers.  If running in
// parallel, i first reduce all of the accumulated times for each
// timer and get the minimum, maximum, and average.  Then i go thru
// the list of timers and make a list of parent timers.  From there i
// can make for loops that step thru the tree and print out the
// results.
void Timer::TimerSummary(void) {

  // all ranks go in, but only rank0 writes the file.
  TimerSummary_();

  // all ranks need to clean up
#ifdef PAPI

  // not sure how important this cleanup and shutdown is...
  //for(int i=0; i<ID_counter; i++) {
  //printf(" cleaning up event %d\n", i);
  //PAPI_cleanup(&i);
  //}

  PAPI_shutdown();
#endif

  // clean up
  delete TimerList;

}

void Timer::TimerSummary_(void) {

  int rank=0, number_procs=1;
#ifdef CH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &number_procs);
#endif

  //printf(" rank%3d in Timer::TimerSummary_\n", rank);

  if (rank==0) {
#ifdef PAPI
    printf("defined: PAPI\n");
#endif
#ifdef CH_AIX
    printf("defined: CH_AIX\n");
#endif
#ifdef FOUR_COUNTERS
    printf("defined: FOUR_COUNTERS\n");
#endif
#ifdef WRITE_TITA_FILES
    printf("defined: WRITE_TITA_FILES\n");
#endif
#ifdef CH_MPI
    printf("defined: MPI\n");
#endif
#ifdef MEMORY_USAGE
    printf("defined: MEMORY_USAGE\n");
#endif
#ifdef CATFISH
    printf("defined: CATFISH\n");
#endif
#ifdef NDEBUG
    printf("defined: NDEBUG\n");
#endif
  }

  // don't write these files out all the darn time...
#ifdef WRITE_TITA_FILES
  FILE *TT;
  char ttfilename[80];
  if (number_procs < 999) {
    sprintf(ttfilename, "tita%03d", rank);
  } else {
    sprintf(ttfilename, "tita%04d", rank);
  }
  TT = fopen(ttfilename, "w");
  if(TT == NULL) {
    printf("problem opening tita\n");
    exit(0);
  }

#ifdef PAPI
  for(list<Timer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; ++tli) {
    double pc0 = (double)(*tli)->papi_counter0();
    double pc1 = (double)(*tli)->papi_counter1();
#ifdef FOUR_COUNTERS
    double pc2 = (double)(*tli)->papi_counter2();
    double pc3 = (double)(*tli)->papi_counter3();
    double dc =  computeDerivedCounter(pc0, pc1, pc2);
    fprintf(TT, " %-16s P: %-16s d%-3d n=%14.0f wc=%13.3f pc0=%17.0f pc1=%17.0f pc2=%17.0f pc3=%17.0f dc=%13.3f\n",
            (*tli)->Name().c_str(), (*tli)->m_Parent.Name().c_str(),
            (*tli)->diagnostic_table, (double)((*tli)->Count()), (*tli)->wc_time(),
            pc0, pc1, pc2, pc3, dc);
#else
    double dc =  computeDerivedCounter(pc0, pc1, 0);
    fprintf(TT, " %-16s P: %-16s d%-3d n=%14.0f wc=%13.3f pc0=%17.0f pc1=%17.0f dc=%13.3f\n",
            (*tli)->Name().c_str(), (*tli)->m_Parent.Name().c_str(),
            (*tli)->diagnostic_table, (double)((*tli)->Count()), (*tli)->wc_time(),
            pc0, pc1, dc);
#endif // FOUR_COUNTERS
  }
#else
  for(list<Timer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; ++tli) {
    fprintf(TT, " %-16s P: %-16s d%-3d n=%14.0f  wc=%13.3f\n",
            (*tli)->Name().c_str(), (*tli)->m_Parent.Name().c_str(),
            (*tli)->diagnostic_table, (double)((*tli)->Count()), (*tli)->wc_time());
  }
#endif // PAPI
  fclose(TT);
#endif // WRITE_TITA_FILES

  double largest_time = 0.0;

  // go thru all Timers in static list and obtain the avg/min/max of
  // all processors onto rank0 -- obviously trivial for serial runs.
  for(list<Timer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; ++tli) {

    //printf(" Timer: %-16s Parent: %-16s diagnostic = %d  n=%d  wc=%10.3f\n",
    // (*tli)->Name().c_str(), (*tli)->m_Parent.Name().c_str(),
    // (*tli)->diagnostic_table), (*tli)->Count(), (*tli)->wc_time();

    double wc = (*tli)->wc_time();
    long long int nc = (*tli)->Count();

    if (number_procs == 1) {
      (*tli)->m_avgWC = wc;
      (*tli)->m_minWC = wc;
      (*tli)->m_maxWC = wc;
      (*tli)->m_totalCount = nc;
      (*tli)->m_avgCount = (double)nc;
    } else {

#ifdef CH_MPI

      // after here, rank0 will be the only guy with the right answers
      MPI_Reduce(&wc, &(*tli)->m_minWC, 1, MPI_DOUBLE,
                 MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&wc, &(*tli)->m_maxWC, 1, MPI_DOUBLE,
                 MPI_MAX, 0, MPI_COMM_WORLD);

      double temp;
      MPI_Reduce(&wc, &temp, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);

      if(rank==0) (*tli)->m_avgWC = temp/number_procs;

      MPI_Reduce(&nc, &(*tli)->m_totalCount, 1, MPI_LONG_LONG,
                 MPI_SUM, 0, MPI_COMM_WORLD);

      if(rank==0) {
        (*tli)->m_avgCount = (double)(*tli)->m_totalCount/number_procs;
      }
#endif

    }

    if ((*tli)->m_avgWC > largest_time) largest_time = (*tli)->m_avgWC;
    // serial or parallel.  assumes serial has rank set to 0.
    //if(rank==0) {
      // in parallel, totalCount is the sum from all processors
      // but only exists on rank0.
      //NtotalStartStops += (double)(*tli)->totalCount;
    //}
  }

#ifdef PAPI
  // go thru all Timers in static list and obtain the sum of
  // PAPI counters onto rank0 -- obviously trivial for serial runs.
  for(list<Timer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; ++tli) {

    //printf(" Timer: %-16s Parent: %-16s diagnostic = %d  pc0=%14.5e pc1=%14.5e\n",
    ///     (*tli)->Name().c_str(), (*tli)->m_Parent.Name().c_str(),
    //    (*tli)->diagnostic_table,
    //     (double)(*tli)->papi_counter0(),
    //     (double)(*tli)->papi_counter1());

    double pc0 = (double)(*tli)->papi_counter0();
    double pc1 = (double)(*tli)->papi_counter1();
#ifdef FOUR_COUNTERS
    double pc2 = (double)(*tli)->papi_counter2();
    double pc3 = (double)(*tli)->papi_counter3();
#endif

    if (number_procs == 1) {
      (*tli)->m_totalPapiCounter0 = pc0;
      (*tli)->m_totalPapiCounter1 = pc1;
#ifdef FOUR_COUNTERS
      (*tli)->m_totalPapiCounter2 = pc2;
      (*tli)->m_totalPapiCounter3 = pc3;
#endif

    } else {

#ifdef CH_MPI
      MPI_Reduce(&pc0, &(*tli)->m_totalPapiCounter0, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&pc1, &(*tli)->m_totalPapiCounter1, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);
#ifdef FOUR_COUNTERS
      MPI_Reduce(&pc2, &(*tli)->m_totalPapiCounter2, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&pc3, &(*tli)->m_totalPapiCounter3, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);
#endif
#endif
    }

    //if (rank==0) {
    // printf(" Timer: %-16s after reduce pc0=%14.5e pc1=%14.5e  pc0=%14.5e pc1=%14.5e\n",
    //        (*tli)->Name().c_str(),
    //        (*tli)->m_totalPapiCounter0,
    //        (*tli)->m_totalPapiCounter1,
    //        (*tli)->total_papi_counter0(),
    //        (*tli)->total_papi_counter1());
    // }
  }
#endif  // PAPI

#ifdef CH_MPI
  double combinedTotalStartStops=0;
  MPI_Reduce(&NtotalStartStops, &combinedTotalStartStops, 1, MPI_DOUBLE,
             MPI_SUM, 0, MPI_COMM_WORLD);
  double avgStartStops = combinedTotalStartStops/number_procs;
#else
  double avgStartStops = NtotalStartStops;
#endif

  // go away if you aren't rank 0
  if (rank != 0) { return; }

  //double largest_time = 0.0;
  //for(list<Timer*>::iterator pti=ParentList.begin() ; pti != ParentList.end() ; pti++) {
  //double pa = (*pti)->m_avgWC*TIME_UNIT_FACTOR;
  // if (pa > largest_time) largest_time = pa;
  //}

  Timer TimerLoop;
  Timer iLooper0;
  const int Nloop = 10000;
  TimerLoop.start();
  for(int i=0; i<Nloop; i++) {
    iLooper0.start();    iLooper0.stop();
  }
  TimerLoop.stop();
  double TimerCost = TimerLoop.wc_time()/(double)iLooper0.Count();
  if (largest_time <= 0.0) largest_time = 0.01;

#ifdef CATFISH
  cout << "Single Timer Start/Stop Cost = " << TimerCost << " sec." << endl;
  cout << "Total Timer Start/Stop Calls on rank0= " << NtotalStartStops << endl;
  cout << "Avg Timer Start/Stop Calls over all procs= " << avgStartStops << endl;
  cout << "Avg Estimated Timer Cost per Proc= " << avgStartStops*TimerCost
         << " sec. (" << (avgStartStops*TimerCost)/largest_time*100.0
         << "%) " << endl;
#else
  pout() << "Single Timer Start/Stop Cost = " << TimerCost << " sec." << endl;
  pout() << "Total Timer Start/Stop Calls on rank0= " << NtotalStartStops << endl;
  pout() << "Avg Timer Start/Stop Calls over all procs= " << avgStartStops << endl;
  pout() << "Avg Estimated Timer Cost per Proc= " << avgStartStops*TimerCost
         << " sec. (" << (avgStartStops*TimerCost)/largest_time*100.0
         << "%) " << endl;
#endif

  cout << " rank" << rank << " writing time.table " << endl;

  FILE *OUT;

  if(TIMER_COUNTER == 0) {
    OUT = fopen("time.table", "w");
  } else if (TIMER_COUNTER == 1) {
    OUT = fopen("time.table1", "w");
  } else {
    OUT = fopen("time.table2", "w");
  }
  if(OUT == NULL) {
    printf("problem opening output file in Timer\n");
    exit(0);
  }

  fprintf(OUT, "Number of Nodes: %d\n", number_procs );
  fprintf(OUT, "\n");

  writeParentTables(OUT, TimerCost);
  writeDiagnosticTables(OUT, TimerCost);

  fprintf(OUT, "Single Timer Start/Stop Cost = %15.7e sec.\n", TimerCost);
  fprintf(OUT, "Total Timer Start/Stop Calls = %20.10e \n", NtotalStartStops);
  fprintf(OUT, "Avg Timer Start/Stop per Proc= %20.10e \n", NtotalStartStops/number_procs);
  fprintf(OUT, "Avg Est Timer Cost per Proc  = %15.3f sec.  (%7.2f%%)\n",
          NtotalStartStops*TimerCost/number_procs,
          (NtotalStartStops*TimerCost/number_procs)/largest_time*100.0);

  fclose(OUT);

  //FILE *QOUT;
  //char qfilename[80];
  //sprintf(qfilename, "q%03d.dat", rank);
  //QOUT = fopen(qfilename, "w");
  //if(QOUT == NULL) {
  //printf("problem opening pex out\n");
  //exit(0);
  //}

  //  for(unsigned int i=0; i<NQ; i++) {
  //fprintf(OUT,"%20.10e  s0=%20.10e s1=%20.10e s2=%20.10e\n",
  //        (double)i, (double)qsubbox_len0[i],
  //        (double)qsubbox_len1[i], (double)qsubbox_len2[i]);
  //}
  //fprintf(OUT,"iter= %20.10e N=%20.10e\n",
  //      (double)qiter, (double)qNN);

  //for(unsigned int i=0; i<qiter; i++) {
  //fprintf(QOUT, "%12d %4d %4d %4d %4d\n",
  //        i, qNN[i], qsubbox_len2[i], qsubbox_len1[i], qsubbox_len0[i]);
  //}
  //fclose(QOUT);

}

void Timer::writeParentTables(FILE *out, const double TimerCost) {

  int number_procs=1;
#ifdef CH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &number_procs);
#endif

  // Create a list of Parent Timers from the list of Timers
  list<Timer*> ParentList;

  //cout << " size of TimerList = " << TimerList.size() << endl;

  double largest_time=0.0;

  // for each timer in entire Timer List
  for(list<Timer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; tli++) {

    if ((*tli)->m_avgWC > largest_time) largest_time = (*tli)->m_avgWC;

    //printf(" Timer #%3d: %-16s Parent: %-16s tableid = %d\n",
    //     (*tli)->m_ID,
    //     (*tli)->Name().c_str(),
    //     (*tli)->m_Parent.Name().c_str(),
    //     (*tli)->tableID());

    // add the Parent of this Timer to the Parent List
    bool add = false;
    if ((*tli)->m_Parent.tableID() >= 0) add = true;
    if ((*tli)->m_diagnostic) add = false;

    // but first, make sure it isn't already in the Parent List
    // for each parent currently in Parent List
    for(list<Timer*>::iterator pti=ParentList.begin() ; pti != ParentList.end() ; pti++) {
      if(*pti == &((*tli)->m_Parent)) {
        add = false;
        break;
      }
    }

    // add them in ascending order based on the table ID
    if(add) {
      //cout << "  considering Parent:" << (*tli)->m_Parent.Name()
      //   << " tid=" << (*tli)->m_Parent.tableID() << endl;

      if (ParentList.size() == 0) {
        ParentList.push_back( &((*tli)->m_Parent) );
        //cout << " first one" << endl;
      } else {
        bool inserted=false;
        for(list<Timer*>::iterator pti=ParentList.begin() ; pti != ParentList.end() ; pti++) {
          if ((*tli)->m_Parent.tableID() < (*pti)->tableID()) {
            ParentList.insert(pti, &((*tli)->m_Parent) );
            inserted = true;
            //cout << " inserted" << endl;
            break;
          }
        }
        if (!inserted) {
          ParentList.push_back( &((*tli)->m_Parent) );
          //cout << " push_back" << endl;
        }
      }
      //cout << "  Adding Parent:" << (*tli)->m_Parent.Name()
      //   << " tableid = " << (*tli)->tableID()
      //   << " wc=" << (*tli)->m_Parent.m_avgWC << endl;
    }
    //cout << " ParentList size = " << ParentList.size() << endl;
  }

  //for(list<Timer*>::iterator pti=ParentList.begin() ; pti != ParentList.end() ; pti++) {
  //printf(" Parent #%3d: %-16s Parent: %-16s tableid = %d\n",
  //       (*pti)->m_ID,
  //       (*pti)->Name().c_str(),
  //       (*pti)->m_Parent.Name().c_str(),
  //       (*pti)->tableID());
  //}

  // for every Parent Timer -- make a new table
  for(list<Timer*>::iterator pti=ParentList.begin() ; pti != ParentList.end() ; pti++) {

    double table_avg_sum     = 0.0;
    double table_min_sum     = 0.0;
    double table_max_sum     = 0.0;
    double table_count_sum   = 0.0;
    double table_pc0_sum     = 0.0;
    double table_pc1_sum     = 0.0;
    double table_pc2_sum     = 0.0;
    double table_pc3_sum     = 0.0;
    double table_percent;

    double wc_avg = (*pti)->m_avgWC*TIME_UNIT_FACTOR;
    double wc_min = (*pti)->m_minWC*TIME_UNIT_FACTOR;
    double wc_max = (*pti)->m_maxWC*TIME_UNIT_FACTOR;
    double count_avg = (*pti)->m_avgCount;

    double parent_avg = wc_avg;

    if (count_avg == 0 && wc_avg == 0.0 && wc_min == 0.0) {
      // prolly don't want this table printed cuz it will prolly be all zeros
      continue;
    }

    fprintf(out, "\n");

    writeTableHeader(out);

    // Parent for each table
#ifdef PAPI
#ifdef FOUR_COUNTERS
    writeLineOfData(out, (*pti)->Name().c_str(), count_avg, 100*wc_avg/largest_time,
                    wc_avg, wc_min, wc_max,
                    (*pti)->total_papi_counter0(),
                    (*pti)->total_papi_counter1(),
                    (*pti)->total_papi_counter2(),
                    (*pti)->total_papi_counter3(),
                    number_procs);
#else
    writeLineOfData(out, (*pti)->Name().c_str(), count_avg, 100*wc_avg/largest_time,
                    wc_avg, wc_min, wc_max,
                    (*pti)->total_papi_counter0(),
                    (*pti)->total_papi_counter1(),
                    number_procs);
#endif

#else
    writeLineOfData(out, (*pti)->Name().c_str(), count_avg, 100*wc_avg/largest_time,
                    wc_avg, wc_min, wc_max, number_procs);
#endif

    writeSingleLineSeparator(out);

    // for every Timer (looking for children)
    for(list<Timer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; tli++) {

      if ((*tli)->m_diagnostic) continue;

      // if this Timer's Parent is the current Parent and
      //  if this Timer's Parent isn't equal to itself (Everything)
      if( &((*tli)->m_Parent) == *pti && *tli != &((*tli)->m_Parent) ) {

        wc_avg = (*tli)->m_avgWC*TIME_UNIT_FACTOR;
        wc_min = (*tli)->m_minWC*TIME_UNIT_FACTOR;
        wc_max = (*tli)->m_maxWC*TIME_UNIT_FACTOR;
        count_avg = (*tli)->m_avgCount;

        if(parent_avg  > 0.) {
          table_percent  = wc_avg/parent_avg*100.0;
        } else {
          table_percent = 0.;
        }

        table_avg_sum  += wc_avg;
        table_min_sum  += wc_min;
        table_max_sum  += wc_max;
        table_count_sum += count_avg;

        // Children
#ifdef PAPI
        table_pc0_sum += (*tli)->total_papi_counter0();
        table_pc1_sum += (*tli)->total_papi_counter1();
#ifdef FOUR_COUNTERS
        table_pc2_sum += (*tli)->total_papi_counter2();
        table_pc3_sum += (*tli)->total_papi_counter3();
        writeLineOfData(out, (*tli)->Name().c_str(), count_avg,
                        table_percent, wc_avg, wc_min, wc_max,
                        (*tli)->total_papi_counter0(),
                        (*tli)->total_papi_counter1(),
                        (*tli)->total_papi_counter2(),
                        (*tli)->total_papi_counter3(),
                        number_procs);
#else
        writeLineOfData(out, (*tli)->Name().c_str(), count_avg,
                        table_percent, wc_avg, wc_min, wc_max,
                        (*tli)->total_papi_counter0(),
                        (*tli)->total_papi_counter1(),
                        number_procs);
#endif
#else
        writeLineOfData(out, (*tli)->Name().c_str(), count_avg,
                        table_percent, wc_avg, wc_min, wc_max, number_procs);
#endif
      }

    } // end of table loop for this parent

    writeTableTotals(out, (*pti)->tableID(),
                     table_count_sum, parent_avg,
                     table_avg_sum, table_min_sum,  table_max_sum,
                     table_pc0_sum, table_pc1_sum,
                     table_pc2_sum, table_pc3_sum, TimerCost);

    fprintf(out,"\n");
  } // end of parent loop
}

void Timer::writeDiagnosticTables(FILE *out, const double TimerCost) {

  // Create a set of Diagnostic Tables -- ie all of the Timers
  // with unique values of member data "diagnostic_table".
  set<int> diagTableSet;

  for(list<Timer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; tli++) {
    //printf(" Timer #%3d: %-16s Parent: %-16s tableid = %d diag=%d\n",
    //     (*tli)->m_ID,
    //     (*tli)->Name().c_str(),
    //     (*tli)->m_Parent.Name().c_str(),
    //     (*tli)->tableID(), (*tli)->m_diagnostic);

    if ((*tli)->m_diagnostic) {
      diagTableSet.insert((*tli)->tableID());
    }
  }

  //cout << " Number of Diag Tables = " << diagTableSet.size() << endl;
  //for (set<int>::iterator si=diagTableSet.begin(); si != diagTableSet.end(); si++) {
  // cout << " si = " << *si << endl;
  //}

  int number_procs=1;
#ifdef CH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &number_procs);
#endif

  // For each unique value of diagnostic_table, print a table.
  for (set<int>::iterator si=diagTableSet.begin(); si != diagTableSet.end(); si++) {

    fprintf(out, "\n\n Diagnostic Table %d\n", *si);

    // crude way of getting the Parent for this diagnostic table.
    // part of the crudeness is that i'm not checking that all
    // Timers in this diagnostic table have the same parent.
    Timer *DiagnosticParent = NULL;
    for(list<Timer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; tli++) {
      if ((*tli)->tableID() != *si) continue;
      DiagnosticParent = &((*tli)->m_Parent);
    }

    writeTableHeader(out);

    double table_avg_sum   = 0.0;
    double table_min_sum   = 0.0;
    double table_max_sum   = 0.0;
    double table_count_sum = 0.0;
    double table_pc0_sum   = 0.0;
    double table_pc1_sum   = 0.0;
    double table_pc2_sum   = 0.0;
    double table_pc3_sum   = 0.0;
    double table_percent;

    double wc_avg = DiagnosticParent->m_avgWC*TIME_UNIT_FACTOR;
    double wc_min = DiagnosticParent->m_minWC*TIME_UNIT_FACTOR;
    double wc_max = DiagnosticParent->m_maxWC*TIME_UNIT_FACTOR;
    double count_avg = DiagnosticParent->m_avgCount;

    double parent_avg = wc_avg;

    // Parent for each table
#ifdef PAPI
#ifdef FOUR_COUNTERS
    writeLineOfData(out, DiagnosticParent->Name().c_str(),
                    count_avg, -1,
                    wc_avg, wc_min, wc_max,
                    DiagnosticParent->total_papi_counter0(),
                    DiagnosticParent->total_papi_counter1(),
                    DiagnosticParent->total_papi_counter2(),
                    DiagnosticParent->total_papi_counter3(),
                    number_procs);
#else
    writeLineOfData(out, DiagnosticParent->Name().c_str(),
                    count_avg, -1,
                    wc_avg, wc_min, wc_max,
                    DiagnosticParent->total_papi_counter0(),
                    DiagnosticParent->total_papi_counter1(),
                    number_procs);
#endif
#else
    writeLineOfData(out, DiagnosticParent->Name().c_str(),
                    count_avg, -1,
                    wc_avg, wc_min, wc_max, number_procs);
#endif

    writeSingleLineSeparator(out);

    for(list<Timer*>::iterator tli=TimerList->begin() ; tli != TimerList->end() ; tli++) {

      // skip all Timers except ones belonging in this table
      if ((*tli)->tableID() != *si) continue;

      wc_avg = (*tli)->m_avgWC*TIME_UNIT_FACTOR;
      wc_min = (*tli)->m_minWC*TIME_UNIT_FACTOR;
      wc_max = (*tli)->m_maxWC*TIME_UNIT_FACTOR;
      count_avg = (*tli)->m_avgCount;

      if(parent_avg  > 0.) {
        table_percent  = wc_avg/parent_avg*100.0;
      } else {
        table_percent = 0.;
      }

      table_avg_sum  += wc_avg;
      table_min_sum  += wc_min;
      table_max_sum  += wc_max;
      table_count_sum += count_avg;

      // Children
#ifdef PAPI
      table_pc0_sum += (*tli)->total_papi_counter0();
      table_pc1_sum += (*tli)->total_papi_counter1();

#ifdef FOUR_COUNTERS
      table_pc2_sum += (*tli)->total_papi_counter2();
      table_pc3_sum += (*tli)->total_papi_counter3();
      writeLineOfData(out, (*tli)->Name().c_str(), count_avg,
                      table_percent, wc_avg, wc_min, wc_max,
                      (*tli)->total_papi_counter0(),
                      (*tli)->total_papi_counter1(),
                      (*tli)->total_papi_counter2(),
                      (*tli)->total_papi_counter3(),
                      number_procs);
#else
      writeLineOfData(out, (*tli)->Name().c_str(), count_avg,
                      table_percent, wc_avg, wc_min, wc_max,
                      (*tli)->total_papi_counter0(),
                      (*tli)->total_papi_counter1(),
                      number_procs);
#endif

#else
      writeLineOfData(out, (*tli)->Name().c_str(), count_avg,
                      table_percent, wc_avg, wc_min, wc_max, number_procs);
#endif

    } // end of table loop for this parent

    writeTableTotals(out, -1*(*si),
                     table_count_sum, parent_avg,
                     table_avg_sum, table_min_sum,  table_max_sum,
                     table_pc0_sum, table_pc1_sum,
                     table_pc2_sum, table_pc3_sum, TimerCost);

    fprintf(out,"\n\n");
  }
}

static void writeSingleLineSeparator(FILE *out) {
#ifdef PAPI
  fprintf(out, "-------------------------------------------------");
  fprintf(out, "-------------------------------------------------\n");
#else
  fprintf(out, "-------------------------------------------------");
  fprintf(out, "---------------------------\n");
#endif
}

static void writeDoubleLineSeparator(FILE *out) {
#ifdef PAPI
  fprintf(out, "=================================================");
  fprintf(out, "=================================================\n");
#else
  fprintf(out, "=================================================");
  fprintf(out, "===========================\n");
#endif
}

static void writeTableHeader(FILE *out) {
#ifdef PAPI
  //            12345678901234567890
  fprintf(out, "      Totals    ");
  fprintf(out, "  WC per.");
  fprintf(out, "      count");

  fprintf(out, "    avg  ");
  fprintf(out, "     min  ");
  fprintf(out, "     max  ");

#ifdef FOUR_COUNTERS
  fprintf(out, "    ");
#else
  if(TIMER_COUNTER==0) {
    fprintf(out, "    cycles     FP INS     mflops");
  } else if(TIMER_COUNTER==1) {
    fprintf(out, "    L1 TCM     L2 TCM     L2/L1");
  } else {
    fprintf(out, "    TOT INS    BR INS     BR/TOT");
  }
#endif

#else
  //            12345678901234567890
  fprintf(out, "      Totals    ");
  fprintf(out, "  WC per.");
  fprintf(out, "      count");

  fprintf(out, "        avg  ");
  fprintf(out, "       min  ");
  fprintf(out, "       max  ");
#endif

  fprintf(out, "\n");

  writeDoubleLineSeparator(out);

}

static void writeTableTotals(FILE *out,
                             const int table_number,
                             const double count_sum,
                             const double parent_avg,
                             const double avg_sum,
                             const double min_sum,
                             const double max_sum,
                             const double pc0_sum,
                             const double pc1_sum,
                             const double pc2_sum,
                             const double pc3_sum,
                             const double TimerCost) {

//   const char tcs[16] = "Timer Cost";

//   double totalTimerCost = (double)count_sum * TimerCost;
//   double percentTimerCost;
//   if(parent_avg  > 0.) {
//     percentTimerCost = totalTimerCost/parent_avg * 100.0;
//   } else {
//     percentTimerCost = 0.;
//   }

// #ifdef PAPI
//   fprintf(out, "%16s (%6.2f%%) %9ld %8.2f [%8.2f,%8.2f]  %10.3e %10.3e %8.2f\n",
//           tcs, percentTimerCost, (long)count_sum, totalTimerCost, 0.0, 0.0, 0.0, 0.0, 0.0);
// #else
//   fprintf(out, "%16s (%6.2f%%) %9ld %10.2f [%10.2f, %10.2f]\n",
//           tcs, percentTimerCost, (long)count_sum, totalTimerCost, 0.0, 0.0);
// #endif

  writeDoubleLineSeparator(out);

  double percent;
  if(parent_avg  > 0.) {
    percent = avg_sum/parent_avg * 100.0;
  } else {
    percent = 0.;
  }

  char stuff2[30];
  if (table_number >= 0) {
    sprintf(stuff2, "T%-4d table tots: ", table_number);
  } else {
    sprintf(stuff2, "DT%-4dtable tots: ", -1*table_number);
  }
#ifdef PAPI

#ifdef FOUR_COUNTERS
  double derived_counter = computeDerivedCounter(pc0_sum, pc1_sum, pc2_sum);
  fprintf(out, "%s%6.2f%%  %9ld %8.2f [%8.2f,%8.2f]  %10.3e %10.3e %10.3e %10.3e %8.2f\n",
          stuff2, percent, (long)count_sum,
          avg_sum, min_sum, max_sum,
          pc0_sum, pc1_sum,  pc2_sum, pc3_sum, derived_counter);
#else
  double derived_counter = computeDerivedCounter(pc0_sum, pc1_sum, 0);
  fprintf(out, "%s%6.2f%%  %9ld %8.2f [%8.2f,%8.2f]  %10.3e %10.3e %8.2f\n",
          stuff2, percent, (long)count_sum,
          avg_sum, min_sum, max_sum,
          pc0_sum, pc1_sum, derived_counter);
#endif
#else
  fprintf(out, "%s%6.2f%%  %9ld %10.2f [%10.2f, %10.2f]\n",
          stuff2, percent, (long)count_sum,
          avg_sum, min_sum, max_sum);
#endif

}

#ifdef PAPI
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const double counter0,
                            const double counter1,
                            const int number_procs) {
  // PAPI counters should be coming in as sum totals of all
  //  processors.
  double derived_counter = computeDerivedCounter(counter0, counter1, 0);

  if(percent < 0) {
    fprintf(out, "%16s           %9.0f %8.2f [%8.2f,%8.2f]  %10.3e %10.3e %8.2f\n",
            name, count, wc_avg, wc_min, wc_max,
            counter0, counter1,
            derived_counter);
  } else {
    fprintf(out, "%16s  %6.2f%%  %9.0f %8.2f [%8.2f,%8.2f]  %10.3e %10.3e %8.2f\n",
            name, percent, count, wc_avg, wc_min, wc_max,
            counter0, counter1,
            derived_counter);
  }

}
#endif

#if defined(PAPI) && defined(FOUR_COUNTERS)
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const double counter0,
                            const double counter1,
                            const double counter2,
                            const double counter3,
                            const int number_procs) {
  // PAPI counters should be coming in as sum totals of all
  //  processors.
  double derived_counter = computeDerivedCounter(counter0, counter1, counter2);

  if(percent < 0) {
    fprintf(out, "%16s           %9.0f %8.2f [%8.2f,%8.2f]  %10.3e %10.3e %10.3e %10.3e %8.2f\n",
            name, count, wc_avg, wc_min, wc_max,
            counter0, counter1, counter2, counter3,
            derived_counter);
  } else {
    fprintf(out, "%16s  %6.2f%%  %9.0f %8.2f [%8.2f,%8.2f]  %10.3e %10.3e %10.3e %10.3e %8.2f\n",
            name, percent, count, wc_avg, wc_min, wc_max,
            counter0, counter1, counter2, counter3,
            derived_counter);
  }
}
#endif

// regular, no PAPI
static void writeLineOfData(FILE *out,
                            const char *name,
                            const double count,
                            const double percent,
                            const double wc_avg,
                            const double wc_min,
                            const double wc_max,
                            const int number_procs) {

  if(percent < 0) {
    fprintf(out, "%16s           %9.0f %10.2f [%10.2f, %10.2f]\n",
            name, count, wc_avg, wc_min, wc_max);
  } else {
    fprintf(out, "%16s  %6.2f%%  %9.0f %10.2f [%10.2f, %10.2f]\n",
            name, percent, count, wc_avg, wc_min, wc_max);
  }

}

#ifdef PAPI

double computeDerivedCounter(double pc0, double pc1, double pc2) {

  const double VALID_PROC_TIME = 1.0e-4;

  double proc_time = 1.0;
  double derived_counter = 0.0;

  if(TIMER_COUNTER == 0) {
    if(CPU_MHZ <= 0.0) CPU_MHZ = 1.0;
    proc_time = pc0/(CPU_MHZ*1.0e6);
    if(proc_time > VALID_PROC_TIME) {
#ifdef FOUR_COUNTERS
      // FP_INS + FP_FMA
      derived_counter = (pc1+pc2)/(proc_time*1.0e6);
#else
      derived_counter = (pc1)/(proc_time*1.0e6);
#endif
    } else {
      derived_counter = 0.0;
    }
  } else if (TIMER_COUNTER == 1) {
    if(pc0 > 0) {
      derived_counter = pc1/pc0;
    } else {
      derived_counter = 0.0;
    }
  } else {
    if(pc0 > 0) {
      derived_counter = pc1/pc0;
    } else {
      derived_counter = 0.0;
    }
  }

  return derived_counter;
}

static int PAPI_initSingleEvent(const int event, int *tag) {

  int retval;
  retval = PAPI_query_event(event);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_create_eventset(tag);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_add_event(tag, event);
  CH_assert(retval == PAPI_OK);
  return retval; // to get rid of warning mesg
}

static int PAPI_initDualEvent(const int event1, const int event2,
                               int *tag) {
  int retval;

  retval = PAPI_query_event(event1);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_query_event(event2);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_create_eventset(tag);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_add_event(tag, event1);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_add_event(tag, event2);
  CH_assert(retval == PAPI_OK);
  return retval;  // to get rid of warning mesg
}

static void PAPI_initMultiEvent(const int event1, const int event2,
                                const int event3, const int event4,
                                int *tag) {

  int retval;

  retval = PAPI_query_event(event1);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_query_event(event2);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_query_event(event3);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_query_event(event4);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_create_eventset(tag);
  CH_assert(retval == PAPI_OK);

  retval = PAPI_add_event(tag, event1);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_add_event(tag, event2);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_add_event(tag, event3);
  CH_assert(retval == PAPI_OK);
  retval = PAPI_add_event(tag, event4);
  CH_assert(retval == PAPI_OK);

  //retval = PAPI_set_multiplex(tag);
  // CH_assert(retval == PAPI_OK);
}

static int PAPI_checkEvent(const int event) {
  int tag = PAPI_NULL;
  long long int value = -1;
  int retval;

  if( (retval = PAPI_query_event(event)) != PAPI_OK ) {
    printf(" problem with PAPI_query_event retval=%d\n", retval);
    return retval;
  }
  if( (PAPI_create_eventset(&tag) != PAPI_OK) ) {
    printf(" problem with PAPI_create_eventset\n");
    return retval;
  }

  if( (PAPI_add_event(&tag, event) != PAPI_OK) ) {
    printf(" problem with PAPI_add_event\n");
    return retval;
  }
  if( (PAPI_start(tag) != PAPI_OK) ) {
    printf(" problem with PAPI_start\n");
    return retval;
  }
  if( (PAPI_stop(tag, &value) != PAPI_OK) ) {
    printf(" problem with PAPI_stop\n");
    return retval;
  }
  printf(" event looks ok.  value=%f\n", (double)value);
  return retval; // to get rid of warning mesg
}

static int PAPI_init(void) {
  int retval;
  retval = PAPI_library_init(PAPI_VER_CURRENT);
  CH_assert(retval == PAPI_VER_CURRENT);

  //retval = PAPI_start_counters(&Events,1))

  return retval; // to get rid of warning mesg
}

static void PAPI_hwinfo(const int rank, float &mhz) {

  const PAPI_hw_info_t *hwinfo = NULL;
  hwinfo = PAPI_get_hardware_info();
  CH_assert( hwinfo != NULL );

  if (rank==0) {
    printf("Vendor string and code   : %s (%d)\n",hwinfo->vendor_string,hwinfo->vendor);
    printf("Model string and code    : %s (%d)\n",hwinfo->model_string,hwinfo->model);
    printf("CPU revision             : %f\n",hwinfo->revision);
    printf("CPU Megahertz            : %f\n",hwinfo->mhz);
    printf("CPU's in an SMP node     : %d\n",hwinfo->ncpu);
    printf("Nodes in the system      : %d\n",hwinfo->nnodes);
    printf("Total CPU's in the system: %d\n",hwinfo->totalcpus);

    int num_counters=PAPI_num_counters();
    printf("Number PAPI counters     : %d\n", num_counters);
    CH_assert (num_counters >= PAPI_OK);
  }

  mhz = hwinfo->mhz;
}

static void PAPI_initMultiplex(void) {
  int retval = PAPI_multiplex_init();
  CH_assert(retval == PAPI_OK);
}

static void PAPI_cleanup(int *i) {
#ifdef NDEBUG  // no debug
  PAPI_cleanup_eventset(i);
#else
  CH_assert(PAPI_cleanup_eventset(i) == PAPI_OK);
#endif
}

#endif // on PAPI

#else
#include "Timer.H"

#endif  // TIMER
