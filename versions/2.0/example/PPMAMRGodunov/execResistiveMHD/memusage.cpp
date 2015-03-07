#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include <cstdio>
// note:  i copied the following #include crud from elsewhere.
// assumed it was still relevant (ndk)
#include <sys/time.h>

#include "memusage.H"
#include "MayDay.H"
#include "ParmParse.H"
#include "CH_HDF5.H"

#ifdef __SOLARIS2__
/* This is a problem.  Version 2.5 of Solaris (at least) does not have
 * rusage.h, preferring to include the definition of struct rusage in
 * resource.h.  Version 2.3 of Solaris (at least) has resource.h but
 * it does not include struct rusage, which is rather defined in rusage.h
 *
 * The _right_ fix
 */
#  ifdef USE_RUSAGE_H
#    include <sys/rusage.h>
#  else
#    include <sys/resource.h>
#  endif
#  include <sys/procfs.h>
#else
#  include <sys/resource.h>
#endif

#ifdef __hpux
#include <sys/syscall.h>
#define getrusage(a,b) syscall(SYS_getrusage,a,b)
#endif

// returns memory used in MB
Real get_memory_usage(void) {

  Real memusage = 0.0;

  // if linux, use a different method.
#ifdef CH_Linux
  FILE* f=fopen("/proc/self/statm","r");
  int size, resident, shared;
  if(fscanf(f, "%d %d %d", &size, &resident, &shared)==3) {
    //printf("%10.2f MB total size, %10.2f MB resident, %10.2f MB shared\n",
    //     (size*4)/1000.0, (resident*4)/1000.0, (shared*4)/1000.0);
    memusage = (size*4)/1000.0;  // ??
  }
  fclose(f);
#else
  static struct rusage rus;
  getrusage(RUSAGE_SELF, &rus);

  // ru_maxrss is a long int and is the "Maximum resident set size (in kilobytes)"

  //printf(" ru_maxrss= %d  %10.2fMB\n",
  // rus.ru_maxrss, rus.ru_maxrss/1000.0);
  //printf(" ru_ixrss=%d ru_idrss=%d ru_isrs=%d ru_maxrss=%d\n",
  //     rus.ru_ixrss, rus.ru_idrss, rus.ru_isrss, rus.ru_maxrss);
  memusage = rus.ru_maxrss/1000.0;
#endif

  return memusage;

}

#if defined(TIMER) && defined(CH_MPI)

#include <mpi.h>

// input is the "end" memory (like the peak) and the output
// is the avg, min, and max memory across all processors
void gather_memory_from_procs(Real end_memory,
                              Real &avg_memory,
                              Real &min_memory,
                              Real &max_memory) {

  // Depending on the mechanism used in get_memory_usage() to report
  // memory use, the variable end_memory will either be the peak
  // memory used at any point in the code execution (as getrusage()
  // works on OSF machine such as halem) or an estimate of the amount
  // of memory currently being used (which will certainly not be the
  // peak).  If end_memory is the peak memory used so far, then it is
  // useful in parallel jobs to gather these peak memory usage values
  // from all processors and obtain the minimum, maximum, and average
  // values.  This is what the following code does. (ndk)
  int result;
  Real sum_memory;
  result = MPI_Reduce(&end_memory, &sum_memory, 1, MPI_CH_REAL,
                      MPI_SUM, 0, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in getting memorys");
    }

  result = MPI_Reduce(&end_memory, &min_memory, 1, MPI_CH_REAL,
                      MPI_MIN, 0, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in getting memorys");
    }

  result = MPI_Reduce(&end_memory, &max_memory, 1, MPI_CH_REAL,
                      MPI_MAX, 0, Chombo_MPI::comm);
  if (result != MPI_SUCCESS)
    {
      MayDay::Error("Sorry, but I had a communication error in getting memorys");
    }

  int rank, number_procs;
  MPI_Comm_rank(Chombo_MPI::comm, &rank);
  MPI_Comm_size(Chombo_MPI::comm, &number_procs);

  // results only need to be on rank0 (and are)
  if (rank==0)
    {
      Real avg_memory = sum_memory/(Real)number_procs;

      pout() << "Gather end memory from procs:  avg: " << avg_memory
             << "  min: " << min_memory
             << "  max: " << max_memory << " (MB)" << endl;

    }
}
#endif
