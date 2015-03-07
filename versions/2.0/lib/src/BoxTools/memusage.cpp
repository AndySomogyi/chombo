#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "memusage.H"
#include "memtrack.H"
#include "parstream.H"
#include "MayDay.H"
#include "ParmParse.H"
#include "CH_HDF5.H"
#include <cstdio>
#include "SPMD.H"

#ifdef CH_MPI
#include <mpi.h>
#endif

// note:  i copied the following #include crud from elsewhere.
// assumed it was still relevant (ndk)
#include <sys/time.h>
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

#include "NamespaceHeader.H"

void print_memory_line(void) 
{
  std::string s="";
  print_memory_line(s.data());
}

void print_memory_line(const char *s)
{
  char temp[180];
#ifdef CH_USE_MEMORY_TRACKING
  Real memtrackCurrentMemory;
  Real memtrackPeakMemory;
  memtrackStamp(memtrackCurrentMemory, memtrackPeakMemory);

  sprintf(temp, "%26s|Mem Usage: OS=%8.3f  MT_peak=%8.3f  MT_current=%8.3f (MB)\n",
          s, get_memory_usage_from_OS(), memtrackPeakMemory, memtrackCurrentMemory);
#else
  sprintf(temp, "%26s|Mem Usage: OS=%8.3f (MB)  MT is off\n",
          s, get_memory_usage_from_OS());
#endif
  pout() << temp;
}

void getMemoryUsageFromOS(Real& residentSetSize, Real& size)
{
  size=0.0;
  residentSetSize=0.0;
#ifndef CH_DISABLE_SIGNALS
  // if Linux(except for catamount), look at the file /proc/self/statm
#ifdef CH_Linux
  FILE* f=fopen("/proc/self/statm","r");
  int isize, iresident, ishared;
  if(fscanf(f, "%d %d %d", &isize, &iresident, &ishared)==3) 
    {
      //printf("%10.2f MB total size, %10.2f MB resident, %10.2f MB shared\n",
      //     (isize*4)/1024.0, (iresident*4)/1024.0, (ishared*4)/1024.0);

      // linux on intel x86 has 4KB page size...
      //return ($size * 4, $share * 4);

      size = (isize*4)/1024.0;
      residentSetSize = (iresident*4)/1024.0;
    }
  fclose(f);
#else
  static struct rusage rus;
  getrusage(RUSAGE_SELF, &rus);

  // ru_maxrss is a long int and is the "Maximum resident set size (in kilobytes)"

  //printf(" ru_maxrss= %d  %10.2fMB\n",
  // rus.ru_maxrss, rus.ru_maxrss/1024.0);
  //printf(" ru_ixrss=%d ru_idrss=%d ru_isrs=%d ru_maxrss=%d\n",
  //     rus.ru_ixrss, rus.ru_idrss, rus.ru_isrss, rus.ru_maxrss);
  residentSetSize = rus.ru_maxrss/1024.0;
#endif
#endif
}

// Maintain backward compatibility in code.
// Return the residentSetSize of process from either /proc/self/statm or getrusage(RUSAGE_SELF, &rus)
// Units should be MB.
Real get_memory_usage_from_OS(void) 
{
  Real residentSetSize=0.0;
  Real size;
  getMemoryUsageFromOS(residentSetSize, size);  // units MB
  return residentSetSize;
}


#ifdef CH_MPI

// input is the "end" memory (like the peak) and the output
// is the avg, min, and max memory across all processors
void gather_memory_from_procs(Real end_memory,
                              Real &avg_memory,
                              Real &min_memory,
                              Real &max_memory) 
{

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

#include "NamespaceFooter.H"
