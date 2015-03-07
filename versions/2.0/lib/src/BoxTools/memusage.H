#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _MEMUSAGE_H_
#define _MEMUSAGE_H_

#include "REAL.H"
#include "NamespaceHeader.H"

///
/**
 */
void print_memory_line(const char *);

///
/**
 */
void print_memory_line(void);

///
/**
 Return the residentSetSize of process from either /proc/self/statm or getrusage(RUSAGE_SELF, &rus)
 Units should be MB.
 If /proc/self/statm is used (Linux only), then size is also returned.
 It's not clear what the difference is right now, but size is generally larger.
 */
void getMemoryUsageFromOS(Real& residentSetSize, Real& size);

///
/**
 Maintain backward compatibility in code.
 Return the residentSetSize of process from either /proc/self/statm or getrusage(RUSAGE_SELF, &rus)
 Units should be MB.
 */
Real get_memory_usage_from_OS(void);

///
/**
 */
void gather_memory_from_procs(Real end_memory,
                              Real &avg_memory,
                              Real &min_memory,
                              Real &max_memory);

#include "NamespaceFooter.H"

#endif
