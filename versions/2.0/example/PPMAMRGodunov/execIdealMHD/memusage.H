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

Real get_memory_usage(void);

void gather_memory_from_procs(Real end_memory,
                              Real &avg_memory,
                              Real &min_memory,
                              Real &max_memory);
#endif
