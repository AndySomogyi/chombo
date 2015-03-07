#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#ifndef _AMR_LEVELTOOLS_H_
#define _AMR_LEVELTOOLS_H_

#include "BoxLayout.H"
#include "BoxLayoutData.H"
#include "FArrayBox.H"

// These helper functions are used by VCAMRSolver (Yuri Omelchenko)

/// Calculate a dot product on valid regions only (not covered by finer grids)
extern Real 
AMRDotProduct(BoxLayoutData<FArrayBox>& a_dataOne, 
              BoxLayoutData<FArrayBox>& a_dataTwo,
              const BoxLayout& a_grids,
              const BoxLayout *a_gridsCFPtr,
              const Interval&  a_comps);


///  Zero out data on invalid regions  (not covered by finer grids)
extern void 
zeroAMRData(BoxLayoutData<FArrayBox>& a_data, 
            const BoxLayout *a_gridsCFPtr,
            const Interval&  a_comps);

#endif


