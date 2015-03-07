#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
#include "DotProdF_F.H"
#include "AMRLevelTools.H"
#include "LayoutIterator.H"
#include "SPMD.H"

/// These helper functions are used by VCAMRSolver (Yuri Omelchenko)
Real 
AMRDotProduct(BoxLayoutData<FArrayBox>& a_dataOne, 
              BoxLayoutData<FArrayBox>& a_dataTwo,
              const BoxLayout& a_grids,
              const BoxLayout *a_gridsCFPtr,
              const Interval& a_comps)
{
  CH_assert(a_dataOne.nComp() == a_dataTwo.nComp());
  DataIterator dit= a_dataOne.dataIterator();
  bool self_product = (&a_dataOne == &a_dataTwo);
  bool has_finer_level = (a_gridsCFPtr != NULL);
  const int startcomp = a_comps.begin();
  const int endcomp = a_comps.end();
  const int ncomp = a_comps.size();

  Real rhodot = 0.0;
  // zero out data in invalid regions 
  for(dit.begin(); dit.ok(); ++dit)   // over grid data (parallel)
    {
     FArrayBox& onefab = a_dataOne[dit()];
     FArrayBox& twofab = a_dataTwo[dit()];
     Box fabbox= a_grids.get(dit());
     CH_assert(onefab.box().contains(fabbox));
     CH_assert(twofab.box().contains(fabbox));

     if(has_finer_level) 
     {
       LayoutIterator lit= a_gridsCFPtr->layoutIterator();
       for(lit.begin(); lit.ok(); ++lit) 
       {
        Box overlayBox = a_gridsCFPtr->get(lit());
        overlayBox &= fabbox;
        if(!overlayBox.isEmpty()) // "invalid data" region
          {
          onefab.setVal(0.,overlayBox,startcomp,ncomp);
          if(!self_product) twofab.setVal(0.,overlayBox,startcomp,ncomp);
          }
       } // endfor lit
     } // endif has_finer_level

     Real dotgrid= 0;
     FORT_DOTPRODUCT(CHF_REAL(dotgrid),
                     CHF_CONST_FRA(onefab), 
                     CHF_CONST_FRA(twofab),
                     CHF_BOX(fabbox),
                     CHF_CONST_INT(startcomp),
                     CHF_CONST_INT(endcomp));
     rhodot += dotgrid;
    } // endfor dit

  // now for the multi-processor fandango

  //gather all the rhodots onto a vector and add them up
  int baseProc = 0;
  Vector<Real> dotVec;
  gather(dotVec, rhodot, baseProc);

  Real rhodotTot = 0.0;
  if(procID() == baseProc)
    {
      CH_assert(dotVec.size() == numProc());
      for(int ivec = 0; ivec < dotVec.size(); ivec++)
        {
          rhodotTot += dotVec[ivec];
        }
    }
  //broadcast the sum to all processors.
  broadcast(rhodotTot, baseProc);

  //return the total
  return rhodotTot;
}

///
/// zero out data in "invalid" data regions 
/// (covered by finer grids)
///

void 
zeroAMRData(BoxLayoutData<FArrayBox>& a_data, 
            const BoxLayout *a_gridsCFPtr,
            const Interval& a_comps)
{
  if(a_gridsCFPtr == NULL) return; // no finer grid

  DataIterator dit= a_data.dataIterator();
  const BoxLayout& grids = a_data.boxLayout();
  const int startcomp = a_comps.begin();
  const int ncomp = a_comps.size();

  LayoutIterator lit= a_gridsCFPtr->layoutIterator();
  for(dit.begin(); dit.ok(); ++dit)   // over grid data (parallel)
    {
    FArrayBox& fab = a_data[dit()];
    Box fabbox= grids.get(dit());
    for(lit.begin(); lit.ok(); ++lit) 
      {
      Box overlayBox = a_gridsCFPtr->get(lit());
      overlayBox &= fabbox;
      if(!overlayBox.isEmpty()) fab.setVal(0.0,overlayBox,startcomp,ncomp);
      } // endfor lit
    } // endfor dit
}
