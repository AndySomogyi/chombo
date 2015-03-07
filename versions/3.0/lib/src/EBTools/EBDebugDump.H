#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBDEBUGDUMP_H_
#define _EBDEBUGDUMP_H_

#include "EBDebugOut.H"
#include "NamespaceHeader.H"

void EBDontYouDareCallMe()
{
  dumpVVoF(NULL);
  dumpIFFAB(NULL);
  dumpIVFAB(NULL);
  dumpLDBIVF(NULL);
  dumpFace(NULL);
  dumpVFace(NULL);
  dumpVoFSten(NULL);
  dumpFaceSten(NULL);
  dumpEBLDDBL(NULL);
  dumpLDEBCF(NULL);
  getMaxEBLevel(NULL);
  getMaxEBFAB(NULL);
  dumpEBFAB(NULL);
  dumpEBFABIrreg(NULL);
  dumpEBLevelIrreg(NULL);
  dumpEBLevel(NULL);
  dumpEBFABIVS(NULL, NULL);
}
#include "NamespaceFooter.H"
#endif
