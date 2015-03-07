#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _FABVIEW_H_
#define _FABVIEW_H_

#include "AMRIO.H"
#include "NamespaceHeader.H"

int InitWriteFAB()
{
#ifdef CH_USE_HDF5
  writeFAB(NULL);

  viewFAB(NULL);

  browseFAB(NULL);

  browseLevel(NULL);

  writeBFR(NULL);

  viewBFR(NULL);

  writeFABname(NULL, "label");

  writeLevel(NULL);

  viewLevel(NULL);

  viewLevelNoGhost(NULL);

  writeLevelname(NULL, "label");

  writeDBL(NULL);

  viewDBL(NULL);

  writeDBLname(NULL, "label");

  viewIVS(NULL, NULL);

  viewVectorBox(NULL, NULL);
#endif // CH_USE_HDF5

  return 1;
}

#include "NamespaceFooter.H"
#endif
