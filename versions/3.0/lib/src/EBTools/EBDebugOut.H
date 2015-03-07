#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBDEBUGOUT_H_
#define _EBDEBUGOUT_H_

#include "DisjointBoxLayout.H"
#include "BoxLayout.H"
#include "BaseIFFAB.H"
#include "BaseIVFAB.H"
#include "FaceIterator.H"
#include "Vector.H"
#include "Box.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "LayoutIterator.H"
#include "BoxIterator.H"
#include "VolIndex.H"
#include "Stencils.H"
#include "NamespaceHeader.H"

class EBCellFAB;

//extern "C"
//{
///
void
dumpLDEBCF(const LevelData<EBCellFAB>*  ldptr);

///
void
getMaxEBLevel(const LevelData<EBCellFAB>*  ldptr);

///
void
getMaxEBFAB(const EBCellFAB*  ldptr);


/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpVVoF(const Vector<VolIndex>* a_vectPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpVFace(const Vector<FaceIndex>* a_vectPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpFace(const FaceIndex* a_vectPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpFaceSten(const FaceStencil* a_stenPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpVoFSten(const   VoFStencil* a_stenPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpIFFAB(const BaseIFFAB<Real>* a_stenPtr);

///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpIVFAB(const BaseIVFAB<Real>* a_stenPtr);


///
/**
   Dump inputs to standard out.  For use inside debuggers.
*/
void dumpLDBIVF(const LayoutData< BaseIVFAB<Real> >* a_stenPtr);

///
/**
 */
void dumpEBLDDBL(const LevelData<EBCellFAB>* a_stenPtr);

///
/**
   dumps the fab to std out. 
*/
void dumpEBFAB(const EBCellFAB* a_fab);

///
/**
   dumps the irreg bits of the fab to std out. 
*/
void dumpEBFABIrreg(const EBCellFAB* a_fab);

///
/**
   dumps the irreg bits of the level to std out. 
*/
void dumpEBLevelIrreg(const LevelData<EBCellFAB>* a_fab);

///
/**
*/
void dumpEBFABIVS(const EBCellFAB* a_fab, const IntVectSet* a_ivs);

///
/**
*/
void dumpEBLevel(const LevelData<EBCellFAB>* a_level);
  
//}

#include "NamespaceFooter.H"
#endif
