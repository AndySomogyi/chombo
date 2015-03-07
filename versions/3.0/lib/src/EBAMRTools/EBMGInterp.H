#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// dtgraves fri aug 17, 2001

#ifndef _EBMGINTERP_H_
#define _EBMGINTERP_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "EBIndexSpace.H"
#include "EBStencil.H"
#include "EBLevelGrid.H"
#include "NamespaceHeader.H"

/// Piecewise constant interpolation.
/**
   Finedata += coarse data.
  */
class EBMGInterp
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBMGInterp();

  ///
  ~EBMGInterp();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     dblFine, dblCoar: The fine and coarse layouts
     of the data.\\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
  */
  EBMGInterp(const DisjointBoxLayout& dblFine,
             const DisjointBoxLayout& dblCoar,
             const EBISLayout& ebislFine,
             const EBISLayout& ebislCoar,
             const ProblemDomain& domainCoar,
             const int& nref,
             const int& nvar,
             const EBIndexSpace* ebisPtr,
             const IntVect& ghostCellsPhi);


  ///
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     dblFine, dblCoar: The fine and coarse layouts
     of the data.\\
     ebislFine, ebislCoar: The fine and coarse layouts
     of the geometric description.\\
     nref: The refinement ratio between the two levels. \\
     nvar: The number of variables contained in the data
     at each VoF.
  */
  void define(const DisjointBoxLayout& dblFine,
              const DisjointBoxLayout& dblCoar,
              const EBISLayout& ebislFine,
              const EBISLayout& ebislCoar,
              const ProblemDomain& domainCoar,
              const int& nref,
              const int& nvar,
              const EBIndexSpace* ebisPtr,
              const IntVect& ghostCellsPhi);

  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /**
     Piecewise constant interpolation function
     FineData += coarseData
  */
  void
  pwcInterp(LevelData<EBCellFAB>&       a_fineData,
            const LevelData<EBCellFAB>& a_coarseData,
            const Interval&             a_variables);

  ///
  /**
     Piecewise constant interpolation function for unchanged layouts (no agglomeration)
     there are no copies here
     FineData += coarseData
  */
  void
  pwcInterpMG(LevelData<EBCellFAB>&       a_fineData,
              const LevelData<EBCellFAB>& a_coarseData,
              const Interval&             a_variables);

protected:
  void
  pwcInterpFAB(EBCellFAB&       a_refCoar,
//                EBCellFAB&       a_refCoarSave,
               const EBCellFAB& a_fine,
               const DataIndex& a_datInd,
               const Interval&  a_variables) const;

  void
  setDefaultValues();

  bool m_isDefined;
  IntVect m_ghost;
  DisjointBoxLayout m_coarGrids;
  DisjointBoxLayout m_fineGrids;
  DisjointBoxLayout m_refinedCoarseGrids;
  Copier m_copierFtoRC;//this goes from m_fineGrids to m_refinedCoarseGrids
  Copier m_copierRCtoF;//this goes from m_refinedCoarseGrids to m_fineGrids
  ProblemDomain m_coarDomain;
  ProblemDomain m_fineDomain;

  EBISLayout m_coarEBISL;
  EBISLayout m_fineEBISL;

  EBISLayout m_refinedCoarseEBISL;

  int m_refRat;
  int m_nComp;
  LevelData<EBCellFAB> m_refinedCoarseData;

  //stencils for prolongation
  LayoutData<RefCountedPtr<EBStencil> >  m_interpEBStencil;

private:
    //disallowed for all the usual reasons
  EBMGInterp(const EBMGInterp& ebcin)
  {
    MayDay::Error("ebmgi 2 invalid operator");
  }
  void operator=(const EBMGInterp& fabin)
  {
    MayDay::Error("ebmgi 3 invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
