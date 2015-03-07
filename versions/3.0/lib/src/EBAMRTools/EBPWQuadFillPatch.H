#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBPWQUADFILLPATCH_H_
#define _EBPWQUADFILLPATCH_H_

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "DisjointBoxLayout.H"
#include "EBISLayout.H"
#include "EBCellFAB.H"
#include "Interval.H"
#include "Stencils.H"
#include "BaseIVFAB.H"
#include "EBIndexSpace.H"
#include "EBArith.H"
#include "CornerCopier.H"
#include "NamespaceHeader.H"

///
/**
   Fills one layer of fine level ghost cells
   by piecewise quadratic interpolation
   from the coarse level.
 */
class EBPWQuadFillPatch
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBPWQuadFillPatch();

  ///
  ~EBPWQuadFillPatch();

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
     at each VoF.\\
     radius: the total number of
     ghost cells filled.
  */
  EBPWQuadFillPatch(const DisjointBoxLayout& a_dblFine,
                 const DisjointBoxLayout& a_dblCoar,
                 const EBISLayout& a_ebislFine,
                 const EBISLayout& a_ebislCoar,
                 const ProblemDomain& a_domainCoar,
                 const int& a_nref,
                 const int& a_nvar,
                 const int& a_radius);

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
     radius: the total number of
     ghost cells filled.
  */
  void define(const DisjointBoxLayout& a_dblFine,
              const DisjointBoxLayout& a_dblCoar,
              const EBISLayout& a_ebislFine,
              const EBISLayout& a_ebislCoar,
              const ProblemDomain& a_domainCoar,
              const int& a_nref,
              const int& a_nvar,
              const int& a_radius,
              const EBIndexSpace* const a_eb = Chombo_EBIS::instance());

  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /**
     Interpolate the fine data from the coarse data
     over the intersection
     of the fine layout with the refined coarse layout. \\
     {\bf Arguments:}\\
     coarDataOld: The data over the coarse layout at coarse time old.\\
     coarDatanew: The data over the coarse layout at coarse time new.\\
     fineData: The data over the fine layout.  \\
     Error occurs unless  coarTimeOld <= fineTime <= coarTimeNew \\
     Fine and coarse data must
     have the same number of variables.\\
     variables:  The variables to interpolate.  Those not
     in this range will be left alone.   This range of variables
     must be in both the coarse and fine data.
  */
  void
  interpolate(LevelData<EBCellFAB>& a_fineData,
              const LevelData<EBCellFAB>& a_coarDataOld,
              const LevelData<EBCellFAB>& a_coarDataNew,
              const Real& a_coarTimeOld,
              const Real& a_coarTimeNew,
              const Real& a_fineTime,
              const Interval& a_variables);

protected:

  //internal use only functions
  void
  interpolateFAB(EBCellFAB& a_fine,
                 const EBCellFAB& a_coarOld,
                 const EBCellFAB& a_coarNew,
                 const Real& a_coarTimeOld,
                 const Real& a_coarTimeNew,
                 const Real& a_fineTime,
                 const DataIndex& a_datInd,
                 const Interval& a_variables) const;

  void
  makeStencils();

  void setDefaultValues();

  bool m_isDefined;

  bool extractFromSten(Real& a_value,
                       const VoFStencil& a_vofsten,
                       const EBCellFAB& a_coarData,
                       const int& ivar) const;

  Real computeDMinMod(const BaseIVFAB<VoFStencil >& a_lostenBF,
                      const BaseIVFAB<VoFStencil >& a_histenBF,
                      const VolIndex& a_coarVoF,
                      const EBCellFAB& a_coarData,
                      const int& ivar,
                      const IntVectSet& loInterpSet,
                      const IntVectSet& hiInterpSet,
                      const IntVectSet& ceInterpSet) const;

  //bloated class data
  ProblemDomain m_coarDomain;
  DisjointBoxLayout m_coarsenedFineGrids;
  DisjointBoxLayout m_fineGrids;
  DisjointBoxLayout m_coarGrids;
  EBISLayout m_coarsenedFineEBISL;
  EBISLayout m_fineEBISL;

  LevelData<EBCellFAB> m_coarsenedFineData;
  LevelData<EBCellFAB> m_coarOnFDataOld;
  LevelData<EBCellFAB> m_coarOnFDataNew;

  int m_refRat;
  int m_nComp;
  int m_radius;
  int m_coarGhostRad;

  //places to do interpolation
  //at fine refinement level
  LayoutData<IntVectSet> m_irregRegionsFine;

  //diffs high and low---need to make the choice
  //about which to use based on the data.   one
  //for each direction of the derivative.
  LayoutData<BaseIVFAB<VoFStencil> >  m_hiStencils[SpaceDim];
  LayoutData<BaseIVFAB<VoFStencil> >  m_loStencils[SpaceDim];

  // per-grid coarse locations that you interpolate from, by type of
  // interpolation in the specified coordinate direction.
  //this specifies constraints due to domain and union of rectangles
  //at coarse level
  LayoutData<IntVectSet> m_coarCeInterp[SpaceDim];
  LayoutData<IntVectSet> m_coarLoInterp[SpaceDim];
  LayoutData<IntVectSet> m_coarHiInterp[SpaceDim];

  //corner copier for the corner exchange case
  CornerCopier m_cornerCopier;

private:
  //disallowed for all the usual reasons
  EBPWQuadFillPatch(const EBPWQuadFillPatch& ebcin)
  {
    MayDay::Error("ebpwl 2 invalid operator");
  }
  void operator=(const EBPWQuadFillPatch& fabin)
  {
    MayDay::Error("ebpwl 3 invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
