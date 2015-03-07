#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _EBLEVELTGA_H_
#define _EBLEVELTGA_H_

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include "REAL.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBLevelGrid.H"
#include "EBFluxFAB.H"
#include "EBFluxRegister.H"
#include "ProblemDomain.H"
#include "BaseLevelTGA.H"
#include "NamespaceHeader.H"


/// Implements TGA algorithm to update solution to heat equation
/** The LevelTGA class implements the Runge-Kutta-based approach
    to solving the heat equation due to Twizell, Gumel, and Arigu.
*/
class EBLevelTGA : public BaseLevelTGA<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>
{

public:

  /// full constructor
  /**
   */
  EBLevelTGA(const Vector<DisjointBoxLayout>&                               a_grids,
             const Vector<int>&                                             a_refRat,
             const ProblemDomain&                                           a_level0Domain,
             RefCountedPtr<AMRLevelOpFactory<LevelData<EBCellFAB> > >&      a_opFact,
             const RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >&     a_solver)
    :BaseLevelTGA<LevelData<EBCellFAB>, EBFluxFAB, EBFluxRegister>(a_grids, a_refRat, a_level0Domain, a_opFact, a_solver)
  {
    m_isEBLGSet = false;
  }

  /// destructor
  virtual ~EBLevelTGA()
  {
  }

  /// computes time-centered diffusive term for explicit updating
  /** compute time-centered L(phi) for use in subsequent update
      operations.  In this case, we do solve for phiNew, then
      subtract source and old phi back out to get L(phi).
  */
  void computeDiffusion(LevelData<EBCellFAB>& a_DiffusiveTerm,
                        LevelData<EBCellFAB>& a_phiOld,
                        LevelData<EBCellFAB>& a_src,
                        EBFluxRegister* a_FineFluxRegPtr,
                        EBFluxRegister* a_crseFluxRegPtr,
                        const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                        const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                        Real a_oldTime,
                        Real a_crseOldTime,
                        Real a_crseNewTime,
                        Real a_dt,
                        int a_level);

  /// do solve
  /** update phi to a_phiNew
   */
  void updateSoln(LevelData<EBCellFAB>& a_phiNew,
                  LevelData<EBCellFAB>& a_phiOld,
                  LevelData<EBCellFAB>& a_src,
                  EBFluxRegister* a_fineFluxRegPtr,
                  EBFluxRegister* a_crseFluxRegPtr,
                  const LevelData<EBCellFAB>* a_crsePhiOldPtr,
                  const LevelData<EBCellFAB>* a_crsePhiNewPtr,
                  Real a_oldTime,
                  Real a_crseOldTime,
                  Real a_crseNewTime,
                  Real a_dt,
                  int a_level,
                  bool a_zeroPhi = true,
                  int a_fluxStartComponent = 0);

  void setSourceGhostCells(LevelData<EBCellFAB>&    a_src,
                           const DisjointBoxLayout& a_grids);

  void setEBLG(Vector<EBLevelGrid>& a_eblg)
  {
    m_isEBLGSet = true;
    m_eblg      = a_eblg;
  }
protected:
  bool m_isEBLGSet;
  Vector<EBLevelGrid> m_eblg;

};

#include "NamespaceFooter.H"
#endif


