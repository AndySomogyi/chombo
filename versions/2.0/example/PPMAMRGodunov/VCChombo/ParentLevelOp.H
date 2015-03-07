#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif
// DTGraves, Weds, July 21, 1999

#ifndef _PARENTLEVELOP_H_
#define _PARENTLEVELOP_H_

#include <iostream>
#include <cmath>
#include <assert.h>
#include <cstdlib>
#include <REAL.H>
#include <IntVect.H>
#include <Box.H>
#include <DisjointBoxLayout.H>
#include <LevelData.H>
#include <BaseFab.H>
#include <FArrayBox.H>
#include "GhostBC.H"
#include "LevelOp.H"
#include "QuadCFInterp.H"
#include "CFIVS.H"
#include "Copier.H"
#include "BaseBottomSmoother.H"

///
/** this class contains stuff needed by derived operators 
    Yuri Omelchenko */
class ParentLevelOp: public LevelOp
{
public:

  ///
  /**
     virtual constructor workaround
  */
  virtual
  LevelOp* new_levelop() const = 0;

  ///
  virtual
  bool isDefined() const;

  ///
  /**
     Full define function.  makes all coarse-fine
     information and sets internal variables 
  */
  virtual
  void define(
              const DisjointBoxLayout& Ba,
              const DisjointBoxLayout* base_ba,
              Real  DxLevel, 
              int refratio,
              const Box& domf,
              bool a_homogeneousOnly=false,
              int ncomp = 1);

  ///
  /**
     Full define function.  makes all coarse-fine
     information and sets internal variables 
  */
  virtual
  void define(
              const DisjointBoxLayout& Ba,
              const DisjointBoxLayout* base_ba,
              Real  DxLevel, 
              int refratio,
              const ProblemDomain& domf,
              bool a_homogeneousOnly=false,
              int ncomp = 1);

  ///
  /**
     Full define function.  define
     from finer levelop.
  */
  virtual
  void define(
              const LevelOp* opfine,
              int refratio);

  ///
  ParentLevelOp();

  ///
  virtual ~ParentLevelOp();

  /**
     Coarse / Fine interpolation operator.
  */
  virtual
  void CFInterp(LevelData<FArrayBox>& a_phi,
                const LevelData<FArrayBox>& a_phiCoarse
                );

  /// 
  /**
     homogeneous coarse/fine interpolation operator
  */
  virtual
  void homogeneousCFInterp(LevelData<FArrayBox>& a_phi);

  ///
  virtual
  void setDomainGhostBC(const DomainGhostBC& a_dombcin);

  ///
  virtual
  void setBottomSmoother(const BaseBottomSmoother& a_bottomSmoother);

  ///
  /**
     does homogeneous coarse/fine interpolation
  */
  virtual
  void homogeneousCFInterp(LevelData<FArrayBox>& a_phif, 
                           const DataIndex& a_datInd, 
                           int a_idir, 
                           Side::LoHiSide a_hiorlo);
  ///
  /**
     This does a GSRB Pre/Conditioned BiCGStab on a level
     for the bottom solver.  
  */
  virtual
  void bottomSmoother(
                      LevelData<FArrayBox>& a_phi,
                      const LevelData<FArrayBox>& a_rhs);

  /// set convergence metrics (passing through to bottom smoother)
  void setConvergenceMetric(Real a_metric, int a_comp);
  
protected:

  ///
  virtual
  void setDefaultValues();

  ///
  virtual
  void clearMemory();

  /// boundary conditons
  DomainGhostBC m_domghostbc;

  /// domain of fine grid
  ProblemDomain m_domain;

  /// the grids at the current level
  DisjointBoxLayout m_grids;

  /// 
  /** Copier object that knows how to perform exchange operations on 
      LevelData objects defined on m_grids */
  Copier m_exchangeCopier;

  //
  /** DisjointBoxLayout for next coarser level
    this is == NULL if there is no coarser level */
  DisjointBoxLayout m_baseBA;

  /// number of components
  int m_ncomp;

  /// refinement ratio between this and the next coarser level
  int m_refRatio;

  /// mesh spacing at this level
  Real m_dxLevel;

  /// mesh spacing at next coarser level
  Real m_dxCrse;

  ///
  /**
    coarse-fine information between this level
    and the next coarser level
  */
  QuadCFInterp m_quadCFI;

  ///
  /* bottom smoother object
   */
  BaseBottomSmoother *m_bottom_smoother_ptr;

  /// has full define function been called?
  bool m_isDefined;

  /// has bc been defined?
  bool m_isBCDefined;

  /// is inhomogeneous CF interpolation possible?
  bool m_ihcfiEnabled;

  ///
  LayoutData<CFIVS> m_loCFIVS[SpaceDim];
  ///
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

private:
  ///
  virtual
  void interpOnIVS(
                  LevelData<FArrayBox>& a_phif, 
                  const FArrayBox& a_phistar, 
                  const DataIndex& dFine,
                  const int a_idir,
                  const Side::LoHiSide a_hiorlo,
                  const IntVectSet& a_interpIVS);

};

#endif
