#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// AMRNodeLevelMG.H
// adapted from AMRLevelMG by DTGraves, Tues, July 6, 1999
// petermc, Tues, Nov 28, 2000

#ifndef AMRNODELEVELMG_H
#define AMRNODELEVELMG_H

#include <stdlib.h>
#include <iostream>
#include "SPACE.H"
#include <math.h>
#include "REAL.H"
#include "IntVect.H"
#include "ProblemDomain.H"
#include "NodeFArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "NodeLevelMG.H"
#include "NodeMGInterp.H"

#include "UsingNamespace.H"

// forward declaration of AMRNodeSolver
// C++ism:  don't need to #include AMRNodeSolver.H unless you have
// a data field that's an AMRNodeSolver
class AMRNodeSolver;

/// Class which manages data and operations on a level for AMRNodeSolver
class AMRNodeLevelMG
/** AMRNodeLevelMG manages the operations on a level for AMRNodeSolver.
    AMRNodeSolver contains a vector of AMRNodeLevelMG's (one for each level).
    This class should be considered internal to AMRNodeSolver and should
    not be considered part of the Chombo API.
 */
{
public:
  friend class AMRNodeSolver;

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.  User must subsequently call define().
   */
  AMRNodeLevelMG();

  ///
  /** Constructor to define AMRNodeLevelMG from the inputs.
      Calls define() function with identical arguments.
   */
  AMRNodeLevelMG(const AMRNodeSolver* const a_parent,
                 int a_level,
                 const NodeLevelOp* const a_opin);

  ///
  /** Destructor.
   */
  ~AMRNodeLevelMG();

  ///
  /** Define AMRNodeLevelMG from the inputs.

      The AMRNodeLevelMG will only work for the
      input AMRNodeSolver.  A NodeLevelMG is created.
  */
  void define(/// parent AMRNodeSolver
              const AMRNodeSolver* const a_parent,
              /// level in the AMR hierarchy
              int a_level,
              /// pointer to operator
              const NodeLevelOp* const a_opin);

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /** Returns <tt>true</tt> if this object was created with the defining
      constructor or if define() has been called.
  */
  bool isDefined() const;

  ///
  /** Compute and return norm of internal data <i>m_resid</i>.
   */
  Real computeResidualNorm(/// which norm to take, or 0 for max norm
                           int a_normType) const;

  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  ///
  /** Set number of multigrid smoothings on way up V-cycle;
      Default is 4.
  */
  void setnumSmoothUp(int a_numSmoothUp);

  ///
  /** Set number of multigrid smoothings on way down V-cycle;
      Default is 4.
  */
  void setnumSmoothDown(int a_numSmoothDown);

  ///
  /** Set whether the solver gives output while solving.
      Default is <tt>true</tt>.
   */
  void setVerbose(bool a_verbose);

  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  ///
  /** Complete, 3-level operator.  Applies operator with
      inhomogeneous coarse-fine boundary conditions and
      inhomogeneous domain boundary conditions.

      As a side effect, this function modifies
      <i>a_phiLevel</i>[<i>m_level</i>] by:
      - projecting down from interior coarse nodes of
      <i>a_phiLevel</i>[<i>m_level</i>+1], if this finer level exists;
      - interpolating from <i>a_phiLevel</i>[<i>m_level</i>-1]
      at the coarse/fine boundary, if this coarser level exists;
      - applying the physical boundary conditions.
  */
  void applyAMROperator(/// L(phi) result of evaluation at level <i>m_level</i>
                        LevelData<NodeFArrayBox>& a_Lofphi,
                        /// pointers to solution, indexed by level; max Vector index must be at least <i>m_level</i>
                        Vector<LevelData<NodeFArrayBox> *>& a_phiLevel);

  ///
  /** Complete, 3-level gradient.  Finds gradient with
      inhomogeneous coarse-fine boundary conditions and
      inhomogeneous domain boundary conditions.

      As a side effect, this function modifies
      <i>a_phiLevel</i>[<i>m_level</i>] by:
      - projecting down from interior coarse nodes of
      <i>a_phiLevel</i>[<i>m_level</i>+1], if this finer level exists;
      - interpolating from <i>a_phiLevel</i>[<i>m_level</i>-1]
      at the coarse/fine boundary, if this coarser level exists;
      - applying the physical boundary conditions.
  */
  void applyAMRGradient(/// gradient of phi at level <i>m_level</i>
                        LevelData<NodeFArrayBox>& a_gradPhi,
                        /// pointers to solution, indexed by level; max Vector index must be at least <i>m_level</i>
                        Vector<LevelData<NodeFArrayBox> *>& a_phiLevel);

  ///
  /** Compute complete, 3-level residual (rhs - operator(phi))
      and put it into local data <i>m_resid</i>.

      <i>m_resid</i> is zeroed out and its interior nodes filled with
      the residual.
  */
  void computeAMRResidual(/// pointers to solution, indexed by level; max Vector index must be at least <i>m_level</i>
                          Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
                          /// pointers to right-hand side, indexed by level; max Vector index must be at least <i>m_level</i>
                          const Vector<LevelData<NodeFArrayBox> *>& a_rhsLevel);

  ///
  /** Sweep up multigrid V-cycle.

      We modify <i>a_phiLevel</i>[<i>m_level</i>] by replacing it with
      <i>m_phiSave</i> plus a correction.

      In succession:
      - <i>m_corr</i> is interpolated from next coarser level's <i>m_corr</i>;
      - <i>m_resid</i> is decremented by operator on <i>m_corr</i>;
      - <i>m_dcorr</i> is set to zero and smoothed with <i>m_resid</i> as right-hand side;
      - <i>m_corr</i> is incremented by <i>m_dcorr</i>;
      - <i>m_phiSave</i> is incremented by <i>m_corr</i>.
  */
  void upSweep(/// pointers to solution, indexed by level; max Vector index must be at least <i>m_level</i>
               Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
               /// pointers to right-hand side, indexed by level; max Vector index must be at least <i>m_level</i>
               const Vector<LevelData<NodeFArrayBox> *>& a_rhsLevel);

  ///
  /** Sweep down multigrid V-cycle.

      We modify <i>a_phiLevel</i>[<i>m_level</i>] by incrementing it
      by a correction.  In succession:
      - <i>m_phiSave</i> is set to the original input
      <i>a_phiLevel</i>[<i>m_level</i>];
      - <i>m_corr</i> is set to zero and smoothed with <i>m_resid</i>
      as right-hand side;
      - <i>a_phiLevel</i>[<i>m_level</i>] is incrementing by
      <i>m_corr</i>;
      - next coarser level's <i>m_corr</i> is initialized to zero;
      - next coarser level's <i>m_resid</i> is set to the residual
      at the coarser level, and then at the coarse interior nodes
      is overwritten with average of this level's
      <i>m_resid</i> - operator(<i>m_corr</i>).

      As a side effect, if there is a coarser level, we modify
      <i>a_phiLevel</i>[<i>m_level</i>-1] by:
      - projecting down from interior coarse nodes of
      <i>a_phiLevel</i>[<i>m_level</i>];
      - interpolating from <i>a_phiLevel</i>[<i>m_level</i>-2]
      at the coarse/fine boundary, if this coarser level exists;
      - applying the physical boundary conditions.
   */
  void downSweep(/// pointers to solution, indexed by level; max Vector index must be at least <i>m_level</i>
                 Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
                 /// pointers to right-hand side, indexed by level; max Vector index must be at least <i>m_level</i>
                 const Vector<LevelData<NodeFArrayBox> *>& a_rhsLevel);

  /*@}*/

protected:

  /** Returns norm of input.
  */
  Real computeNorm(/// input variable
                   const LevelData<NodeFArrayBox>& a_mfinput,
                   /// type of norm, or 0 for max norm
                   int a_normType) const;

  /** Use short V-cycle or NodeLevelMG's smoother to reduce
      high wave numbers of residual. */
  void smooth(/// solution at this level
              LevelData<NodeFArrayBox>& a_phi,
              /// right-hand side at this level
              const LevelData<NodeFArrayBox>& a_rhs);

  /** Project <i>a_phiFine</i> at coarse interior nodes to <i>a_phi</i>.
   */
  void projectFineInterior(LevelData<NodeFArrayBox>& a_phi,
                           const LevelData<NodeFArrayBox>& a_phiFine);

  /** Delete internal memory and set pointers to NULL.
   */
  void clearMemory();

  /** Set default values.  Does not deal with memory.
   */
  void setDefaultValues();

  /** the boss
   */
  const AMRNodeSolver* m_parent;

  /** NodeLevelMG level solver object to relax on this level
   */
  NodeLevelMG m_levelMG;

  /** for interpolation
   */
  NodeMGInterp m_mginterp;

  /** has this AMRNodeLevelMG been defined yet?
   */
  bool m_isDefined;

  /** flag for verbose output
   */
  bool m_verbose;

  /** the grids at this level, coarsened
   */
  DisjointBoxLayout m_coarsenedGrids;

  /** residual
   */
  LevelData<NodeFArrayBox> m_resid;

  // LevelData<NodeFArrayBox> m_residInterior;

  /** correction
   */
  LevelData<NodeFArrayBox> m_corr;

  /** the infamous correction to the correction
   */
  LevelData<NodeFArrayBox> m_dcorr;

  /** L(phi) on this level
   */
  LevelData<NodeFArrayBox> m_LofPhi;

  /** extra copy of solution
   */
  LevelData<NodeFArrayBox> m_phiSave;

  /** coarsened residual, on m_coarsenedGrids
   */
  LevelData<NodeFArrayBox> m_residCoarsened;

  /** index of this level
   */
  int m_level;

  /** mesh spacing at this level
   */
  Real m_dx;

  /** pointer to operator
   */
  NodeLevelOp* m_levelOpPtr;

  /** pointer to coarsening operator
   */
  NodeCoarseAverage m_averageOp;

  /** CELL-centered physical domain at this level
   */
  ProblemDomain m_domain;

  /** the grids at the next finer level
   */
  DisjointBoxLayout m_finerGrids;

  /** the grids at the next finer level, coarsened
   */
  DisjointBoxLayout m_coarsenedFineGrids;

  /** interior boundary nodes of grids at this level
   */
  LayoutData< Vector<IntVectSet> > m_IVSV;

  /** exterior boundary nodes of grids at this level
   */
  LayoutData< Vector<IntVectSet> > m_IVSVext;

  /** interior boundary nodes of the coarsened grids at this level
   */
  LayoutData< Vector<IntVectSet> > m_IVSVcoarsened;

  /** interior boundary nodes of the coarsened grids at next finer level
   */
  LayoutData< Vector<IntVectSet> > m_IVSVcoarsenedFine;

private:

  AMRNodeLevelMG& operator = (const AMRNodeLevelMG&);
  AMRNodeLevelMG(const AMRNodeLevelMG&);

};

#endif
