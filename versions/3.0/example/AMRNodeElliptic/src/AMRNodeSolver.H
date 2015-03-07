#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// AMRNodeSolver.H
// adapted from AMRSolver by DTGraves, Tues, July 6, 1999
// petermc, Wed, Nov 28, 2000

#ifndef AMRNODESOLVER_H
#define AMRNODESOLVER_H

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "SPACE.H"
#include "REAL.H"
#include "IntVect.H"
#include "ProblemDomain.H"
#include "NodeFArrayBox.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "AMRNodeLevelMG.H"
#include "LevelNodeSolver.H"

#include "UsingNamespace.H"

/// Class which manages grid hierarchy and composite elliptic solution
/**
   AMRNodeSolver manages the AMR / multigrid solution to the elliptic
   equation on a set of grids on multiple levels satisfying certain
   proper-nesting conditions.

   It can be used either as a solver, or to apply the AMR / multigrid
   V-cycle as a preconditioner for some other iterative method, such as
   a Krylov method.  */
class AMRNodeSolver
{
public:

  friend class AMRNodeLevelMG;

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Creates an AMRNodeSolver whose internal state is undefined.
      Need to call define() function to use any of
      the functionality of the class.
  */
  AMRNodeSolver();

  ///
  /** Creates a fully-defined AMRNodeSolver.   Calls
      define() function with identical arguments.
  */
  AMRNodeSolver(const Vector<DisjointBoxLayout>& a_gridsLevel,
                const Vector<ProblemDomain>&     a_domainLevel,
                const Vector<Real>&              a_dxLevel,
                const Vector<int>&               a_refRatio,
                int                              a_numLevels,
                int                              a_lBase,
                const NodeLevelOp* const         a_opin,
                int                              a_minLength = 1);

  ///
  /** Creates a fully-defined AMRNodeSolver.   Calls
      define() function with identical arguments.
  */
  AMRNodeSolver(const Vector<DisjointBoxLayout>& a_gridsLevel,
                const Vector<Box>&               a_domainLevel,
                const Vector<Real>&              a_dxLevel,
                const Vector<int>&               a_refRatio,
                int                              a_numLevels,
                int                              a_lBase,
                const NodeLevelOp* const         a_opin,
                int                              a_minLength = 1);

  ///
  /** Destructor.
   */
  ~AMRNodeSolver();

  ///
  /** Defines AMRNodeSolver's internal state.

      Except for <i>a_refRatio</i>, each Vector argument contains
      a component for each level, and Vector index corresponds to
      level number.
  */
  void define(/// the grids (CELL-centered) at all levels; Vector length at least <i>a_numLevels</i>
              const Vector<DisjointBoxLayout>& a_gridsLevel,
              /// the CELL-centered physical domains at all levels; Vector length at least <i>a_numLevels</i>
              const Vector<ProblemDomain>&     a_domainLevel,
              /// the grid spacing at all levels; Vector length at least <i>a_numLevels</i>
              const Vector<Real>&              a_dxLevel,
              /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
              const Vector<int>&               a_refRatio,
              /// the number of AMR levels in the calculation
              int                              a_numLevels,
              /// index of coarsest level on which solution is to be computed; must be set at time of definition in order to build the bottom LevelNodeSolver
              int                              a_lBase,
              /// pointer to the NodeLevelOp to use in the solution
              const NodeLevelOp* const         a_opin,
              /// minimum length of maximally coarsened box in LevelNodeSolver, or 0 for no coarsening
              int                              a_minLength = 1);


  ///
  /** Defines AMRNodeSolver's internal state.

      Except for <i>a_refRatio</i>, each Vector argument contains
      a component for each level, and Vector index corresponds to
      level number.
  */
  void define(/// the grids (CELL-centered) at all levels; Vector length at least <i>a_numLevels</i>
              const Vector<DisjointBoxLayout>& a_gridsLevel,
              /// the CELL-centered physical domains at all levels; Vector length at least <i>a_numLevels</i>
              const Vector<Box>&               a_domainLevel,
              /// the grid spacing at all levels; Vector length at least <i>a_numLevels</i>
              const Vector<Real>&              a_dxLevel,
              /// refinement ratios between adjacent levels, starting with <i>a_refRatio</i>[0], the refinement ratio between levels 0 and 1; Vector length at least <i>a_numLevels</i>-1
              const Vector<int>&               a_refRatio,
              /// the number of AMR levels in the calculation
              int                              a_numLevels,
              /// index of coarsest level on which solution is to be computed; must be set at time of definition in order to build the bottom LevelNodeSolver
              int                              a_lBase,
              /// pointer to the NodeLevelOp to use in the solution
              const NodeLevelOp* const         a_opin,
              /// minimum length of maximally coarsened box in LevelNodeSolver, or 0 for no coarsening
              int                              a_minLength = 1);

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /** Returns <tt>true</tt> if full define() function has been called.
   */
  bool isDefined() const;

  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  ///
  /** Set number of multigrid smoothings on way up V-cycle;
      Default is 4.
  */
  void setNumSmoothUp(int a_numSmoothUp);

  ///
  /** Set number of multigrid smoothings on way down V-cycle;
      Default is 4.
  */
  void setNumSmoothDown(int a_numSmoothDown);

  ///
  /** Set tolerance of iterative solution.  Default is 1.0e-10.
   */
  void setTolerance(Real a_tolerance);

  ///
  /** Set whether bottom solver does smoothing at bottom level.
      Default is set in LevelNodeSolver.
   */
  void setBottomSmoothing(bool a_doBottomSmooth);

  ///
  /** Set tolerance of bottom solver.  Default is set in LevelNodeSolver.
   */
  void setBottomTolerance(Real a_tolerance);

  ///
  /** Set "operator tolerance" of iterative solution.  Iteration will
      stop if (new_residual/old_residual) > 1-<i>a_operatorTolerance</i>
      (and at least <i>m_minIter</i> iterations have been performed).
      Default is 1.0e-5.
  */
  void setOperatorTolerance(Real a_operatorTolerance);

  ///
  /** Set max number of iterations.  Default is 42.
   */
  void setMaxIter(int a_maxIter);

  ///
  /** Set min number of iterations.  Only relevant when residual is
      not decreasing fast enough to satisfy the "operator tolerance".
      Default is 4.
  */
  void setMinIter(int a_minIter);

  ///
  /** Set max number of iterations for bottom solver.
      Default is set in LevelNodeSolver.
   */
  void setBottomMaxIter(int a_maxIter);

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
  /** Solves the elliptic equation over the hierarchy of levels
      <i>m_lBase</i> ... <i>m_finestLevel</i>
      where <i>m_finestLevel</i> = <i>m_numLevels</i>-1.

      If <i>m_lBase</i> > 0, then the data at level <i>m_lBase</i> - 1 is used
      to interpolate boundary conditions at boundary nodes that are not
      adjacent to the domain boundary.

      Solves to tolerance <i>m_tolerance</i>.
  */
  void solveAMR(/// pointers to current guess at the solution values for levels (lMin = max(<i>m_lBase</i>-1,0)) ... <i>m_finestLevel</i>.  Vector index corresponds to level number.  These values are updated in place.
                Vector<LevelData<NodeFArrayBox> *>&       a_phiLevel,
                /// pointers to right-hand side for levels lMin ... <i>m_finestLevel</i>.  Not modified.
                const Vector<LevelData<NodeFArrayBox> *>& a_rhsLevel);

  ///
  /** Does one relaxation V-cycle using a multigrid solver.
  */
  void AMRVCycleMG(/// pointers to current guess at the solution values for levels (lMin = max(<i>m_lBase</i>-1,0)) ... <i>m_finestLevel</i>.  Vector index corresponds to level number.  These values are updated in place.
                   Vector<LevelData<NodeFArrayBox> *>&       a_phiLevel,
                   /// pointers to right-hand side for levels lMin ... <i>m_finestLevel</i>.  Not modified.
                   const Vector<LevelData<NodeFArrayBox> *>& a_rhsLevel);

  ///
  /** Calculate norm of multilevel residual on levels
      <i>m_lBase</i> to <i>m_finestLevel</i>.
      Does not include data covered by finer levels (not valid data).
      Assumes that the residual has already been computed in
      <i>m_amrmgLevel</i>[<i>ilev</i>]-><i>m_resid</i> for all levels
      <i>ilev</i> from <i>m_lBase</i> to <i>m_finestLevel</i>.
  */
  Real computeResidualNorm(/// which norm to take, or 0 for max norm
                           int a_normType);


  ///
  /** Calculate residual on level <i>a_ilev</i> using the multilevel operator.
      The residual is rhs - L(phi).

      As a side effect, this function modifies <i>a_phiLevel</i>[<i>a_ilev</i>] by:
      - projecting down from interior coarse nodes of
      <i>a_phiLevel</i>[<i>a_ilev</i>+1], if this finer level exists;
      - interpolating from <i>a_phiLevel</i>[<i>a_ilev</i>-1]
      at the coarse/fine boundary, if this coarser level exists;
      - applying the physical boundary conditions.

      This function also stores the residual at level <i>a_ilev</i> in
      <i>m_amrmgLevel</i>[<i>a_ilev</i>]-><i>m_resid</i>.
   */
  void computeAMRResidual(/// residual, rhs - L(phi), at level <i>a_ilev</i>
                          LevelData<NodeFArrayBox>&                 a_res,
                          /// pointers to solution, indexed by level; max Vector index must be at least <i>a_ilev</i>
                          Vector<LevelData<NodeFArrayBox> *>&       a_phiLevel,
                          /// pointers to right-hand side, indexed by level; max Vector index must be at least <i>a_ilev</i>
                          const Vector<LevelData<NodeFArrayBox> *>& a_rhsLevel,
                          /// level at which the residual is evaluated
                          int                                       a_ilev);

  ///
  /** Calculate multilevel L(phi) at level <i>a_ilev</i>.
      This also depends on data from the next coarser level
      and the next finer level.

      As a side effect, this function modifies <i>a_phiLevel</i>[<i>a_ilev</i>] by:
      - projecting down from interior coarse nodes of
      <i>a_phiLevel</i>[<i>a_ilev</i>+1], if this finer level exists;
      - interpolating from <i>a_phiLevel</i>[<i>a_ilev</i>-1]
      at the coarse/fine boundary, if this coarser level exists;
      - applying the physical boundary conditions.
  */
  void applyAMROperator(/// L(phi) result of operator evaluation at level <i>a_ilev</i>
                        LevelData<NodeFArrayBox>&           a_lofPhi,
                        /// pointers to solution, indexed by level; max Vector index must be at least <i>a_ilev</i>
                        Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
                        /// level at which operator is evaluated
                        int                                 a_ilev);

  ///
  /** Calculate multilevel gradient at level <i>a_ilev</i>.
      This also depends on data from the next coarser level
      and the next finer level.

      As a side effect, this function modifies <i>a_phiLevel</i>[<i>a_ilev</i>] by:
      - projecting down from interior coarse nodes of
      <i>a_phiLevel</i>[<i>a_ilev</i>+1], if this finer level exists;
      - interpolating from <i>a_phiLevel</i>[<i>a_ilev</i>-1]
      at the coarse/fine boundary, if this coarser level exists;
      - applying the physical boundary conditions.
  */
  void applyAMRGradient(/// gradient of phi at level <i>a_ilev</i>
                        LevelData<NodeFArrayBox>&           a_gradPhi,
                        /// pointers to solution, indexed by level; max Vector index must be at least <i>a_ilev</i>
                        Vector<LevelData<NodeFArrayBox> *>& a_phiLevel,
                        /// level at which gradient is evaluated
                        int                                 a_ilev);

  /*@}*/

protected:

  /** Set data to default values.  Not for external use.
   */
  void setDefaultValues();

  /** Delete <i>m_amrmgLevel</i> Vector.
   */
  void clear();

  /** solver tolerance
   */
  Real m_tolerance;

  /** operator tolerance:  iteration will
      stop if (new_residual/old_residual) > 1-<i>m_operatorTolerance</i>
   */
  Real m_operatorTolerance;

  /** index of coarsest level on which solution is to be computed
   */
  int m_lBase;

  /** number of levels allowed
   */
  int m_numLevels;

  /** index of finest level used, equal to <i>m_numLevels</i>-1
   */
  int m_finestLevel;

  /** refinement ratio between levels; <i>m_refRatio</i>[ilev] is refinement ratio between levels ilev and ilev+1
   */
  Vector<int> m_refRatio;

  /** grid spacing at each level
   */
  Vector<Real> m_dxLevel;

  /** grids, CELL-centered (disjoint) at each level
   */
  Vector<DisjointBoxLayout> m_gridsLevel;

  /** CELL-centered physical domain at each level
   */
  Vector<ProblemDomain> m_domainLevel;

  /** maximum number of solver iterations
   */
  int m_maxIter;

  /** minimum number of solver iterations before "hung convergence" criterion applies
   */
  int m_minIter;

  /** number of smoother iterations on sweep up
   */
  int m_numSmoothUp;

  /** number of smoother iterations on sweep down
   */
  int m_numSmoothDown;

  /** data at each refinement level
   */
  Vector<AMRNodeLevelMG*> m_amrmgLevel;

  /** whether the full define() function has been called
   */
  bool m_isDefined;

  /** multigrid solver for coarsest level
   */
  LevelNodeSolver m_levelSolver;

  /** whether to write verbose output; default is <tt>true</tt>
   */
  bool m_verbose;

private:
  AMRNodeSolver(const AMRNodeSolver&);
  AMRNodeSolver& operator=(const AMRNodeSolver&);


};
#endif
