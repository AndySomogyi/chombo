#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// LevelNodeSolver.H
// adapted from LevelSolver by DTGraves, Fri, July 23, 1999
// petermc, 12 Apr 2001

#ifndef LEVELNODESOLVER_H
#define LEVELNODESOLVER_H

#include <cstdlib>
#include <iostream>
#include "SPACE.H"
#include <cmath>
#include "REAL.H"
#include "IntVectSet.H"
#include "ProblemDomain.H"
#include "NodeFArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "NodeLevelMG.H"
#include "NodeLevelOp.H"

#include "UsingNamespace.H"

/// Solves elliptic equations on a level.
class LevelNodeSolver
/** LevelNodeSolver solves elliptic equations on a level using multigrid.
*/
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves LevelNodeSolver undefined.
  */
  LevelNodeSolver();

  ///
  /**
     Full constructor.
     Calls full define() function with identical arguments.
  */
  LevelNodeSolver(const DisjointBoxLayout& a_grids,
                  const DisjointBoxLayout* a_gridsCoarsePtr,
                  const ProblemDomain&     a_domain,
                  Real                     a_dx,
                  int                      a_refToCoarse,
                  const NodeLevelOp* const a_opin,
                  int                      a_minLength = 1);

  ///
  /**
     Full constructor.
     Calls full define() function with identical arguments.
  */
  LevelNodeSolver(const DisjointBoxLayout& a_grids,
                  const DisjointBoxLayout* a_gridsCoarsePtr,
                  const Box&               a_domain,
                  Real                     a_dx,
                  int                      a_refToCoarse,
                  const NodeLevelOp* const a_opin,
                  int                      a_minLength = 1);

  ///
  /** Destructor.
   */
  virtual ~LevelNodeSolver();

  ///
  /**
     Full define function.
  */
  virtual void define(/// grids on which solution will be computed
                      const DisjointBoxLayout& a_grids,
                      /// pointer to next coarser AMR level's grids, or NULL if there is no coarser level
                      const DisjointBoxLayout* a_gridsCoarsePtr,
                      /// domain of solution level
                      const ProblemDomain&     a_domain,
                      /// grid spacing at solution level
                      Real                     a_dx,
                      /// refinement ratio to next coarser level (if any)
                      int                      a_refToCoarse,
                      /// operator L used in solution of L(phi)=rho
                      const NodeLevelOp* const a_opin,
                      /// minimum side length of maximally coarsened box, or 0 if no coarsening; set to 1 for max coarsening
                      int                      a_minLength = 1);

  ///
  /**
     Full define function.
  */
  virtual void define(/// grids on which solution will be computed
                      const DisjointBoxLayout& a_grids,
                      /// pointer to next coarser AMR level's grids, or NULL if there is no coarser level
                      const DisjointBoxLayout* a_gridsCoarsePtr,
                      /// domain of solution level
                      const Box&               a_domain,
                      /// grid spacing at solution level
                      Real                     a_dx,
                      /// refinement ratio to next coarser level (if any)
                      int                      a_refToCoarse,
                      /// operator L used in solution of L(phi)=rho
                      const NodeLevelOp* const a_opin,
                      /// minimum side length of maximally coarsened box, or 0 if no coarsening; set to 1 for max coarsening
                      int                      a_minLength = 1);

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


  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  ///
  /** Set number of smoothings on way up V-cycle; Default is 4.
      Passed to NodeLevelMG.
   */
  void setnumSmoothUp(int a_numSmoothUp);

  ///
  /** Set number of smoothings on way down V-cycle; Default is 4.
      Passed to NodeLevelMG.
   */
  void setnumSmoothDown(int a_numSmoothDown);

  ///
  /** Sets whether smoothing is done at bottom level.
      Default is <tt>true</tt>.
   */
  void setBottomSmoothing(bool a_bottomSolveFlag);

  ///
  /** Sets solution tolerance.  Default is 1.0e-10.
   */
  void setTolerance(Real a_tolerance);

  ///
  /**
     Set "operator tolerance" of iterative solution.  Iteration will
     stop if
     (new_residual/old_residual) > 1-<i>a_operatorTolerance</i>
     (and at least <i>m_minIter</i> iterations have been performed).
     Default is 1.0e-4.
  **/
  void setOperatorTolerance(Real a_operatorTolerance);

  ///
  /** Sets whether the solver gives output while solving.
      Default is <tt>false</tt>.
   */
  void setVerbose(bool a_verbose);

  ///
  /** Sets max number of multigrid V-cycles to use in solution;
      Default is 33.
   */
  void setMaxIter(int a_maxIter);

  ///
  /**
     Set min number of iterations.  Only relevant when residual is
     not decreasing fast enough to satisfy the "operator tolerance".
     Default is 4.
  */
  void setMinIter(int a_minIter);

  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  /** Delete internal memory and set pointers to NULL.
   */
  void clearMemory();

  /** Set default values.  Does not deal with memory.
   */
  void setDefaultValues();

  ///
  /** Does a level solve on this level using multigrid
      and homogeneous bcs at coarse/fine interfaces.
  */
  void levelSolveH(/// input initial guess at solution, output solution
                   LevelData<NodeFArrayBox>&       a_phi,
                   /// right-hand side
                   const LevelData<NodeFArrayBox>& a_rhs,
                   /// whether to initialize guess of phi to zero
                   bool                            a_initializePhiToZero = true);

  ///
  /** Does a level solve on this level using multigrid
      and inhomogeneous bcs at coarse/fine interfaces.
  */
  void levelSolve(/// input initial guess at solution, output solution
                  LevelData<NodeFArrayBox>&       a_phi,
                  /// pointer to solution at next coarser level, or NULL if there is no coarser level
                  const LevelData<NodeFArrayBox>* a_phiCoarse,
                  /// right-hand side
                  const LevelData<NodeFArrayBox>& a_rhs,
                  /// whether to initialize guess of phi to zero
                  bool                            a_initializePhiToZero = true);

  /*@}*/

protected:

  /** level solver object to relax on this level
   */
  NodeLevelMG m_levelMG;

  /** pointer to operator
   */
  NodeLevelOp* m_levelOpPtr;

  /** has this LevelNodeSolver been defined yet?
   */
  bool m_isDefined;

  /** CELL-centered grids at this level
   */
  DisjointBoxLayout m_grids;

  /** CELL-centered physical domain at this level
   */
  ProblemDomain m_domain;

  /** residual
   */
  LevelData<NodeFArrayBox> m_resid;

  /** useful scratch space
   */
  LevelData<NodeFArrayBox> m_scratch;

  /** correction
   */
  LevelData<NodeFArrayBox> m_corr;

  /** refinement ratio to next coarser level
   */
  int m_refToCoarse;

  /** mesh spacing at this level
   */
  Real m_dx;

  /** maximum number of iterations
   */
  int m_maxIter;

  /** minimum number of iterations
   */
  int m_minIter;

  /** solver tolerance
   */
  Real m_tolerance;

  /** operator tolerance:  iteration will
      stop if (new_residual/old_residual) > 1-<i>m_operatorTolerance</i>
   */
  Real m_operatorTolerance;

  /** whether to do smoothing at bottom level
   */
  bool m_bottomSolveFlag;

  /** whether to write verbose output; default is <tt>true</tt>
   */
  bool m_verbose;

  /** interior boundary nodes of grids at this level
   */
  LayoutData< Vector<IntVectSet> > m_IVSV;

  /** exterior boundary nodes of grids at this level
   */
  LayoutData< Vector<IntVectSet> > m_IVSVext;

private:

  int countCoarserLevels(int a_minLength);

  void operator = (const LevelNodeSolver&){;}
  LevelNodeSolver(const LevelNodeSolver&){;}

};

#endif
