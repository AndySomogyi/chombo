#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// AMRNodeSolverAlt.H
// adapted from AMRSolver by DTGraves, Tues, July 6, 1999
// AMRNodeSolver by petermc, 28 Nov 2000
// AMRNodeSolverAlt by petermc, 22 May 2002

#ifndef AMRNODESOLVERALT_H
#define AMRNODESOLVERALT_H

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "SPACE.H"
#include "REAL.H"
#include "ProblemDomain.H"
#include "NodeFArrayBox.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "LevelNodeSolver.H"

#include "UsingNamespace.H"

/// Class which manages grid hierarchy and composite elliptic solution, solving level-by-level
/**
   AMRNodeSolverAlt manages the AMR/multigrid solution to the elliptic
   equation on a set of grids on multiple levels satisfying certain
   proper-nesting conditions.

  It can be used either as a solver, or to apply the AMR / multigrid
  V-cycle as a preconditioner for some other iterative method, such as
  a Krylov method.  */
class AMRNodeSolverAlt
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Creates a AMRNodeSolverAlt whose internal state is undefined.
      Need to call define() function to use any of
      the functionality of the class.
  */
  AMRNodeSolverAlt();

  ///
  /** Creates a fully-defined AMRNodeSolverAlt.   Calls
      define() function with identical arguments.
  */
  AMRNodeSolverAlt(const Vector<DisjointBoxLayout>& a_gridsLevel,
                   const Vector<ProblemDomain>&     a_domainLevel,
                   const Vector<Real>&              a_dxLevel,
                   const Vector<int>&               a_refRatio,
                   int                              a_numLevels,
                   int                              a_lBase,
                   const NodeLevelOp* const         a_opin,
                   int                              a_minLength = 1);

  ///
  /** Creates a fully-defined AMRNodeSolverAlt.   Calls
      define() function with identical arguments.
  */
  AMRNodeSolverAlt(const Vector<DisjointBoxLayout>& a_gridsLevel,
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
  ~AMRNodeSolverAlt();

  ///
  /** Defines AMRNodeSolverAlt's internal state.

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
  /** Defines AMRNodeSolverAlt's internal state.

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
  /** Set tolerance of solvers.  Default is set in LevelNodeSolver.
   */
  void setTolerance(Real a_tolerance);

  ///
  /** Set whether solver does smoothing at bottom level.
      Default is set in LevelNodeSolver.
   */
  void setBottomSmoothing(bool a_doBottomSmooth);

  ///
  /** Set max number of iterations for solvers at every level.
      Default is set in LevelNodeSolver.
   */
  void setMaxIter(int a_maxIter);

  ///
  /** Set whether the solver gives output while solving.  Default is true.
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

  /*@}*/

protected:

  /** Set data to default values.  Not for external use.
   */
  void setDefaultValues();

  /** Delete <i>m_levelSolver</i> Vector.
   */
  void clear();

  /** solver tolerance
   */
  Real m_tolerance;

  /** index of coarsest level on which solution is to be computed.
   */
  int m_lBase;

  /** number of levels allowed
   */
  int m_numLevels;

  /** index of finest level used, equal to <i>m_numLevels</i>-1
   */
  int m_finestLevel;

  /** refinement ratio between levels; <i>m_refRatio</i>[ilev] is
   refinement ratio between levels ilev and ilev+1
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

  /** number of smoother iterations on sweep up
   */
  int m_numSmoothUp;

  /** number of smoother iterations on sweep down
   */
  int m_numSmoothDown;

  /** whether the full define() function has been called
   */
  bool m_isDefined;

  /** multigrid solver for each level
   */
  Vector<LevelNodeSolver*> m_levelSolver;

  /** whether to write verbose output; default is <tt>true</tt>
   */
  bool m_verbose;

private:
  AMRNodeSolverAlt(const AMRNodeSolverAlt&);
  AMRNodeSolverAlt& operator=(const AMRNodeSolverAlt&);


};
#endif
