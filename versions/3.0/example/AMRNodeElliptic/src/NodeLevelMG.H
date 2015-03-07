#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodeLevelMG.H
// adapted from LevelMG by DTGraves, Fri, July 16, 1999
// petermc, 30 May 2001

#ifndef NODELEVELMG_H
#define NODELEVELMG_H

#include <iostream>
#include <cmath>
#include "SPACE.H"
#include <cstdlib>
#include "REAL.H"
#include "ProblemDomain.H"
#include "DisjointBoxLayout.H"
#include "NodeFArrayBox.H"
#include "LevelData.H"
#include "NodeCoarseAverage.H"
#include "NodeLevelOp.H"

#include "UsingNamespace.H"

/// Multigrid solver on a level
/**
   Multigrid solver on a level.
   This class is to be considered internal
   to LevelNodeSolver and not a part of the Chombo API.
 */
class NodeLevelMG
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.  User must subsequently call define().
   */
  NodeLevelMG();

  ///
  /** Constructor to define NodeLevelMG from the inputs.
      Calls full define() function with identical arguments.
   */
  NodeLevelMG(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_gridsCoarsePtr,
              const ProblemDomain&     a_domain,
              Real                     a_dx,
              int                      a_refToCoarse,
              const NodeLevelOp* const a_opin,
              int                      a_nCoarserLevels);

  ///
  /** Constructor to define NodeLevelMG from the inputs.
      Calls full define() function with identical arguments.
   */
  NodeLevelMG(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_gridsCoarsePtr,
              const Box&               a_domain,
              Real                     a_dx,
              int                      a_refToCoarse,
              const NodeLevelOp* const a_opin,
              int                      a_nCoarserLevels);

  ///
  /** Destructor.
   */
  ~NodeLevelMG();

  ///
  /** Constructor for coarsened version of object.
      Calls define() function with identical arguments.
   */
  NodeLevelMG(NodeLevelMG&       a_L,
              int                a_refToCoarse,
              const NodeLevelOp* a_opin);

  ///
  /** Full define function.
  */
  void define(/// CELL-centered grids on which solution will be computed
              const DisjointBoxLayout& a_grids,
              /// pointer to next coarser AMR level's CELL-centered grids, or NULL if there is no coarser AMR level
              const DisjointBoxLayout* a_gridsCoarsePtr,
              /// CELL-centered physical domain at this level
              const ProblemDomain&     a_domain,
              /// mesh spacing
              Real                     a_dx,
              // refinement ratio to next coarser level (if any)
              int                      a_refToCoarse,
              /// operator L used in solution of L(phi) = rho
              const NodeLevelOp* const a_opin,
              /// number of coarser multigrid levels
              int                      a_nCoarserLevels);

  ///
  /** Full define function.
  */
  void define(/// CELL-centered grids on which solution will be computed
              const DisjointBoxLayout& a_grids,
              /// pointer to next coarser AMR level's CELL-centered grids, or NULL if there is no coarser AMR level
              const DisjointBoxLayout* a_gridsCoarsePtr,
              /// CELL-centered physical domain at this level
              const Box&               a_domain,
              /// mesh spacing
              Real                     a_dx,
              // refinement ratio to next coarser level (if any)
              int                      a_refToCoarse,
              /// operator L used in solution of L(phi) = rho
              const NodeLevelOp* const a_opin,
              /// number of coarser multigrid levels
              int                      a_nCoarserLevels);

  ///
  /** Function to define coarsened version of this NodeLevelMG.
   */
  void define(/// original object
              const NodeLevelMG& a_L,
              /// refinement ratio between <i>a_L</i> and this object
              int                a_refToCoarse,
              /// operator L used in solution of L(phi) = rho at same level as <i>a_L</i>
              const NodeLevelOp* a_opin);

  // void clear();

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /** Returns true if this object was created with the defining
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
   */
  void setnumSmoothUp(int a_numSmoothUp);

  ///
  /** Set number of smoothings on way down V-cycle; Default is 4.
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
  /** Invoke relaxation step. Default is pure multigrid V-cycle, suitable
      for use in multilevel solver application; otherwise, use approximate
      solver such as CG at bottom level. It is assumed that the problem has
      already been put in residual-correction form. In particular, only
      the homogeneous form of the physical and coarse/fine boundary
      conditions need be invoked.

      Modified are <i>m_resid</i>, <i>m_crseResid</i>, and <i>m_crseCorr</i>.
  */
  void mgRelax(/// solution on this level
               LevelData<NodeFArrayBox>&       a_phi,
               /// right-hand side on this level
               const LevelData<NodeFArrayBox>& a_rhs,
               /// whether to use bottom solver at coarsest level
               bool                            a_bottomsolveflag);

  /*@}*/

  /** this is a dangerous access function that should not generally be used.
   */
  NodeLevelOp* levelOpPtr();

  /** this is another access function that is kinda bad
   */
  NodeLevelMG* lCoarsePtr();

protected:

  /** Set default values.  Does not deal with memory.
   */
  void setDefaultValues();

  /** Delete internal memory and set pointers to NULL.
   */
  void clearMemory();

  /** has this NodeLevelMG been defined yet?
   */
  bool m_isDefined;

  /** flag for verbose output
   */
  bool m_verbose;

  //these are owned by NodeLevelMG

  /** residual at this level
   */
  LevelData<NodeFArrayBox> m_resid;

  /** residual at next coarser level
   */
  LevelData<NodeFArrayBox> m_crseResid;

  /** correction at next coarser level
   */
  LevelData<NodeFArrayBox> m_crseCorr;

  /** grids at this level
   */
  DisjointBoxLayout m_grids;

  /** coarsened version of CELL-centered grids at this level
   */
  DisjointBoxLayout m_coarsenedGrids;

  /** pointer to CELL-centered grids at next coarser AMR level, or NULL if no coarser level
   */
  const DisjointBoxLayout* m_gridsCoarsePtr;

  /** mesh spacing at this level
   */
  Real m_dx;

  /** CELL-centered physical domain at this level
   */
  ProblemDomain m_domain;

  /** number of coarser multigrid levels
   */
  int m_nCoarserLevels;

  /** pointer to coarser level
   */
  NodeLevelMG* m_lCoarsePtr;

  /** the operator
   */
  NodeLevelOp* m_levelopPtr;

  /** number of CELLs in physical domain
   */
  int m_numBottomPoints;

  /** number of smoother iterations on sweep up
   */
  int m_numSmoothUp;

  /** number of smoother iterations on sweep down
   */
  int m_numSmoothDown;

  /** refinement ratio to next coarser multigrid level
   */
  int m_refToCoarsened;

  /** operator for averaging from this level to next coarser multigrid level
   */
  NodeCoarseAverage m_averageOp;

  /** reference box used in interpolation from next coarser multigrid level
   */
  Box m_boxRef;

  /** weights of corners in interpolation from next coarser multigrid level
   */
  FArrayBox m_weights;

private:

  /** Correct <i>a_phi</i> on intersection with <i>m_crseCorr</i>.
      Should only be called internally because this
      is not written for general LevelData<NodeFArrayBox>es */
  void crseCorrect(LevelData<NodeFArrayBox>& a_phi);

  void operator=(const NodeLevelMG& a_levmgin){};

  NodeLevelMG(const NodeLevelMG& a_levmgin) {};
};

#endif
