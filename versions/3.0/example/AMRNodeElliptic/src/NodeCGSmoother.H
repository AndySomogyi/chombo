#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodeCGSmoother.H
// adapted from CGSmoother by DFMartin, Sun, May 5, 2002
// petermc, 5 June 2002

#ifndef NODECGSMOOTHER_H
#define NODECGSMOOTHER_H

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "NodeFArrayBox.H"
#include "ProblemDomain.H"
#include "NodeBaseBottomSmoother.H"

#include "UsingNamespace.H"

/// Implements CG as a bottom smoother for LevelOp-derived classes
/** The NodeCGSmoother class implements the CG algorithm for
    NodeLevelOp-derived classes
*/
class NodeCGSmoother : public NodeBaseBottomSmoother
{

public:

  /**
     \name Constructor and destructor
  */
  /*@{*/

  ///
  /** Default constructor.
   */
  NodeCGSmoother();

  ///
  /** Destructor.
   */
  virtual ~NodeCGSmoother();

  ///
  /** this gets around the "no virtual constructor" rule
   */
  virtual NodeBaseBottomSmoother* new_bottomSmoother() const;

  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  ///
  /** set max number of CG iterations
   */
  void setMaxIter(int a_max_iter);

  ///
  /** set solver tolerance (amount residual must be decreased, etc)
   */
  void setSolverTol(Real a_solverTol);

  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  ///
  /** Perform CG smoothing of L(<i>a_phi</i>) = <i>a_rhs</i>, where
      L is the operator defined in <i>a_levelop_ptr</i>,
      and the preconditioner is
      <i>a_levelop_ptr</i>->levelPreconditioner().
   */
  virtual void doBottomSmooth(/// data at this level, to be smoothed
                              LevelData<NodeFArrayBox>& a_phi,
                              /// right-hand side at this level
                              const LevelData<NodeFArrayBox>& a_rhs,
                              /// pointer to the operator
                              NodeLevelOp* a_levelop_ptr
                              );
  /*@}*/

protected:

  /** max number of iterations
   */
  int m_maxIter;

  /** solver tolerance
   */
  Real m_solverTol;

  /** small parameter to check for division by zero
   */
  Real m_small;

  /** parameter to check for hanging
   */
  Real m_converge_small;
};

#endif
