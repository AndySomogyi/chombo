#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodeNoOpSmoother.H
// adapted from NoOpSmoother by DFMartin, Sun, May 5, 2002
// petermc, 5 June 2002

#ifndef NODENOOPSMOOTHER_H
#define NODENOOPSMOOTHER_H

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

/// Implements NoOp as a bottom smoother for LevelOp-derived classes
/** The NodeNoOpSmoother class implements the NoOp algorithm for
    NodeLevelOp-derived classes
*/
class NodeNoOpSmoother : public NodeBaseBottomSmoother
{

public:

  /**
     \name Constructor and destructor
  */
  /*@{*/

  ///
  /** Default constructor.
   */
  NodeNoOpSmoother();

  ///
  /** Destructor.
   */
  virtual ~NodeNoOpSmoother();

  ///
  /** this gets around the "no virtual constructor" rule
   */
  virtual NodeBaseBottomSmoother* new_bottomSmoother() const;

  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  /** Do nothing.
   */
  virtual void doBottomSmooth(LevelData<NodeFArrayBox>& a_phi,
                              const LevelData<NodeFArrayBox>& a_rhs,
                              NodeLevelOp* a_levelop_ptr);
  /*@}*/
};

#endif
