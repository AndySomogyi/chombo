#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodeBaseBottomSmoother.H
// adapted from BaseBottomSmoother by DFMartin, Sun, May 5, 2002
// petermc, 5 June 2002

#ifndef NODEBASEBOTTOMSMOOTHER_H
#define NODEBASEBOTTOMSMOOTHER_H

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
#include "NodeLevelOp.H"

#include "UsingNamespace.H"

/// Pure base class to encapsulate bottom smoother for NodeLevelOp-derived objects.
/** The NodeBaseBottomSmoother class is a virtual base class for implementing
    different bottomSmoother algorithms for NodeLevelOp-derived classes. It
    takes a NodeLevelOp pointer as an argument during its definition, in order
    to access the problem-dependent (equation-dependent) functionality of
    the specific NodeLevelOp derivative class.
*/

class NodeBaseBottomSmoother
{

public:

  /**
     \name Constructor and destructor
  */
  /*@{*/

  ///
  /** Default constructor.
   */
  NodeBaseBottomSmoother() {};

  ///
  /** Destructor.
   */
  virtual ~NodeBaseBottomSmoother() {};

  ///
  /** This gets around the "no virtual constructor" rule.
   */
  virtual NodeBaseBottomSmoother* new_bottomSmoother() const = 0;

  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  ///
  /** Control output to stdout.
   */
  void setVerbose( bool a_verbose ){ m_verbose = a_verbose ; } ;

  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  ///
  /** Perform smoothing of L(<i>a_phi</i>) = <i>a_rhs</i>, where
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
                              ) = 0;

  /*@}*/

  /** enable output to stdout if <tt>true</tt>
   */
  bool m_verbose;
};

#endif
