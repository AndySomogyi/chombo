#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodeLevelOp.H
// adapted from LevelOp by DTGraves, Mon, July 26, 1999
// petermc, Tue, Nov 28, 2000
#ifndef NODELEVELOP_H
#define NODELEVELOP_H

#include <iostream>
#include <math.h>
#include "SPACE.H"
#include <stdlib.h>
#include "REAL.H"
#include "ProblemDomain.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "NodeFArrayBox.H"
#include "IntVectSet.H"

#include "UsingNamespace.H"

/// Pure base class to encapsulate level operations API for node-centered elliptic solvers.

class NodeLevelOp
/**
    NodeLevelOp is a pure base class to encapsulate level operations API
    for node-centered elliptic solvers.  All elliptic operators used by
    AMRNodeSolver inherit the NodeLevelOp interface.
*/
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.  User must subsequently call define().
   */
  NodeLevelOp() {};

  ///
  /** This gets around the "no virtual constructor" rule.
  */
  // "= 0" is a C++ism meaning that this is a _pure_ virtual function:
  // every derived class must override it.
  virtual NodeLevelOp* new_levelop() const = 0;

  ///
  /** Destructor.
   */
  virtual ~NodeLevelOp(){};

  ///
  /** Full define function.  Makes all coarse/fine
      information and sets internal variables.  The current level
      is taken to be the fine level.
  */
  virtual void define(/// grids at this level
                      const DisjointBoxLayout& a_grids,
                      /// pointer to grids at next coarser AMR level, or NULL if none
                      const DisjointBoxLayout* a_gridsCoarsePtr,
                      /// mesh spacing at this level
                      Real                     a_dx,
                      /// refinement ratio between this and next coarser level
                      int                      a_refToCoarse,
                      /// physical domain at this level
                      const ProblemDomain&     a_domain,
                      /// only homogeneous coarse/fine interpolation?
                      bool                     a_homogeneousOnly = false,
                      /// number of components
                      int                      a_ncomp = 1,
                      /// verbose output?
                      bool                     a_verbose = false
                      ) = 0;

  ///
  /** Full define function.  Makes all coarse/fine
      information and sets internal variables.  The current level
      is taken to be the fine level.
  */
  virtual void define(/// grids at this level
                      const DisjointBoxLayout& a_grids,
                      /// pointer to grids at next coarser AMR level, or NULL if none
                      const DisjointBoxLayout* a_gridsCoarsePtr,
                      /// mesh spacing at this level
                      Real                     a_dx,
                      /// refinement ratio between this and next coarser level
                      int                      a_refToCoarse,
                      /// physical domain at this level
                      const Box&               a_domain,
                      /// only homogeneous coarse/fine interpolation?
                      bool                     a_homogeneousOnly = false,
                      /// number of components
                      int                      a_ncomp = 1,
                      /// verbose output?
                      bool                     a_verbose = false
                      ) = 0;

  ///
  /** Full define function, based on a coarsening of the finer NodeLevelOp.
      Coarse-level grids are not set because only homogeneous interpolation
      is used.
  */
  virtual void define(/// pointer to NodeLevelOp at next finer level
                      const NodeLevelOp* a_opfine,
                      /// refinement ratio between <i>a_opfine</i>'s and this level
                      int                a_refToFine
                      ) = 0;

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /** Returns <tt>true</tt> if this object was created with the defining
      constructor or if define() has been called.
  */
  virtual bool isDefined() const = 0;

  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  /** Control output to stdout.
   */
  virtual void setVerbose( bool a_verbose ){ m_verbose = a_verbose ; } ;

  /*@}*/

  /**
     \name Coarse/fine interpolation functions
  */
  /*@{*/

  ///
  /** Coarse / Fine (inhomogeneous) interpolation operator.
      Fill the nodes of <i>a_phi</i> on the coarse/fine interface with
      interpolated data from <i>a_phiCoarse</i>.
  */
  virtual void CFInterp(/// the data at this level
                        LevelData<NodeFArrayBox>&       a_phi,
                        /// the data at the next coarser level
                        const LevelData<NodeFArrayBox>& a_phiCoarse,
                        /// inhomogeneous physical boundary condition?
                        bool                            a_inhomogeneous
                        ) = 0;

  ///
  /** Homogeneous coarse/fine interpolation operator.
      Fill the nodes of <i>a_phi</i> on the coarse/fine interface with zeroes.
  */
  virtual void homogeneousCFInterp(/// the data at this level
                                   LevelData<NodeFArrayBox>& a_phi
                                   ) = 0;

  /*@}*/

  /**
     \name Evaluation of residual
  */
  /*@{*/

  ///
  /** Evaluate residual of operator,<br>
      <i>a_resid</i> = <i>a_rhs</i> - operator(<i>a_phi</i>),<br>
      with inhomogeneous coarse/fine boundary conditions,
      inhomogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after interpolating
      from *<i>a_phiCoarsePtr</i> at the coarse/fine boundary and applying
      the physical boundary conditions.
  */
  virtual void residualI(/// the residual, at interior nodes
                         LevelData<NodeFArrayBox>&       a_resid,
                         /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                         LevelData<NodeFArrayBox>&       a_phi,
                         /// pointer to data at next coarser level
                         const LevelData<NodeFArrayBox>* a_phiCoarsePtr,
                         /// right-hand side at this level
                         const LevelData<NodeFArrayBox>& a_rhs) = 0;


  ///
  /** Evaluate residual of operator,<br>
      <i>a_resid</i> = <i>a_rhs</i> - operator(<i>a_phi</i>),<br>
      with homogeneous coarse/fine boundary conditions,
      homogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after setting to zero
      at the coarse/fine boundary and applying homogeneous
      physical boundary conditions.
  */
  virtual void residualH(/// the residual, at interior nodes
                         LevelData<NodeFArrayBox>&       a_resid,
                         /// the data at this level, set to zero at coarse/fine boundary
                         LevelData<NodeFArrayBox>&       a_phi,
                         /// right-hand side at this level
                         const LevelData<NodeFArrayBox>& a_rhs) = 0;


  ///
  /** Evaluate residual of operator,<br>
      <i>a_resid</i> = <i>a_rhs</i> - operator(<i>a_phi</i>),<br>
      with inhomogeneous coarse/fine boundary conditions,
      homogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after interpolating
      from *<i>a_phiCoarsePtr</i> at the coarse/fine boundary and applying
      homogeneous physical boundary conditions.
  */
  virtual void residualIcfHphys(/// the residual, at interior nodes
                                LevelData<NodeFArrayBox>&       a_resid,
                                /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                                LevelData<NodeFArrayBox>&       a_phi,
                                /// pointer to data at next coarser level
                                const LevelData<NodeFArrayBox>* a_phiCoarsePtr,
                                /// right-hand side at this level
                                const LevelData<NodeFArrayBox>& a_rhs) = 0;


  ///
  /** Evaluate residual of operator,<br>
      <i>a_resid</i> = <i>a_rhs</i> - operator(<i>a_phi</i>),<br>
      with homogeneous coarse/fine boundary conditions,
      inhomogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after setting to zero
      at the coarse/fine boundary and applying the
      physical boundary conditions.
  */
  virtual void residualHcfIphys(/// the residual, at interior nodes
                                LevelData<NodeFArrayBox>&       a_resid,
                                /// the data at this level, set to zero at coarse/fine boundary
                                LevelData<NodeFArrayBox>&       a_phi,
                                /// right-hand side at this level
                                const LevelData<NodeFArrayBox>& a_rhs) = 0;

  /*@}*/

  /**
     \name Evaluation of operator
  */
  /*@{*/

  ///
  /** Evaluate operator,
      inhomogeneous coarse/fine boundary conditions,
      inhomogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i> after interpolating
      from *<i>a_phiCoarsePtr</i> at the coarse/fine boundary and applying
      the physical boundary conditions.
  */
  virtual void applyOpI(/// value of the operator, at interior nodes
                        LevelData<NodeFArrayBox>&       a_LofPhi,
                        /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                        LevelData<NodeFArrayBox>&       a_phi,
                        /// pointer to data at next coarser level
                        const LevelData<NodeFArrayBox>* a_phiCoarsePtr
                        ) = 0;

  ///
  /** Evaluate operator,
      homogeneous coarse/fine boundary conditions,
      homogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after setting to zero
      at the coarse/fine boundary and applying homogeneous
      physical boundary conditions.
  */
  virtual void applyOpH(/// value of the operator, at interior nodes
                        LevelData<NodeFArrayBox>&       a_LofPhi,
                        /// the data at this level, set to zero at coarse/fine boundary
                        LevelData<NodeFArrayBox>&       a_phi
                        ) = 0;

  ///
  /** Evaluate operator,
      inhomogeneous coarse/fine boundary conditions,
      homogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after interpolating
      from *<i>a_phiCoarsePtr</i> at the coarse/fine boundary and applying
      homogeneous physical boundary conditions.
  */
  virtual void applyOpIcfHphys(/// value of the operator, at interior nodes
                               LevelData<NodeFArrayBox>&       a_LofPhi,
                               /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                               LevelData<NodeFArrayBox>&       a_phi,
                               /// pointer to data at next coarser level
                               const LevelData<NodeFArrayBox>* a_phiCoarsePtr
                               ) = 0;

  ///
  /** Evaluate operator,
      homogeneous coarse/fine boundary conditions,
      inhomogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after setting to zero
      at the coarse/fine boundary and applying the
      physical boundary conditions.
  */
  virtual void applyOpHcfIphys(/// value of the operator, at interior nodes
                               LevelData<NodeFArrayBox>&       a_LofPhi,
                               /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                               LevelData<NodeFArrayBox>&       a_phi
                               ) = 0;

  /*@}*/

  /**
     \name Evaluation of gradient
  */
  /*@{*/

  ///
  /** Evaluate gradient,
      inhomogeneous coarse/fine boundary conditions,
      inhomogeneous physical boundary conditions.

      The gradient is evaluated on <i>a_phi</i> after interpolating
      from *<i>a_phiCoarsePtr</i> at the coarse/fine boundary and applying
      the physical boundary conditions.
  */
  virtual void gradient(/// value of the gradient, at interior nodes
                        LevelData<NodeFArrayBox>&       a_gradPhi,
                        /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                        LevelData<NodeFArrayBox>&       a_phi,
                        /// pointer to data at next coarser level
                        const LevelData<NodeFArrayBox>* a_phiCoarsePtr
                        ) = 0;

  /*@}*/

  /**
     \name Smoothing functions
  */
  /*@{*/

  ///
  /** Smoother.
      Assumes that problem has already been put in
      residual-correction form,
      so that coarse/fine boundary conditions are homogeneous.
  */
  virtual void smooth(/// the data at this level
                      LevelData<NodeFArrayBox>&       a_phi,
                      /// right-hand side at this level
                      const LevelData<NodeFArrayBox>& a_rhs
                      ) = 0;

  ///
  /** Apply preconditioning to the solution on this level.
      Given <i>a_rhshat</i>, returns <i>a_phihat</i>
      from M(<i>a_phihat</i>) = <i>a_rhshat</i> .
  */
  virtual void levelPreconditioner(/// data at this level
                                   LevelData<NodeFArrayBox>&       a_phihat,
                                   /// right-hand side at this level
                                   const LevelData<NodeFArrayBox>& a_rhshat
                                   ) = 0;

  ///
  /** Smoother at bottom level.
  */
  virtual void bottomSmoother(/// data at this level
                              LevelData<NodeFArrayBox>&       a_phi,
                              /// right-hand side at this level
                              const LevelData<NodeFArrayBox>& a_rhs) = 0;

  /*@}*/

  // these have to be public because smoother needs them for norms

  /** CELL-centered physical domain of current level
   */
  ProblemDomain m_domain;

  /** exterior boundary nodes
   */
  LayoutData< Vector<IntVectSet> > m_IVSVext;

  /** enable output to stdout if <tt>true</tt>
   */
  bool m_verbose;
};

#endif
