#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodePoissonOp.H
// adapted from PoissonOp by DTGraves, Weds, July 21, 1999
// petermc, Tue, Dec 12, 2000
#ifndef NODEPOISSONOP_H
#define NODEPOISSONOP_H

#include "NodeLevelOp.H"
#include <cmath>
#include <cstdlib>
#include "BaseFab.H"
#include "NodeCFIVS.H"
#include "NodeBC.H"
#include "NodeQuadCFInterp.H"
#include "NodeBaseBottomSmoother.H"
#include "LoHiSide.H"

#include "UsingNamespace.H"

///  Node-centered Poisson operator
class NodePoissonOp: public NodeLevelOp
/** This class contains stuff to solve div (grad phi) = rhs at one level.
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
  NodePoissonOp();

  ///
  /** Virtual constructor workaround.
      Sets these to be the same as this:
      - m_dombc, boundary condition;
      - m_interpolationDegree, degree of c/f interpolation;
      - m_bottomSmootherPtr, pointer to bottom smoother.
   */
  virtual NodeLevelOp* new_levelop() const;

  ///
  /** Destructor.
   */
  ~NodePoissonOp();

  ///
  /** Full define function.  Makes all coarse/fine
      information and sets internal variables.  The current level
      is taken to be the fine level.
  */
  void define(/// grids at this level
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
              bool                     a_verbose = false);


  ///
  /** Full define function.  Makes all coarse/fine
      information and sets internal variables.  The current level
      is taken to be the fine level.
  */
  void define(/// grids at this level
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
              bool                     a_verbose = false);

  ///
  /** Full define function, based on a coarsening of the finer NodeLevelOp.
      Coarse-level grids are not set because only homogeneous interpolation
      is used.
  */
  void define(/// pointer to NodeLevelOp at next finer level
              const NodeLevelOp* a_opfine,
              /// refinement ratio between <i>a_opfine</i>'s and this level
              int                a_refToFine);

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /** Returns <tt>true</tt> if this object was created with the defining
      constructor or if define() has been called.
  */
  virtual bool isDefined() const;

  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  // /** set max number of bottomSmoother iterations
  // */
  // void setMaxIter(int a_maxIter);

  ///
  /** Control output to stdout.
   */
  virtual void setVerbose( bool a_verbose );

  ///
  /** Set bottom smoother.
   */
  void setBottomSmoother(const NodeBaseBottomSmoother& a_bottomSmoother);

  ///
  /** Set physical boundary conditions.
   */
  void setDomainNodeBC(const DomainNodeBC& a_dombcin);

  ///
  /** Set interpolation degree in NodeQuadCFInterp:
      - 1 for linear in 2-D or bilinear in 3-D;
      - 2 (default) for quadratic in 2-D or biquadratic in 3-D.

      Like setDomainNodeBC, this function must be called
      before define().  You cannot change the interpolation degree later.
   */
  void setInterpolationDegree(int a_interpolationDegree);

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
  void CFInterp(/// the data at this level
                LevelData<NodeFArrayBox>&       a_phi,
                /// the data at the next coarser level
                const LevelData<NodeFArrayBox>& a_phiCoarse,
                /// inhomogeneous physical boundary condition?
                bool                            a_inhomogeneous);

  ///
  /** Homogeneous coarse/fine interpolation operator.
      Fill the nodes of <i>a_phi</i> on the coarse/fine interface with zeroes.
  */
  void homogeneousCFInterp(/// the data at this level
                           LevelData<NodeFArrayBox>& a_phi);

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
  void residualI(/// the residual, at interior nodes
                 LevelData<NodeFArrayBox>&       a_resid,
                 /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                 LevelData<NodeFArrayBox>&       a_phi,
                 /// pointer to data at next coarser level
                 const LevelData<NodeFArrayBox>* a_phiCoarsePtr,
                 /// right-hand side at this level
                 const LevelData<NodeFArrayBox>& a_rhs);

  ///
  /** Evaluate residual of operator,<br>
      <i>a_resid</i> = <i>a_rhs</i> - operator(<i>a_phi</i>),<br>
      with homogeneous coarse/fine boundary conditions,
      homogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after setting to zero
      at the coarse/fine boundary and applying homogeneous
      physical boundary conditions.
  */
  void residualH(/// the residual, at interior nodes
                 LevelData<NodeFArrayBox>&       a_resid,
                 /// the data at this level, set to zero at coarse/fine boundary
                 LevelData<NodeFArrayBox>&       a_phi,
                 /// right-hand side at this level
                 const LevelData<NodeFArrayBox>& a_rhs);

  ///
  /** Evaluate residual of operator,<br>
      <i>a_resid</i> = <i>a_rhs</i> - operator(<i>a_phi</i>),<br>
      with inhomogeneous coarse/fine boundary conditions,
      homogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after interpolating
      from *<i>a_phiCoarsePtr</i> at the coarse/fine boundary and applying
      homogeneous physical boundary conditions.
  */
  void residualIcfHphys(/// the residual, at interior nodes
                        LevelData<NodeFArrayBox>&       a_resid,
                        /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                        LevelData<NodeFArrayBox>&       a_phi,
                        /// pointer to data at next coarser level
                        const LevelData<NodeFArrayBox>* a_phiCoarsePtr,
                        /// right-hand side at this level
                        const LevelData<NodeFArrayBox>& a_rhs);

  ///
  /** Evaluate residual of operator,<br>
      <i>a_resid</i> = <i>a_rhs</i> - operator(<i>a_phi</i>),<br>
      with homogeneous coarse/fine boundary conditions,
      inhomogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after setting to zero
      at the coarse/fine boundary and applying the
      physical boundary conditions.
  */
  void residualHcfIphys(/// the residual, at interior nodes
                        LevelData<NodeFArrayBox>&       a_resid,
                        /// the data at this level, set to zero at coarse/fine boundary
                        LevelData<NodeFArrayBox>&       a_phi,
                        /// right-hand side at this level
                        const LevelData<NodeFArrayBox>& a_rhs);

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
  void applyOpI(/// value of the operator, at interior nodes
                LevelData<NodeFArrayBox>&       a_LofPhi,
                /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                LevelData<NodeFArrayBox>&       a_phi,
                /// pointer to data at next coarser level
                const LevelData<NodeFArrayBox>* a_phiCoarsePtr);

  ///
  /** Evaluate operator,
      homogeneous coarse/fine boundary conditions,
      homogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after setting to zero
      at the coarse/fine boundary and applying homogeneous
      physical boundary conditions.
  */
  void applyOpH(/// value of the operator, at interior nodes
                LevelData<NodeFArrayBox>&       a_LofPhi,
                /// the data at this level, set to zero at coarse/fine boundary
                LevelData<NodeFArrayBox>&       a_phi);

  ///
  /** Evaluate operator,
      inhomogeneous coarse/fine boundary conditions,
      homogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after interpolating
      from *<i>a_phiCoarsePtr</i> at the coarse/fine boundary and applying
      homogeneous physical boundary conditions.
  */
  void applyOpIcfHphys(/// value of the operator, at interior nodes
                       LevelData<NodeFArrayBox>&       a_LofPhi,
                       /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                       LevelData<NodeFArrayBox>&       a_phi,
                       /// pointer to data at next coarser level
                       const LevelData<NodeFArrayBox>* a_phiCoarsePtr);

  ///
  /** Evaluate operator,
      homogeneous coarse/fine boundary conditions,
      inhomogeneous physical boundary conditions.

      The operator is evaluated on <i>a_phi</i>, after setting to zero
      at the coarse/fine boundary and applying the
      physical boundary conditions.
  */
  void applyOpHcfIphys(/// value of the operator, at interior nodes
                       LevelData<NodeFArrayBox>&       a_LofPhi,
                       /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                       LevelData<NodeFArrayBox>&       a_phi);

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
  void gradient(/// value of the gradient, at interior nodes
                LevelData<NodeFArrayBox>&       a_gradPhi,
                /// the data at this level, interpolated from *<i>a_phiCoarsePtr</i> at coarse/fine boundary
                LevelData<NodeFArrayBox>&       a_phi,
                /// pointer to data at next coarser level
                const LevelData<NodeFArrayBox>* a_phiCoarsePtr);

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
  void smooth(/// the data at this level
              LevelData<NodeFArrayBox>&       a_phi,
              /// right-hand side at this level
              const LevelData<NodeFArrayBox>& a_rhs);

  ///
  /** Apply preconditioning to the solution on this level.
      Given <i>a_rhshat</i>, returns <i>a_phihat</i>
      from M(<i>a_phihat</i>) = <i>a_rhshat</i> .
  */
  void levelPreconditioner(/// data at this level
                           LevelData<NodeFArrayBox>&       a_phihat,
                           /// right-hand side at this level
                           const LevelData<NodeFArrayBox>& a_rhshat);

  ///
  /** Do a smoothing at the bottom level.  Default is conjugate gradient.
   */
  void bottomSmoother(/// data at this level
                      LevelData<NodeFArrayBox>&       a_phi,
                      /// right-hand side at this level
                      const LevelData<NodeFArrayBox>& a_rhs);

  /*@}*/

protected:

  /** boundary conditons
   */
  DomainNodeBC m_dombc;

  /** degree of interpolation at c/f interfaces:
      1 for (bi)linear, 2 for (bi)quadratic */
  int m_interpolationDegree;

  // moved to NodeLevelOp
  // physical domain of current level
  // ProblemDomain m_domain;

  /** nodes inside physical domain of current level
   */
  ProblemDomain m_domainInteriorNodes;

  /** CELL-centered grids at the current level
   */
  DisjointBoxLayout m_grids;

  /** Copier object that knows how to perform exchange operations on
      LevelData objects defined on <i>m_grids</i> */
  Copier m_exchangeCopier;

  /** CELL-centered grids at next coarser level, if any
   */
  DisjointBoxLayout m_gridsCoarse;

  /** number of components
   */
  int m_ncomp;

  /** refinement ratio between this and the next coarser level
   */
  int m_refToCoarse;

  /** mesh spacing at this level
   */
  Real m_dx;

  /** mesh spacing at next coarser level
   */
  Real m_dxCoarse;

  /** has full define function been called?
   */
  bool m_isDefined;

  /** has boundary condition been defined?
   */
  bool m_isBCDefined;

  /** is inhomogeneous CF interpolation possible?
   */
  bool m_ihcfiEnabled;

  /** used for figuring out which nodes lie on interface with next coarser level
   */
  LayoutData<NodeCFIVS> m_loCFIVS[SpaceDim];

  /** used for figuring out which nodes lie on interface with next coarser level
   */
  LayoutData<NodeCFIVS> m_hiCFIVS[SpaceDim];

  // weights for coarse/fine interpolation
  // Tuple<FArrayBox*, SpaceDim> m_wtLo, m_wtHi, m_wtC;

  /** coarse-fine quadratic interpolation
   */
  NodeQuadCFInterp m_qcfi;

  /** bottom smoother object
   */
  NodeBaseBottomSmoother* m_bottomSmootherPtr;

  /** interior boundary nodes
   */
  LayoutData< Vector<IntVectSet> > m_IVSV;

  /** whether each section of m_IVSV is a complete box
   */
  LayoutData< BitSet > m_IVSVfull;

  // moved to NodeLevelOp
  // exterior boundary nodes
  // LayoutData< Vector<IntVectSet> > m_IVSVext;

private:
  //internally useful functions

  /**
     Get members related to the coarse/fine interface.
  */
  void setCFIVS();

  /**
     Does GSRB on a level.
  */
  void levelGSRB(LevelData<NodeFArrayBox>& a_phi,
                 const LevelData<NodeFArrayBox>& a_rhs);

  /**
     Applies operator, after coarse/fine and physical
     boundary conditions have been set.
  */
  void applyOpOnly(LevelData<NodeFArrayBox>& a_LofPhi,
                   const LevelData<NodeFArrayBox>& a_phi);

  /**
     Evaluates gradient, after coarse/fine and physical
     boundary conditions have been set.
  */
  void gradientOnly(LevelData<NodeFArrayBox>& a_gradPhi,
                    const LevelData<NodeFArrayBox>& a_phi);

  /**
     Does homogeneous coarse/fine interpolation.
  */
  void homogeneousCFInterp(FArrayBox& a_phiFab,
                           const DataIndex& a_datInd,
                           const int a_idir,
                           const Side::LoHiSide a_hiorlo);

  void clearMemory();

  void setDefaultValues();

  void residualDiff(LevelData<NodeFArrayBox>& a_resid,
                    const LevelData<NodeFArrayBox>& a_rhs);

};

#endif
