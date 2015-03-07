#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodePoissonBC.H
// DTGraves, Mon, July 19, 1999
// petermc, 13 Feb 2001
#ifndef NODEPOISSONBC_H
#define NODEPOISSONBC_H

#include "Box.H"
// #include "NodeFArrayBox.H"
#include "Vector.H"
#include "LoHiSide.H"
#include "NodeBC.H"

#include "UsingNamespace.H"

/// Class to enforce homogeneous Neumann boundary conditions
/**
   FaceNodeBC class implementation to enforce homogeneous Neumann
   boundary conditions \\
   Recall that FaceNodeBC is a class to encapsulate the
   operations  of boundary conditions at a face.
   If the solution is phi and the face normal direction is x,
   the boundary conditions usually used can be expressed as \\
   A*dphi/dx + B*phi = C. \\
   For homogeneous Neumann BCs, B=C=0, A=1. \\
   These functions are only meant to be called by DomainNodeBC
 */
class NeumannNodeBC : public FaceNodeBC
{
public:
  friend class DomainNodeBC;

  ///
  ~NeumannNodeBC();

  ///
  NeumannNodeBC();

  ///
  NeumannNodeBC(int dir, Side::LoHiSide sd);

  ///
  NeumannNodeBC(int dir, Side::LoHiSide sd, const Interval& a_comps);

  ///
  /** Set inhomogeneous boundary condition on this face.
      Default is homogeneous.
   */
  void setInhomogeneous(bool a_inhomogeneous);

protected:

  ///
  /**
     fills neumann factor with 1.0
     and other two fields with zero.
     Thus 1*dphi/dn + 0*phi = 0
  */
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real a_dx,
               const ProblemDomain& a_domain) const;

  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real a_dx,
               const Box& a_domain) const;

  ///
  FaceNodeBC* new_boxBC() const;

  // has inhomogeneous component
  bool m_inhomogeneous;

private:

  NeumannNodeBC(const NeumannNodeBC& neuin){;}

  void operator=(const NeumannNodeBC& neuin){;}

};

/// Class to enforce homogeneous Dirichlet boundary conditions
/**
   FaceNodeBC class implementation to enforce homogeneous Dirichlet
   boundary conditions.\\
   Recall that FaceNodeBC is a class to encapsulate the
   operations  of boundary conditions at a face.
   If the solution is phi and the face normal direction is x,
   the boundary conditions usually used can be expressed as \\
   A*dphi/dx + B*phi = C. \\
   For homogeneous Dirichlet BCs, A=C=0, B=1. \\
   These functions are only meant to be called by DomainNodeBC.
 */
class DirichletNodeBC : public FaceNodeBC
{
public:
  friend class DomainNodeBC;

  ///
  DirichletNodeBC(int dir, Side::LoHiSide sd);

  ///
  DirichletNodeBC(int dir, Side::LoHiSide sd, const Interval& a_comps);

  ///
  ~DirichletNodeBC();

  ///
  DirichletNodeBC();

  ///
  /** Set inhomogeneous boundary condition on this face.
      Default is homogeneous.
   */
  void setInhomogeneous(bool a_inhomogeneous);

protected:

  ///
  /**
     fills Dirichlet factor with 1.0
     and other two fields with zero.
     Thus 0*dphi/dn + 1*phi = 0
  */
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real a_dx,
               const ProblemDomain& a_domain) const;

  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real a_dx,
               const Box& a_domain) const;

  ///
  FaceNodeBC* new_boxBC() const;

  // has inhomogeneous component
  bool m_inhomogeneous;

private:

  DirichletNodeBC(const DirichletNodeBC& dirin){;}

  void operator=(const DirichletNodeBC& dirin){;}

};

#endif
