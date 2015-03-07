#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _DIRICHLETPOISSONDOMAINBC_H_
#define _DIRICHLETPOISSONDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"
#include "BaseBCValue.H"
#include "NamespaceHeader.H"

// #define BOXSIZE_HO 2
///
/**
 */
class DirichletPoissonDomainBC: public BaseDomainBC
{
public:
  ///
  /**
   */
  DirichletPoissonDomainBC();

  ///
  /**
   */
  virtual ~DirichletPoissonDomainBC();

  ///
  /**
   */
  virtual void setValue(Real a_value);

  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_func);

  ///
  /**
     0 Neumann
     1 Dirichlet
     2 periodic
     Anything else is invalid
  */
  virtual int whichBC(int                  a_idir,
                      Side::LoHiSide       a_side)
  {
    return 1;
  };

  ///
  /**
   */
  virtual void setEBOrder(int a_ebOrder);

  ///
  /**
   */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  /**
     Elliptic solver flux.
   */
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);
  ///
  /**
     This is used by the projections to get
     grad(phi) at domain faces.
   */
  virtual void getFaceGradPhi(Real&                 a_faceFlux,
                              const FaceIndex&      a_face,
                              const int&            a_comp,
                              const EBCellFAB&      a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_dx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useAreaFrac,
                              const RealVect&       a_centroid,
                              const bool&           a_useHomogeneous);

  ///
  /**
     This is used by the operator to get
     grad(phi) at domain faces.
   */
  virtual void getFaceFluxGradPhi(Real&                 a_faceFlux,
                                  const FaceIndex&      a_face,
                                  const int&            a_comp,
                                  const EBCellFAB&      a_phi,
                                  const RealVect&       a_probLo,
                                  const RealVect&       a_dx,
                                  const int&            a_idir,
                                  const Side::LoHiSide& a_side,
                                  const DataIndex&      a_dit,
                                  const Real&           a_time,
                                  const bool&           a_useAreaFrac,
                                  const RealVect&       a_centroid,
                                  const bool&           a_useHomogeneous);

  ///
  /**
     This is used by the projections to get
     velocity at domain faces.  Sets velocity to a dirichlet  value.
   */
  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_face,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side,
                          const bool&           a_doDivFreeOutflow);


private:
  bool m_onlyHomogeneous;
  bool m_isFunctional;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_func;

  int m_ebOrder;
};

class DirichletPoissonDomainBCFactory: public BaseDomainBCFactory
{
public:
  ///
  /**
   */
  DirichletPoissonDomainBCFactory();

  ///
  /**
   */
  virtual ~DirichletPoissonDomainBCFactory();

  ///
  /**
   */
  virtual DirichletPoissonDomainBC* create(const ProblemDomain& a_domain,
                                           const EBISLayout&    a_layout,
                                           const RealVect&      a_dx);

  ///
  /**
   */
  virtual void setValue(Real a_value);

  ///
  /**
   */
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_func);

  ///
  /**
   */
  virtual void setEBOrder(int a_ebOrder);

private:
  bool m_onlyHomogeneous;
  bool m_isFunctional;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_func;

  int m_ebOrder;
};
#include "NamespaceFooter.H"
#endif
