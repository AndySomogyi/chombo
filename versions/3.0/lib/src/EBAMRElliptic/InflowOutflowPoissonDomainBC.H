#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _INFLOWOUTFLOWPOISSONDOMAINBC_H_
#define _INFLOWOUTFLOWPOISSONDOMAINBC_H_

#include "RefCountedPtr.H"
#include "BaseDomainBC.H"

#include "NamespaceHeader.H"

///
/**
 */
class InflowOutflowPoissonDomainBC: public BaseDomainBC
{
public:

  ///
  /**
   */
  virtual ~InflowOutflowPoissonDomainBC()
  {;}

  ///
  /**
     Elliptic solver flux.
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
     This is used by the projections to get
     velocity at domain faces.  Sets velocity to a InflowOutflow  value.
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

  ///
  InflowOutflowPoissonDomainBC(const int& a_flowDir, const Real& a_inflowVel)
  {
    m_flowDir = a_flowDir;
    m_inflowVel = a_inflowVel;
  }

  ///
  /**
     Return the a face centered velocity that makes a_vof div-free
  */
  static Real getDivFreeFace(const FaceIndex&      a_face,
                             const VolIndex&       a_vof,
                             const int&            a_dir,
                             const Side::LoHiSide& a_side,
                             const EBFluxFAB&      a_vel,
                             const RealVect&       a_dx,
                             const EBISBox&        a_ebisBox);

private:

  ///weak construction forbidden to keep things simple
  InflowOutflowPoissonDomainBC(){;};
  int m_flowDir;
  Real m_inflowVel;

};

class InflowOutflowPoissonDomainBCFactory: public BaseDomainBCFactory
{
public:


  ///
  /**
   */
  InflowOutflowPoissonDomainBCFactory(const int& a_flowDir, const Real& a_inflowVel)
  {
    m_flowDir   = a_flowDir;
    m_inflowVel = a_inflowVel;
  }


  ///
  /**
   */
  virtual ~InflowOutflowPoissonDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual InflowOutflowPoissonDomainBC* create(const ProblemDomain& a_domain,
                                               const EBISLayout&    a_layout,
                                               const RealVect&      a_dx)
  {
    InflowOutflowPoissonDomainBC* newBC = new InflowOutflowPoissonDomainBC(m_flowDir, m_inflowVel);

    return newBC;
  }


private:
  InflowOutflowPoissonDomainBCFactory(){;}

  int m_flowDir;
  Real m_inflowVel;

};



///
/**
 */
class InflowOutflowHelmholtzDomainBC: public BaseDomainBC
{
public:

  ///
  /**
   */
  virtual ~InflowOutflowHelmholtzDomainBC()
  {;}

  ///
  /**
     Elliptic solver flux.
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
     This is used by the projections to get
     velocity at domain faces.  Sets velocity to a InflowOutflow  value.
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


  ///
  InflowOutflowHelmholtzDomainBC(const int& a_flowDir, const Real& a_inflowVel, const bool& a_doSlipWalls)
  {
    m_flowDir   = a_flowDir;
    m_inflowVel = a_inflowVel;
    m_doSlipWalls = a_doSlipWalls;
  }

private:
  //weak construction forbidden
  InflowOutflowHelmholtzDomainBC(){;};

  int m_flowDir;
  bool m_doSlipWalls;
  Real m_inflowVel;
};

///
/**
 **/
class InflowOutflowHelmholtzDomainBCFactory: public BaseDomainBCFactory
{
public:

  ///
  InflowOutflowHelmholtzDomainBCFactory(const int& a_flowDir, const Real& a_inflowVel, const bool& a_doSlipWalls)
  {
    m_flowDir   = a_flowDir;
    m_inflowVel = a_inflowVel;
    m_doSlipWalls = a_doSlipWalls;
  }


  ///
  /**
   */
  virtual ~InflowOutflowHelmholtzDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual InflowOutflowHelmholtzDomainBC* create(const ProblemDomain& a_domain,
                                                 const EBISLayout&    a_layout,
                                                 const RealVect&      a_dx)
  {
    InflowOutflowHelmholtzDomainBC* newBC = new InflowOutflowHelmholtzDomainBC(m_flowDir, m_inflowVel, m_doSlipWalls);

    return newBC;
  }

private:
  //weak construction
  InflowOutflowHelmholtzDomainBCFactory(){;}

  bool m_doSlipWalls;
  int m_flowDir;
  Real m_inflowVel;
};

#include "NamespaceFooter.H"

#endif
