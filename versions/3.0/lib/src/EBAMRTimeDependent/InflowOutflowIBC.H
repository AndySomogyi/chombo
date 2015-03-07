#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#ifndef _INFLOWOUTFLOWIBC_H_
#define _INFLOWOUTFLOWIBC_H_

#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "EBIBC.H"
#include "EBIBCFactory.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "BaseEBBC.H"
#include "DirichletPoissonEBBC.H"
#include "ExtrapAdvectBC.H"

#include "NamespaceHeader.H"

///
/**
   This class is meant to be  a server for initial boundary conditions
   for all the various stages of BCG INS.
*/
class InflowOutflowIBC: public EBIBC
{
public:
  ///
  /**
     Create poisselle flow in a tube
  */
  InflowOutflowIBC(int a_flowDir, Real a_inflowVel, int a_orderEBBC, bool a_doSlipWalls)
  {
    m_flowDir     = a_flowDir;
    m_inflowVel   = a_inflowVel;
    m_orderEBBC   = a_orderEBBC;
    m_doSlipWalls = a_doSlipWalls;
  }

  ///
  virtual ~InflowOutflowIBC()
  {;}

  ///
  virtual void initializeVelocity(LevelData<EBCellFAB>&    a_velocity,
                                  const DisjointBoxLayout& a_grids,
                                  const EBISLayout&        a_ebisl,
                                  const ProblemDomain&     a_domain,
                                  const RealVect&          a_origin,
                                  const Real&              a_time,
                                  const RealVect&          a_dx) const ;

  ///
  virtual void initializePressure(LevelData<EBCellFAB>&    a_pressure,
                                  const DisjointBoxLayout& a_grids,
                                  const EBISLayout&        a_ebisl,
                                  const ProblemDomain&     a_domain,
                                  const RealVect&          a_origin,
                                  const Real&              a_time,
                                  const RealVect&          a_dx) const ;

  ///
  virtual void initializeScalar ( LevelData<EBCellFAB>&    a_scalar,
                                  const DisjointBoxLayout& a_grids,
                                  const EBISLayout&        a_ebisl,
                                  const ProblemDomain&     a_domain,
                                  const RealVect&          a_origin,
                                  const Real&              a_time,
                                  const RealVect&          a_dx) const ;

  ///
  /**
     Return pressure boundary conditions for domain.
  */
  virtual RefCountedPtr<BaseDomainBCFactory> getPressBC() const ;

  ///
  /**
   */
  virtual RefCountedPtr<BaseDomainBCFactory> getMACVelBC() const ;

  ///
  /**
     The initial conditions of this class are not used.
     The advection class needs boundary conditions
  */
  virtual  RefCountedPtr<EBPhysIBCFactory> getVelAdvectBC(int a_velComp) const ;


  ///
  /**
     The initial conditions of this class are not used.
     The advection class needs boundary conditions
  */
  virtual  RefCountedPtr<EBPhysIBCFactory> getScalarAdvectBC(const int&  a_comp) const ;


  ///
  /**
   */
  virtual RefCountedPtr<BaseDomainBCFactory> getVelBC(int a_icomp) const ;

  ///
  /**
     Return velocity boundary conditions for embedded boundary.
  */
  virtual RefCountedPtr<BaseEBBCFactory> getVelocityEBBC(int a_velComp) const ;

  virtual RefCountedPtr<BaseEBBCFactory> getPressureEBBC() const ;

protected:
  int  m_flowDir;
  Real m_inflowVel;
  int  m_orderEBBC;
  bool m_doSlipWalls;

private:
  /// weak construction is bad, mkay?
  InflowOutflowIBC()
  {
    MayDay::Error("invalid operator");
  }
};

///
/**
 */
class InflowOutflowIBCFactory: public EBIBCFactory
{
public:
  ///
  /**
   */
  InflowOutflowIBCFactory(int a_flowDir, Real a_inflowVel, int a_orderEBBC, bool a_doSlipWalls)
  {
    m_flowDir     = a_flowDir;
    m_inflowVel   = a_inflowVel;
    m_orderEBBC   = a_orderEBBC;
    m_doSlipWalls = a_doSlipWalls;
  }

  ///
  virtual ~InflowOutflowIBCFactory() {};

  ///
  /**
   */
  virtual EBIBC* create()  const
  {
    InflowOutflowIBC* retDerived = new InflowOutflowIBC(m_flowDir, m_inflowVel, m_orderEBBC, m_doSlipWalls);
    EBIBC* retBase =  (EBIBC*)retDerived;
    return retBase;
  }

protected:
  int  m_flowDir;
  Real m_inflowVel;
  int  m_orderEBBC;
  bool m_doSlipWalls;

private:
  /// weak construction is bad but i will allow copy construction and assignment
  InflowOutflowIBCFactory()
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"

#endif
