#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _INFLOWOUTFLOWADVECTBC_H_
#define _INFLOWOUTFLOWADVECTBC_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "EBPhysIBC.H"
#include "EBPhysIBCFactory.H"

#include "NamespaceHeader.H"

///for velocitie boundary conditions
/**
   sets = inflow vel at inflow face.
   sets = zero at noflow faces.
   sets = extrapval at outflow face
*/
class InflowOutflowAdvectBC : public EBPhysIBC
{
public:
  virtual ~InflowOutflowAdvectBC()
  {;}

  ///
  InflowOutflowAdvectBC(int a_flowDir, Real a_inflowVel, int a_velComp)
  {
    m_flowDir   = a_flowDir;
    m_inflowVel = a_inflowVel;
    m_velComp   = a_velComp;
    m_isDefined = false;
  }

  ///
  void define(const ProblemDomain&  a_domain,
              const RealVect&       a_dx)
  {
    m_domain    = a_domain;
    m_dx        = a_dx;
    m_isDefined = true;
  }

  ///  For every box in this level, this function is called.
  void fluxBC(EBFluxFAB&            a_flux,
              const EBCellFAB&      a_Wcenter,
              const EBCellFAB&      a_Wextrap,
              const Side::LoHiSide& a_sd,
              const Real&           a_time,
              const EBISBox&        a_ebisBox,
              const DataIndex&      a_dit,
              const Box&            a_box,
              const Box&            a_faceBox,
              const int&            a_dir);

  /// Initialize---no op here
  void initialize(LevelData<EBCellFAB>& a_conState,
                  const EBISLayout& a_ebisl) const
  {
    MayDay::Error("should not be called");
  }

  /// Set boundary slopes --no op here
  void setBndrySlopes(EBCellFAB&       a_deltaPrim,
                      const EBCellFAB& a_primState,
                      const EBISBox&   a_ebisBox,
                      const Box&       a_box,
                      const int&       a_dir)
  {;}


protected:
  bool           m_isDefined;
  int            m_velComp;
  ProblemDomain  m_domain;
  RealVect       m_dx;
  int m_flowDir;
  Real m_inflowVel;

private:
  //weak construction is bad but I will allow copy construction and assignment
  //because there are no pointers in the internal data
  InflowOutflowAdvectBC()
  {
    MayDay::Error("invalid operator");
  }

};


///
/**
 */
class InflowOutflowAdvectBCFactory: public EBPhysIBCFactory
{
public:

  ///
  ~InflowOutflowAdvectBCFactory()
  {;}

  ///
  InflowOutflowAdvectBCFactory(int a_flowDir, Real a_inflowVel, int a_velComp)
  {
    m_velComp = a_velComp;
    m_flowDir = a_flowDir;
    m_inflowVel = a_inflowVel;
  }

  ///
  EBPhysIBC* create() const
  {
    InflowOutflowAdvectBC* retval =  new InflowOutflowAdvectBC(m_flowDir, m_inflowVel, m_velComp);
    return static_cast<EBPhysIBC*>(retval);
  }


protected:

  int  m_velComp;
  int  m_flowDir;
  Real m_inflowVel;

private:
  //weak construction is bad but I will allow copy construction and assignment
  InflowOutflowAdvectBCFactory()
  {
    MayDay::Error("Invalid operator");
  }

};

#include "NamespaceFooter.H"

#endif
