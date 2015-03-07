#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// NodePoissonBC.cpp

#include "Box.H"
#include "REAL.H"
#include "SPACE.H"
#include "Tuple.H"
#include "Vector.H"
#include "LoHiSide.H"
#include "NodePoissonBC.H"
#include "MayDay.H"
#include "testProb_F.H"

// ---------------------------------------------------------
NeumannNodeBC::NeumannNodeBC() : FaceNodeBC()
{
}


// ---------------------------------------------------------
NeumannNodeBC::NeumannNodeBC(int dir, Side::LoHiSide side)
  : FaceNodeBC(dir, side)
{
  m_inhomogeneous = false;
}


// ---------------------------------------------------------
NeumannNodeBC::NeumannNodeBC(int dir, Side::LoHiSide side, const Interval& a_comps)
  : FaceNodeBC(dir, side, a_comps)
{
  m_inhomogeneous = false;
}


// ---------------------------------------------------------
NeumannNodeBC::~NeumannNodeBC()
{
}


// ---------------------------------------------------------
FaceNodeBC*
NeumannNodeBC::new_boxBC() const
{
  NeumannNodeBC* newop = new NeumannNodeBC();
  if (newop == NULL)
    {
      MayDay::Error("Out of Memory in NeumannNodeBC::new_boxBC");
    }
  newop->setInhomogeneous(m_inhomogeneous);
  return static_cast<FaceNodeBC*>(newop);
}


// ---------------------------------------------------------
void
NeumannNodeBC::setInhomogeneous(bool a_inhomogeneous)
{
  m_inhomogeneous = a_inhomogeneous;
}


// ---------------------------------------------------------
void
NeumannNodeBC::fillBCValues(FArrayBox& a_neumfac,
                            FArrayBox& a_dircfac,
                            FArrayBox& a_inhmval,
                            Real a_dx,
                            const Box& a_domain) const
{
  ProblemDomain probdomain(a_domain);
  fillBCValues(a_neumfac, a_dircfac, a_inhmval, a_dx, probdomain);
}


// ---------------------------------------------------------
void
NeumannNodeBC::fillBCValues(FArrayBox& a_neumfac,
                            FArrayBox& a_dircfac,
                            FArrayBox& a_inhmval,
                            Real a_dx,
                            const ProblemDomain& a_domain) const
{
 CH_assert(!a_neumfac.box().isEmpty());
 CH_assert(!a_dircfac.box().isEmpty());
 CH_assert(!a_inhmval.box().isEmpty());
  // 1*dphi/dn + 0*phi = 0
  a_neumfac.setVal(1.0);
  a_dircfac.setVal(0.0);
  a_inhmval.setVal(0.0);
  // if (m_inhomogeneous) do something else
  // a_inhmval.copy(m_inhmval);
}

// ---------------------------------------------------------

// ---------------------------------------------------------
DirichletNodeBC::DirichletNodeBC() : FaceNodeBC()
{
}


// ---------------------------------------------------------
DirichletNodeBC::DirichletNodeBC(int dir, Side::LoHiSide side)
  : FaceNodeBC(dir, side)
{
  m_inhomogeneous = false;
}


// ---------------------------------------------------------
DirichletNodeBC::DirichletNodeBC(int dir, Side::LoHiSide side, const Interval& a_comps)
  : FaceNodeBC(dir, side, a_comps)
{
  m_inhomogeneous = false;
}


// ---------------------------------------------------------
DirichletNodeBC::~DirichletNodeBC()
{
}


// ---------------------------------------------------------
FaceNodeBC*
DirichletNodeBC::new_boxBC() const
{
  DirichletNodeBC* newop = new DirichletNodeBC();
  if (newop == NULL)
    {
      MayDay::Error("Out of Memory in DirichletNodeBC::new_boxBC");
    }
  newop->setInhomogeneous(m_inhomogeneous);
  return static_cast<FaceNodeBC*>(newop);
}


// ---------------------------------------------------------
void
DirichletNodeBC::setInhomogeneous(bool a_inhomogeneous)
{
  m_inhomogeneous = a_inhomogeneous;
}


// ---------------------------------------------------------
void
DirichletNodeBC::fillBCValues(FArrayBox& a_neumfac,
                              FArrayBox& a_dircfac,
                              FArrayBox& a_inhmval,
                              Real a_dx,
                              const Box& a_domain) const
{
  ProblemDomain probdomain(a_domain);
  fillBCValues(a_neumfac, a_dircfac, a_inhmval, a_dx, probdomain);
}


// ---------------------------------------------------------
void
DirichletNodeBC::fillBCValues(FArrayBox& a_neumfac,
                              FArrayBox& a_dircfac,
                              FArrayBox& a_inhmval,
                              Real a_dx,
                              const ProblemDomain& a_domain) const
{
 CH_assert(!a_neumfac.box().isEmpty());
 CH_assert(!a_dircfac.box().isEmpty());
 CH_assert(!a_inhmval.box().isEmpty());
  // 0*dphi/dn + 1*phi = 0
  a_neumfac.setVal(0.0);
  a_dircfac.setVal(1.0);
  if (m_inhomogeneous)
    {
      FORT_GETPHINODEQUADRATIC(CHF_FRA(a_inhmval),
                               CHF_BOX(a_inhmval.box()),
                               CHF_CONST_REAL(a_dx));
    }
  else
    {
      a_inhmval.setVal(0.0);
    }
}
