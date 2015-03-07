#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "Box.H"
#include "FArrayBox.H"
#include "REAL.H"
#include "SPACE.H"
#include "Tuple.H"
#include "LoHiSide.H"
#include "MayDay.H"

#include "scalarDirichletBC.H"

void
scalarDirichletBC::fillBCValues(FArrayBox& a_neumfac,
                                FArrayBox&       a_dircfac,
                                FArrayBox&       a_inhmval,
                                Real             a_dx,
                                const Box&       a_domain) const
{
  ProblemDomain physdomain(a_domain);
  fillBCValues(a_neumfac, a_dircfac, a_inhmval, a_dx, physdomain);
}

void
scalarDirichletBC::fillBCValues(FArrayBox&     a_neumfac,
                          FArrayBox&           a_dircfac,
                          FArrayBox&           a_inhmval,
                          Real                 a_dx,
                          const ProblemDomain& a_domain) const
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());
  a_neumfac.setVal(0.0);
  a_dircfac.setVal(1.0);
  a_inhmval.setVal(m_BCVal);
}

BoxGhostBC* 
scalarDirichletBC::new_boxghostbc() const
{
  scalarDirichletBC* newop = new scalarDirichletBC();
  if(newop == NULL)
    {
      MayDay::Error("Out of Memory in scalarDirichletBC::new_boxghostbc");
    }

  newop->setVal(m_BCVal);
  return static_cast<BoxGhostBC*>(newop);
}


scalarDirichletBC::scalarDirichletBC() : BoxGhostBC()
{
  m_BCVal = 1.0e8;
}

scalarDirichletBC::~scalarDirichletBC()
{
}

scalarDirichletBC::scalarDirichletBC(int            a_dir,
                                     Side::LoHiSide a_side)
    : BoxGhostBC(a_dir,a_side)
{
  m_BCVal = 1.0e8;
}

scalarDirichletBC::scalarDirichletBC(int             a_dir,
                                     Side::LoHiSide  a_side,
                                     const Interval& a_comps)
    : BoxGhostBC(a_dir,a_side,a_comps)
{
  m_BCVal = 1.0e8;

}


void
scalarDirichletBC::setVal(const Real a_bcVal) 
{
  m_BCVal = a_bcVal;
}


scalarInflowBC::scalarInflowBC() : BoxGhostBC()
{
}

scalarInflowBC::scalarInflowBC(int            a_dir, 
                               Side::LoHiSide a_side) 
  : BoxGhostBC(a_dir,a_side)
{
}

scalarInflowBC::scalarInflowBC(int             a_dir, 
                               Side::LoHiSide  a_side,
                               const Interval& a_comps) 
  : BoxGhostBC(a_dir,a_side,a_comps)
{
}

scalarInflowBC::~scalarInflowBC()
{
}

void
scalarInflowBC::setScalarType(int a_scalType)
{
  m_scalType = a_scalType;
}

int
scalarInflowBC::scalarType() const
{
  return m_scalType;
}

void 
scalarInflowBC::isHomogeneous(bool a_isHomogeneous)
{
  m_isHomogeneous = a_isHomogeneous;
}

bool
scalarInflowBC::isHomogeneous() const
{
  return m_isHomogeneous;
}

void
scalarInflowBC::fillBCValues(FArrayBox& a_neumfac,
                             FArrayBox& a_dircfac,
                             FArrayBox& a_inhmval,
                             Real       a_dx,
                             const Box& a_domain) const
{
  ProblemDomain physdomain(a_domain);
  fillBCValues(a_neumfac, a_dircfac, a_inhmval, a_dx, physdomain);
}

void
scalarInflowBC::fillBCValues(FArrayBox&           a_neumfac,
                             FArrayBox&           a_dircfac,
                             FArrayBox&           a_inhmval,
                             Real                 a_dx,
                             const ProblemDomain& a_domain) const
  {
  // this will need to be reworked a little.
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());
  a_neumfac.setVal(0.0);
  a_dircfac.setVal(1.0);
  a_inhmval.setVal(m_scalType);
}

BoxGhostBC*
scalarInflowBC::new_boxghostbc() const
{

  scalarInflowBC* newop = new scalarInflowBC();
  if(newop == NULL)
    {
      MayDay::Error("Out of Memory in scalarInflowBC::new_boxghostbc");
    }

  newop->setScalarType(m_scalType);
  return static_cast<BoxGhostBC*>(newop);
}
