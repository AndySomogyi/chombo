#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ComplementIF.H"

#include "NamespaceHeader.H"

ComplementIF::ComplementIF(const BaseIF& a_impFunc,
                           const bool&   a_complement)
{
  m_impFunc = a_impFunc.newImplicitFunction();
  m_complement = a_complement;
}

ComplementIF::ComplementIF(const ComplementIF& a_inputIF,
                           const bool&         a_complement)
{
  m_impFunc = a_inputIF.m_impFunc->newImplicitFunction();
  m_complement = a_complement;
}

ComplementIF::~ComplementIF()
{
  delete m_impFunc;
}

void ComplementIF::GetParams(bool& a_complement) const
{
  // Copy parameter information over
  a_complement = m_complement;
}

void ComplementIF::SetParams(const bool& a_complement)
{
  // Set parameter information
  m_complement = a_complement;
}

Real ComplementIF::value(const RealVect& a_point) const
{
  Real retval;

  // Implicit function value
  retval = m_impFunc->value(a_point);

  // Return the negative if complement is turned on (true)
  if (m_complement)
  {
    retval = -retval;
  }

  return retval;
}

Real ComplementIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  RealVect point;
  for (int idir = 0 ; idir < GLOBALDIM ; idir++)
    {
      point[idir] = a_point[idir];
    }

  return value(point);
}


IndexTM<Real,GLOBALDIM> ComplementIF::grad(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  IndexTM<Real,GLOBALDIM> retval;

  // Gradient of implicit function
  retval = m_impFunc->grad(a_point);

  // Return the negative if complement is turned on (true)
  if (m_complement)
  {
    for (int idir = 0; idir < GLOBALDIM; idir++)
    {
      retval[idir] = -retval[idir];
    }
  }

  return retval;
}

IndexTM<Real,GLOBALDIM> ComplementIF::normal(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  IndexTM<Real,GLOBALDIM> retval;

  // Normal of implicit function
  retval = m_impFunc->normal(a_point);

  // Return the negative if complement is turned on (true)
  if (m_complement)
  {
    for (int idir = 0; idir < GLOBALDIM; idir++)
    {
      retval[idir] = -retval[idir];
    }
  }

  return retval;
}

Vector<IndexTM <Real,GLOBALDIM> > ComplementIF::gradGrad(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  Vector<IndexTM<Real,GLOBALDIM> > retval;

  // Gradient of the gradient of the implicit function
  retval = m_impFunc->gradGrad(a_point);

  // Return the negative if complement is turned on (true)
  if (m_complement)
  {
    for (int idir = 0; idir < GLOBALDIM; idir++)
    {
      for (int jdir = 0; jdir < GLOBALDIM; jdir++)
      {
        retval[idir][jdir] = -retval[idir][jdir];
      }
    }
  }

  return retval;
}

Vector<IndexTM <Real,GLOBALDIM> > ComplementIF::gradNormal(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  Vector<IndexTM<Real,GLOBALDIM> > retval;

  // Gradient of the normal of the implicit function
  retval = m_impFunc->gradNormal(a_point);

  // Return the negative if complement is turned on (true)
  if (m_complement)
  {
    for (int idir = 0; idir < GLOBALDIM; idir++)
    {
      for (int jdir = 0; jdir < GLOBALDIM; jdir++)
      {
        retval[idir][jdir] = -retval[idir][jdir];
      }
    }
  }

  return retval;
}
GeometryService::InOut ComplementIF::InsideOutside(const RealVect& a_low, const RealVect& a_high) const
{
  GeometryService::InOut r = m_impFunc->InsideOutside(a_low, a_high);
  if(!m_complement) return r;
  if(r == GeometryService::Regular) return GeometryService::Covered;
  if(r == GeometryService::Covered) return GeometryService::Regular;
  return GeometryService::Irregular;
}



BaseIF* ComplementIF::newImplicitFunction() const
{
  ComplementIF* complementPtr = new ComplementIF(*m_impFunc,m_complement);

  return static_cast<BaseIF*>(complementPtr);
}

#include "NamespaceFooter.H"

