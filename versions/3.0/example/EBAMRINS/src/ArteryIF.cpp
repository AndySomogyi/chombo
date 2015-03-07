#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ArteryIF.H"

#include "NamespaceHeader.H"

ArteryIF::ArteryIF(const int&      a_type,
                   const Real&     a_radius,
                   const Real&     a_amplitude,
                   const Real&     a_minX,
                   const Real&     a_maxX,
                   const RealVect& a_center,
                   const bool&     a_inside)
{
  if (a_type != 1 && a_type != 2)
  {
    MayDay::Error("ArteryIF::ArteryIF: Artery type must be 1 or 2");
  }

  // Remember the parameters
  m_type      = a_type;
  m_radius    = a_radius;
  m_amplitude = a_amplitude;
  m_minX      = a_minX;
  m_maxX      = a_maxX;
  m_center    = a_center;
  m_inside    = a_inside;
}

ArteryIF::ArteryIF(const ArteryIF& a_inputIF)
{
  // Remember the parameters
  m_type      = a_inputIF.m_type;
  m_radius    = a_inputIF.m_radius;
  m_amplitude = a_inputIF.m_amplitude;
  m_minX      = a_inputIF.m_minX;
  m_maxX      = a_inputIF.m_maxX;
  m_center    = a_inputIF.m_center;
  m_inside    = a_inputIF.m_inside;
}

ArteryIF::~ArteryIF()
{
}

Real ArteryIF::value(const RealVect& a_point) const
{
  Real retval;

  Real x = a_point[0];

  // Set x to 0.0 at the beginning and 1.0 at the end of half the portion of
  // the artery length in excess of 1.0.  This creates a constant diameter
  // inlet and outlet.
  if ((x - m_minX) < ((m_maxX - m_minX) - 1.0)/2.0)
  {
    x = 0.0;
  }
  else if ((x - m_minX) < ((m_maxX - m_minX) + 1.0)/2.0)
  {
    x = (x - m_minX) - ((m_maxX - m_minX) - 1.0)/2.0;
  }
  else
  {
    x = 1.0;
  }

  Real func = 0.5 * (cos(2*M_PI*    x)
                   - cos(2*M_PI*3.5*x));

  Real perturb = m_amplitude * func;

  if (m_type == 1)
  {
    Real radius = 0.0;

    for (int idir = 1; idir < SpaceDim; idir++)
    {
      radius += (a_point[idir]-0.5)*(a_point[idir]-0.5);
    }

    radius = sqrt(radius);

    retval = radius - (m_radius + perturb);
  }

  if (m_type == 2)
  {
    Real radius = 0.0;

    if (SpaceDim == 2)
    {
      radius += ((a_point[1]-m_center[1]) - perturb*sin(2*M_PI*x))
              * ((a_point[1]-m_center[1]) - perturb*sin(2*M_PI*x));
    }

    if (SpaceDim == 3)
    {
      radius += ((a_point[1]-m_center[1]) - perturb*cos(2*M_PI*x))
              * ((a_point[1]-m_center[1]) - perturb*cos(2*M_PI*x));

      radius += ((a_point[2]-m_center[2]) - perturb*sin(2*M_PI*x))
              * ((a_point[2]-m_center[2]) - perturb*sin(2*M_PI*x));
    }

    radius = sqrt(radius);

    retval = radius - (m_radius + perturb);
  }

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* ArteryIF::newImplicitFunction() const
{
  ArteryIF* arteryPtr = new ArteryIF(m_type,
                                     m_radius,
                                     m_amplitude,
                                     m_minX,
                                     m_maxX,
                                     m_center,
                                     m_inside);

  return static_cast<BaseIF*>(arteryPtr);
}

#include "NamespaceFooter.H"
