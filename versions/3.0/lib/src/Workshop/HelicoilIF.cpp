#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "HelicoilIF.H"

#include "NamespaceHeader.H"

HelicoilIF::HelicoilIF(const Real& a_helixR,
                       const Real& a_helixPitch,
                       const Real& a_circleR,
                       const bool& a_inside)
{
  // Remember the parameters
  m_helixR     = a_helixR;
  m_helixPitch = a_helixPitch;
  m_circleR    = a_circleR;

  // Save the square of the radius
  m_circleR2   = m_circleR * m_circleR;

  // Save inside flag
  m_inside = a_inside;
}

HelicoilIF::HelicoilIF(const HelicoilIF& a_inputIF)
{
  // Remember the parameters
  m_helixR     = a_inputIF.m_helixR;
  m_helixPitch = a_inputIF.m_helixPitch;
  m_circleR    = a_inputIF.m_circleR;

  // Save the square of the radius
  m_circleR2   = m_circleR * m_circleR;

  // Save inside flag
  m_inside = a_inputIF.m_inside;
}

HelicoilIF::~HelicoilIF()
{
}

Real HelicoilIF::value(const RealVect& a_point) const
{
  Real retval;

  Real theta;
  theta = atan2(a_point[1],a_point[0]);

  RealVect helixPt;

#if CH_SPACEDIM == 2
  helixPt[0] = m_helixR * cos(theta);
  helixPt[1] = m_helixR * sin(theta);
#elif CH_SPACEDIM == 3
  Real n = floor((a_point[2] - m_helixPitch * theta) / (2*M_PI * m_helixPitch));

  Real dz1 = fabs(a_point[2] - m_helixPitch * (theta + 2*M_PI *  n   ));
  Real dz2 = fabs(a_point[2] - m_helixPitch * (theta + 2*M_PI * (n+1)));

  Real Theta;

  if (dz1 < dz2)
  {
    Theta = theta + 2*M_PI*n;
  }
  else
  {
    Theta = theta + 2*M_PI*(n+1);
  }

  helixPt[0] = m_helixR * cos(Theta);
  helixPt[1] = m_helixR * sin(Theta);
  helixPt[2] = m_helixPitch * Theta;
#else
  MayDay::Abort("need higher dim in HelicoilIF\n");
#endif

  retval = 0.0;

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    Real delta = a_point[idir] - helixPt[idir];
    retval += delta*delta;
  }

  retval -= m_circleR2;

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

BaseIF* HelicoilIF::newImplicitFunction() const
{
  HelicoilIF* helicoilPtr;

  helicoilPtr = new HelicoilIF(m_helixR,m_helixPitch,m_circleR,m_inside);

  return static_cast<BaseIF*>(helicoilPtr);
}

#include "NamespaceFooter.H"
