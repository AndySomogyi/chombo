#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using std::endl;
#include "RefinementCriterion.H"
#include "FacesRefinement.H"

#include "NamespaceHeader.H"

FacesRefinement::FacesRefinement(const bool & a_refinedOnce,
                                 const bool & a_refinedTwice,
                                 const bool & a_refinedThree,
                                 const bool & a_refinedFour)
{
  m_refinedOnce = a_refinedOnce;
  m_refinedTwice = a_refinedTwice;
  m_refinedThree = a_refinedThree;
  m_refinedFour = a_refinedFour;
  // refinement of faces is isotropic
  // m_refineInZ = false;
}

FacesRefinement::FacesRefinement(const FacesRefinement & a_facesRefinement)
{
  m_refinedOnce = a_facesRefinement.m_refinedOnce;
  m_refinedTwice = a_facesRefinement.m_refinedTwice;
  m_refinedThree = a_facesRefinement.m_refinedThree;
  m_refinedFour = a_facesRefinement.m_refinedFour;
  // m_refineInZ = a_facesRefinement.m_refineInZ;
}

// destructor
FacesRefinement::~FacesRefinement()
{
}

bool FacesRefinement::doRefine(Vector<int>                 & a_refineInDir,
                               const int                   & a_dim,
                               const Vector<Real>          & a_dx,
                               const Vector<Vector<Real> > & a_residual)
{
  CH_assert(a_refineInDir.size() == a_dim);
  CH_assert(a_dx         .size() == a_dim);

  for (int idir = 0; idir < a_dim; idir++)
  {
    a_refineInDir[idir] = 0;
  }

  return false;
}

bool FacesRefinement::doRefine(const int & a_dim,
                               const Real & a_dxRatio)
{
  if (a_dim == 2 && m_refinedOnce == false)
    {
      setRefinedOnce();
      return true;
    }
  else if (a_dim == 2 && m_refinedTwice == false)
    {
      setRefinedTwice();
      return true;
    }
  else if (a_dim == 2 && m_refinedThree == false)
    {
      setRefinedThree();
      return true;
    }
  else if (a_dim == 2 && m_refinedFour == false)
    {
      setRefinedFour();
      return true;
    }
  return false;
}

bool FacesRefinement::doRefine(const Vector< Vector<Real> > & a_residual)
{
  return false;
}

void FacesRefinement::setRefinedOnce()
{
  m_refinedOnce = true;
}

void FacesRefinement::setRefinedTwice()
{
  m_refinedTwice = true;
}

void FacesRefinement::setRefinedThree()
{
  m_refinedThree = true;
}

void FacesRefinement::setRefinedFour()
{
  m_refinedFour = true;
}

#include "NamespaceFooter.H"
