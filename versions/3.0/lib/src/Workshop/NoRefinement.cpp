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

#include "NoRefinement.H"

#include "NamespaceHeader.H"

NoRefinement::NoRefinement()
{
}

NoRefinement::~NoRefinement()
{
}

bool NoRefinement::doRefine(Vector<int>                 & a_refineInDir,
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

#include "NamespaceFooter.H"
