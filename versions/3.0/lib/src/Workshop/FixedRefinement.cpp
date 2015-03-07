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

#include "FixedRefinement.H"

#include "NamespaceHeader.H"

FixedRefinement::FixedRefinement()
{
}

FixedRefinement::FixedRefinement(const int & a_counter)
{
  setCounter(a_counter);
}

FixedRefinement::~FixedRefinement()
{
}

bool FixedRefinement::doRefine(Vector<int>                 & a_refineInDir,
                               const int                   & a_dim,
                               const Vector<Real>          & a_dx,
                               const Vector<Vector<Real> > & a_residual)
{
  CH_assert(a_refineInDir.size() == a_dim);
  CH_assert(a_dx         .size() == a_dim);

  bool retval = false;

  if (m_counter > 0)
  {
    retval = true;

    for (int idir = 0; idir < a_dim; idir++)
    {
      a_refineInDir[idir] = 1;
    }

    m_counter--;
  }

  return retval;
}

void FixedRefinement::setCounter(const int & a_counter)
{
  if (a_counter < 0)
  {
    MayDay::Abort("FixedRefinement::setCounter - counter must be >= 0");
  }

  m_counter = a_counter;
}

int FixedRefinement::getCounter()
{
  return m_counter;
}

#include "NamespaceFooter.H"
