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
using std::cout;

#include "RefinementCriterion.H"
#include "ZDirectionRefinement.H"

#include "NamespaceHeader.H"

ZDirectionRefinement::ZDirectionRefinement(const bool & a_doNoRefinement)
{
  m_numOfRefinement = 0;
  m_doNoRefinement = a_doNoRefinement;
  // m_refineInZ = true;
}

ZDirectionRefinement::ZDirectionRefinement(const ZDirectionRefinement & a_zDirectionRefinement)
{
  //TODO
}

ZDirectionRefinement::~ZDirectionRefinement()
{
}

bool ZDirectionRefinement::doRefine(const int & a_dim,
                                    const Real & a_dxRatio)
{
  if (m_doNoRefinement)
    {
      return false;
    }
  //TODO: change 1.0 to input parameter
  else if (a_dim == 3 && a_dxRatio > 1.0)
    {
      m_numOfRefinement++;
      //cout<<"dim = "<<a_dim<<",number of refinement = "<<m_numOfRefinement<<endl;
      return true;
    }
  //cout<<"dim = "<<a_dim<<",number of refinement = "<<m_numOfRefinement<<endl;
  return false;
}

bool ZDirectionRefinement::doRefine(const Vector< Vector<Real> > & a_residual)
{
  return false;
}

#include "NamespaceFooter.H"
