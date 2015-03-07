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
#include "ResidualBasedRefinement.H"

#include "NamespaceHeader.H"

ResidualBasedRefinement::ResidualBasedRefinement(const Real & a_tolerance)
{
  m_numOfRefinement = 0;
  m_tolerance = a_tolerance;
  // residual based refinement is isotropic
  // m_refineInZ = false;
}

ResidualBasedRefinement::ResidualBasedRefinement(const ResidualBasedRefinement & a_residualBasedRefinement)
{
  //TODO
}

ResidualBasedRefinement::~ResidualBasedRefinement()
{
}

bool ResidualBasedRefinement::doRefine(const int & a_dim,
                                       const Real & a_dxRatio)
{
  return false;
}

bool ResidualBasedRefinement::doRefine(const Vector< Vector<Real> > & a_residual)
{
  int size = a_residual.size();
  CH_assert(size!=0);
  int numRes = 3;
  if (a_residual[0][0] == LARGEREALVAL)
    {
      //for all in or all out
      return false;
    }
  else
    {
      for (int iDegree = 0 ; iDegree < size ; iDegree++)
        {
          for (int iRes = 0 ; iRes < numRes ; iRes++)
            {
              if (a_residual[iDegree][iRes] > m_tolerance)
                //TODO: think about:only need to use max norm
                {
                  m_numOfRefinement++;
                  //cout<<"dim = "<<a_dim<<",number of refinement = "<<m_numOfRefinement<<endl;
                  return true;
                }
            }
        }
    }
  //cout<<"dim = "<<a_dim<<",number of refinement = "<<m_numOfRefinement<<endl;
  return false;
}

#include "NamespaceFooter.H"
