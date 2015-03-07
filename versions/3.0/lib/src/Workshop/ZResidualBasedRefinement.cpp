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
#include "ZResidualBasedRefinement.H"

#include "NamespaceHeader.H"

ZResidualBasedRefinement::ZResidualBasedRefinement(const bool & a_refineInZ,
                                                   const bool & a_refineWithResidual,
                                                   const Real & a_tolerance)
{
  m_numOfRefinement = 0;
  m_tolerance = a_tolerance;
  // m_refineInZ = a_refineInZ;
  m_refineWithResidual = a_refineWithResidual;
}

ZResidualBasedRefinement::ZResidualBasedRefinement(const ZResidualBasedRefinement & a_ZresidualBasedRefinement)
{
}

ZResidualBasedRefinement::~ZResidualBasedRefinement()
{
}

bool ZResidualBasedRefinement::doRefine(Vector<int>                  & a_refineInDir,
                                        const int                    & a_dim,
                                        const Vector<Real>           & a_dx,
                                        const Vector< Vector<Real> > & a_residual)
{
  CH_assert(a_refineInDir.size() == a_dim);
  CH_assert(a_dx         .size() == a_dim);

  for (int idir = 0; idir < a_dim; idir++)
  {
    a_refineInDir[idir] = 0;
  }

#if 0
  if (m_refineInZ)
  {
    if (a_dim == 3)
    {
      if (a_dxRatio > 1.0)
      {
        // we are in 3D and dx is still much smaller than dz => refine in Z
        m_numOfRefinement++;
      }
      else
      {
        // we are in 3D and dx is the same order as dz => stop the refinement
        // in Z, authorize the refinement with residuals
        m_refineInZ = false;
      }
    }
  }

  if (m_refineInZ)
  {
    a_refineInDir[GLOBALDIM] = 1;
  }

  return m_refineInZ;
#else
  return false;
#endif
}

#if 0
bool ZResidualBasedRefinement::doRefine(const Vector< Vector<Real> > & a_residual)
{
  int size = a_residual.size();
  CH_assert(size!=0);
  int numRes = 3;
  if (m_refineInZ || m_refineWithResidual == false ||a_residual[0][0] == LARGEREALVAL)
    {
      return false;
    }
  else
    {
      for (int iDegree = 0 ; iDegree < size ; iDegree++)
        {
          for (int iRes = 0 ; iRes < numRes ; iRes++)
            {
              if (a_residual[iDegree][iRes] > m_tolerance)
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
#endif

#include "NamespaceFooter.H"

