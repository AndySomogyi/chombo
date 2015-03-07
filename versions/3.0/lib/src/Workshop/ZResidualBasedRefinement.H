#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _ZRESIDUALBASEDREFINEMENT_H_
#define _ZRESIDUALBASEDREFINEMENT_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include "MayDay.H"

#include "NamespaceHeader.H"

class ZResidualBasedRefinement: public RefinementCriterion
{
public:
  //constructors
  ZResidualBasedRefinement(const bool & a_refineInZ = false,
                           const bool & a_refineWithResidual = false,
                           const Real & a_tolerance = MACHINEPRECISION);

  ZResidualBasedRefinement(const ZResidualBasedRefinement & a_ZresidualBasedRefinement);

  //destructor
  virtual ~ZResidualBasedRefinement();

  virtual bool doRefine(Vector<int>                  & a_refineInDir,
                        const int                    & a_dim,
                        const Vector<Real>           & a_dx,
                        const Vector< Vector<Real> > & a_residual);

  int  m_numOfRefinement;

protected:
  Real m_tolerance;
  bool m_refineWithResidual;
};

#include "NamespaceFooter.H"

#endif

