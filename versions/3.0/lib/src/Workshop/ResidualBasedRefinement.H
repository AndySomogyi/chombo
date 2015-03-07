#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _RESIDUALBASEDREFINEMENT_H_
#define _RESIDUALBASEDREFINEMENT_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include "MayDay.H"

#include "NamespaceHeader.H"

class ResidualBasedRefinement: public RefinementCriterion
{
public:
  // constructors
  ResidualBasedRefinement(const Real & a_tolerance = REFINEMENTTHRESHOLD);

  ResidualBasedRefinement(const ResidualBasedRefinement & a_residualBasedRefinement);

  // destructor
  virtual ~ResidualBasedRefinement();

  virtual bool doRefine(const int & a_dim,
                        const Real & a_dxRatio);

  virtual bool doRefine(const Vector< Vector<Real> > & a_residual);

protected:
  int m_numOfRefinement;
  Real m_tolerance;
};

#include "NamespaceFooter.H"

#endif

