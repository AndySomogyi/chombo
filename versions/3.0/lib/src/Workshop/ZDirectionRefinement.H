#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _ZDIRECTIONREFINEMENT_H_
#define _ZDIRECTIONREFINEMENT_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include "MayDay.H"

#include "NamespaceHeader.H"

class ZDirectionRefinement: public RefinementCriterion
{
public:
  // constructors
  ZDirectionRefinement(const bool & a_doNoRefinement = false);

  ZDirectionRefinement(const ZDirectionRefinement & a_zDirectionRefinement);

  // destructor
  virtual ~ZDirectionRefinement();

  virtual bool doRefine(const int & a_dim,
                        const Real & a_dxRatio);

  virtual bool doRefine(const Vector< Vector<Real> > & a_residual);

protected:
  int m_numOfRefinement;
  bool m_doNoRefinement;
};

#include "NamespaceFooter.H"

#endif

