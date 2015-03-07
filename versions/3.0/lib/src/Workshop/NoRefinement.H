#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _NOREFINEMENT_H_
#define _NOREFINEMENT_H_

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include "MayDay.H"
#include "REAL.H"
#include "Vector.H"

#include "Notation.H"
#include "RefinementCriterion.H"

#include "NamespaceHeader.H"

///
/**
   This refinement criterion always false.  Used to stop subdivision.
 */
class NoRefinement: public RefinementCriterion
{
public:
  /// Null constructor
  NoRefinement();

  /// Destructor
  virtual ~NoRefinement();

  /// Never refine - return false and "a_refineInDir" zero in all components.
  /**
     Never refine - return false and "a_refineInDir" zero in all components.
   */
  virtual bool doRefine(Vector<int>                 & a_refineInDir,
                        const int                   & a_dim,
                        const Vector<Real>          & a_dxRatio,
                        const Vector<Vector<Real> > & a_residual);

protected:
};

#include "NamespaceFooter.H"

#endif
