#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _FIXEDREFINEMENT_H_
#define _FIXEDREFINEMENT_H_

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
   This refinement criterion refines a fixed number of levels.
 */
class FixedRefinement: public RefinementCriterion
{
public:
  /// Null constructor
  FixedRefinement();

  /// Constructor which sets counter
  FixedRefinement(const int & a_counter);

  /// Destructor
  virtual ~FixedRefinement();

  /// Refine a fixed number of times
  /**
     Refine until the counter gets to zero.  If the counter is non-zero,
     return true, "a_refineInDir" one in all components, and decrement counter.
     If the counter is zero, return false and "a_refineInDir" zero in all
     components.
   */
  virtual bool doRefine(Vector<int>                 & a_refineInDir,
                        const int                   & a_dim,
                        const Vector<Real>          & a_dxRatio,
                        const Vector<Vector<Real> > & a_residual);

  /// Set the counter
  virtual void setCounter(const int & a_counter);

  /// Get the counter
  virtual int getCounter();

protected:
  int m_counter;
};

#include "NamespaceFooter.H"

#endif
